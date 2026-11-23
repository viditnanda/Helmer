// Adapter to the trace interface in msolve 0.7.5.
#include "f4_trace.h"
#include <msolve/neogb/libneogb.h>

struct HelmerF4Trace {
  bs_t* input;
  ht_t* hash;
  md_t* meta;
  int learned;
};

HelmerF4Trace* helmer_trace_create(const int32_t* lens, const int32_t* exps,
    const fmpz* const* coeffs, size_t terms, int vars, int gens) {
  HelmerF4Trace* trace = calloc(1, sizeof(*trace));
  // The characteristic-zero import expects pointers to numerator/denominator pairs.
  mpz_t* integers = malloc(terms * sizeof(*integers));
  mpz_t** pairs = malloc(2 * terms * sizeof(*pairs));
  if (!trace || !integers || !pairs) {
    free(pairs);
    free(integers);
    free(trace);
    return NULL;
  }
  mpz_t one;
  mpz_init_set_ui(one, 1);
  for (size_t i = 0; i < terms; ++i) {
    mpz_init(integers[i]);
    fmpz_get_mpz(integers[i], coeffs[i]);
    pairs[2*i] = &integers[i];
    pairs[2*i+1] = &one;
  }
  int ok = initialize_gba_input_data(&trace->input, &trace->hash, &trace->meta,
      lens, exps, pairs, 0, 0, 0, vars, gens, 0, 12, 1, 0, 0, 2, 0, 1, 0, 0, 0);
  for (size_t i = 0; i < terms; ++i) mpz_clear(integers[i]);
  mpz_clear(one);
  free(integers);
  free(pairs);
  if (ok != 1) {
    free(trace);
    return NULL;
  }
  return trace;
}

int64_t helmer_trace_run(HelmerF4Trace* trace, uint32_t prime, int32_t* size,
    int32_t** lengths, int32_t** exps, void** coeffs) {
  if (trace->learned && (!trace->meta->tr || !trace->meta->tr->ltd)) return 0;
  trace->meta->f4_qq_round = trace->learned ? 2 : 1;
  int32_t error = 0;
  bs_t* basis = core_gba(trace->input, trace->meta, &error, prime);
  if (error || !basis) return 0;
  md_t export_meta = *trace->meta;
  export_meta.fc = prime;
  export_meta.ff_bits = 32;
  md_t* meta = &export_meta;
  int64_t terms = export_results_from_gba(size, lengths, exps, coeffs, malloc,
                                           &basis, &basis->ht, &meta);
  // The learning basis shares its hash table with the input. Replay owns a copy.
  if (trace->learned) free_basis_and_only_local_hash_table_data(&basis);
  else free_basis_without_hash_table(&basis);
  trace->learned = 1;
  return terms;
}

void helmer_trace_free(HelmerF4Trace* trace) {
  if (!trace) return;
  if (trace->meta->tr) free_trace(&trace->meta->tr);
  free_shared_hash_data(trace->hash);
  free_basis(&trace->input);
  free(trace->meta);
  free(trace);
}
