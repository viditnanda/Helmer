#pragma once
#include <flint/fmpz.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif
typedef struct HelmerF4Trace HelmerF4Trace;
#ifdef HELMER_F4_TRACE
HelmerF4Trace* helmer_trace_create(const int32_t* lens, const int32_t* exps,
    const fmpz* const* coeffs, size_t terms, int vars, int gens);
int64_t helmer_trace_run(HelmerF4Trace*, uint32_t prime, int32_t* size,
    int32_t** lengths, int32_t** exps, void** coeffs);
void helmer_trace_free(HelmerF4Trace*);
#endif
#ifdef __cplusplus
}
#endif
