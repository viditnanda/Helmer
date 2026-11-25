#!/bin/sh
# Run every example and stop at the first failure (used by CI). The helmer module must be importable.
HERE="$(cd "$(dirname "$0")" && pwd)"
LOG=${TMPDIR:-/tmp}/helmer_examples.log
for f in "$HERE"/*.py; do
    printf '%s ... ' "$(basename "$f")"
    python3 "$f" > "$LOG" 2>&1 || { echo FAILED; cat "$LOG"; exit 1; }
    echo ok
done
printf 'helmer command line ... '
python3 -m helmer "y^2*z - x^2" > "$LOG" 2>&1 && grep -q "dim 0" "$LOG" || { echo FAILED; cat "$LOG"; exit 1; }
echo ok
echo "all examples ran"
