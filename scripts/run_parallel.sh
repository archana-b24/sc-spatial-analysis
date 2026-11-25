#!/usr/bin/env bash
# Full pipeline with the independent branches run concurrently (DEVELOPMENT_LOG §3.12).
#
# run_all.py is strictly sequential; after `cnv`, four branches only read sc_annotated.h5ad (+ the DE table):
#   A: de → ccc            B: abundance → trajectory            C: deconvolve → spatial → ccc_spatial
# Running them side by side takes the wall time from ~sum to ~max of the branches (cell2location dominates).
# Each branch is an ordinary run_all.py call (same stamps, same per-step logs in logs/run_all/), with capped
# BLAS/torch threads so the branches do not oversubscribe the 16 cores.
#
#   bash scripts/run_parallel.sh [--force]          # --force re-runs every step (ignores completion stamps)
#
# Progress: logs/run_full.log (one line per step start/finish, prefixed by branch); per-step logs as usual.
set -uo pipefail
cd "$(dirname "$0")/.."
FORCE=${1:-}
LOG=logs/run_full.log
mkdir -p logs/run_all

run() {  # run <branch-label> <threads> <steps...>
  local label=$1 threads=$2; shift 2
  OMP_NUM_THREADS=$threads MKL_NUM_THREADS=$threads OPENBLAS_NUM_THREADS=$threads NUMBA_NUM_THREADS=$threads \
    pixi run python scripts/run_all.py --only "$@" $FORCE 2>&1 | sed -u "s/^/[$label] /" | tee -a "$LOG"
  return "${PIPESTATUS[0]}"
}

echo "=== run_parallel $(date -Iseconds) ${FORCE} ===" | tee -a "$LOG"
run upstream 16 fetch_visium spatial_prep select qc fm_embed integrate cnv || { echo "UPSTREAM FAILED" | tee -a "$LOG"; exit 1; }

run A 6 de ccc &
pa=$!
run B 4 abundance trajectory &
pb=$!
run C 6 deconvolve spatial ccc_spatial &
pc=$!
fail=0
wait $pa || fail=1
wait $pb || fail=1
wait $pc || fail=1
if [ $fail -ne 0 ]; then
  echo "BRANCH FAILED — see the ✗ lines above; fix and re-run the failed steps with run_all.py --only … --force" | tee -a "$LOG"
  exit 1
fi
run final 16 report export || exit 1
echo "PARALLEL RUN ALL DONE $(date -Iseconds)" | tee -a "$LOG"
