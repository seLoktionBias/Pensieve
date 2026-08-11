#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

echo "[1/11] CLI help forms"
for flag in -h -help --help; do python bin/pensieve "$flag" >/dev/null; done
python bin/pensieve --help -long > /tmp/pensieve_long_help.$$
grep -q "Pensieve v3.36 - full manual" /tmp/pensieve_long_help.$$
grep -q "diagnostics < alignment < events < asr < integrate < plot" /tmp/pensieve_long_help.$$
rm -f /tmp/pensieve_long_help.$$

echo "[2/11] Python syntax"
python -m py_compile scripts/*.py bin/pensieve tests/*.py

echo "[3/11] Shell syntax"
for f in scripts/*.sh install.sh tests/*.sh; do bash -n "$f"; done

echo "[4/11] Biological/orchestration synthetic tests"
python tests/backend_consistency_test.py

echo "[5/11] Reference-free + safe codeml exit tests"
python tests/paml_exit_and_reference_free_test.py

echo "[6/11] codeml control file"
grep -Eq '^clock[[:space:]]*=[[:space:]]*0' templates/dummy_codon_asr.ctl
grep -Eq '^fix_blength[[:space:]]*=[[:space:]]*0' templates/dummy_codon_asr.ctl
grep -Eq '^cleandata[[:space:]]*=[[:space:]]*0' templates/dummy_codon_asr.ctl
grep -Eq '^RateAncestor[[:space:]]*=[[:space:]]*1' templates/dummy_codon_asr.ctl

echo "[7/11] No legacy authoritative event/integration contracts"
if grep -R -n -E 'candidate_indel_frameshift_events|paml_indelmap_asr_combined' scripts/run_one_gene_00_to_04.sh scripts/04_ancestral_orf_walk.py 2>/dev/null; then
  echo "Legacy authoritative integration reference remains" >&2; exit 1
fi

echo "[8/11] No hard-coded user paths"
if grep -R -n -E '/home/[a-z0-9_]+/|/Users/' bin scripts templates 2>/dev/null; then
  echo "Hard-coded user path found" >&2; exit 1
fi

echo "[9/11] Plot execution smoke test (runs when Rscript is available)"
bash tests/plot_smoke_test.sh

echo "[10/11] No nested Slurm child-job dependency"
if grep -R -n -E 'sbatch.*codeml|child_job_registry|CHILD_JOB_REGISTRY' scripts/02_run_asr_backends.sh scripts/run_one_gene_00_to_04.sh 2>/dev/null; then
  echo "Nested Slurm child-job orchestration remains" >&2; exit 1
fi

echo "[11/11] Portable installation/runtime environment handling"
bash tests/install_portability_test.sh

rm -rf scripts/__pycache__ bin/__pycache__ tests/__pycache__
echo "Smoke test passed."
