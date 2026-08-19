#!/usr/bin/env bash
# Pensieve v3.32 ASR backends.
#
# The parent Pensieve process owns scheduling.  In --mode slurm the *whole*
# pipeline runs inside one Slurm allocation; codeml is an ordinary child process
# in that allocation.  Nested Slurm submission is intentionally disabled because
# it made failure propagation and resume logic unreliable in earlier versions.
set -euo pipefail

GENE=""; WORKDIR="$PWD"; THREADS=4; TIME_HOURS=24
PAML_MODE="local"; ENV_NAME=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --gene) GENE="$2"; shift 2;;
    --workdir) WORKDIR="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --time-hours|--time) TIME_HOURS="$2"; shift 2;;
    --paml-mode) PAML_MODE="$2"; shift 2;;
    --env-name) ENV_NAME="$2"; shift 2;;
    --indelmap|--indelmap-dir) shift 2;;
    --child-job-registry|--slurm-account|--slurm-partition) shift 2;;
    --alignment|--dated|--run_from_step|--run-from-step|--run_up_to|--run-up-to|--dat-dir) shift 2;;
    *) echo "[ERROR] Unknown option: $1" >&2; exit 1;;
  esac
done

[[ -n "$GENE" ]] || { echo "[ERROR] --gene is required" >&2; exit 1; }
if [[ "$PAML_MODE" == "slurm" ]]; then
  echo "[WARN] --paml-mode slurm is deprecated; running codeml inside the current allocation." >&2
  PAML_MODE="local"
fi
[[ "$PAML_MODE" == "local" || "$PAML_MODE" == "same" ]] || { echo "[ERROR] --paml-mode must be local/same/slurm" >&2; exit 1; }

cd "$WORKDIR"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
OUTDIR="results_02/$GENE"
mkdir -p "$OUTDIR"

require_file(){ [[ -s "$1" ]] || { echo "[ERROR] Missing or empty required file ($2): $1" >&2; exit 1; }; }

SEQFILE="$OUTDIR/02_${GENE}.codon_for_paml.phy"
TREEFILE="$OUTDIR/02_${GENE}.tree_for_asr.nwk"
require_file "$SEQFILE" step02_codon_for_paml
require_file "$TREEFILE" step02_tree_for_asr

PAML_DIR="$OUTDIR/paml_codon_asr"
mkdir -p "$PAML_DIR"
cp -f "$SEQFILE" "$PAML_DIR/02_${GENE}.codon_for_paml.phy"
cp -f "$TREEFILE" "$PAML_DIR/tree.nwk"
sed -e "s|ssssss|02_${GENE}.codon_for_paml.phy|" \
    -e "s|tttttt|tree.nwk|" \
    -e "s|oooooo|codon_asr.out|" \
    "$PACKAGE_DIR/templates/dummy_codon_asr.ctl" > "$PAML_DIR/codon_asr.ctl"

# Remove products from any earlier attempt before launching codeml.  This is
# essential for safe recovery: a stale rst must never be mistaken for a valid
# marginal reconstruction from the current invocation.
rm -f "$PAML_DIR/rst" "$PAML_DIR/codon_asr.out" "$PAML_DIR/lnf" \
      "$PAML_DIR/2NG.dN" "$PAML_DIR/2NG.dS" "$PAML_DIR/2NG.t" \
      "$PAML_DIR/4fold.nuc" "$PAML_DIR/rub" "$PAML_DIR/codeml.log" \
      "$OUTDIR/02_${GENE}.paml_marginal_validation.tsv" \
      "$OUTDIR/02_${GENE}.codeml_run_status.tsv"

echo "[$(date)] codeml: marginal codon ASR (clock=0, fix_blength=0, cleandata=0)"

# Pensieve requires only the marginal ASR ("(1) Marginal reconstruction of
# ancestral sequences" in rst), never PAML's optional joint reconstruction
# that follows it.  Joint reconstruction is far more memory-hungry (its own
# conP/space allocations routinely dwarf the marginal pass) and is a common
# source of OOM kills/crashes on large gene trees that Pensieve does not use
# and does not need.
#
# codeml is therefore launched in the background and polled: as soon as its
# own rst shows the literal "Joint reconstruction of ancestral sequences"
# marker, the preceding marginal section is guaranteed to already be flushed
# to disk (a stream is written in order, so that marker cannot appear on disk
# before everything written ahead of it), so codeml is stopped right there
# instead of being left to risk a crash/OOM during work Pensieve discards
# anyway. If codeml instead exits on its own (small genes commonly finish
# both passes before the first poll), it is left alone and handled exactly as
# before.
( cd "$PAML_DIR" && exec codeml codon_asr.ctl ) > "$PAML_DIR/codeml.log" 2>&1 &
CODEML_PID=$!

STOPPED_BEFORE_JOINT=0
POLL_SECONDS="${PENSIEVE_ASR_POLL_SECONDS:-10}"
while kill -0 "$CODEML_PID" 2>/dev/null; do
  if [[ -s "$PAML_DIR/rst" ]] && grep -q "Joint reconstruction of ancestral sequences" "$PAML_DIR/rst" 2>/dev/null; then
    echo "[$(date)] codeml has written the complete marginal ASR section and is entering joint reconstruction, which Pensieve does not use; stopping codeml now." >&2
    kill "$CODEML_PID" 2>/dev/null || true
    for _ in 1 2 3 4 5 6 7 8 9 10; do
      kill -0 "$CODEML_PID" 2>/dev/null || break
      sleep 1
    done
    kill -9 "$CODEML_PID" 2>/dev/null || true
    STOPPED_BEFORE_JOINT=1
    break
  fi
  sleep "$POLL_SECONDS"
done
set +e
wait "$CODEML_PID"
CODEML_STATUS=$?
set -e
if [[ "$STOPPED_BEFORE_JOINT" -eq 1 ]]; then
  echo "[$(date)] codeml stopped by Pensieve after completing marginal ASR (exit status $CODEML_STATUS reflects that stop, not a codeml failure)." >&2
fi

# Validate the actual dependency rather than treating codeml's process status as
# sufficient.  Some PAML builds/datasets complete and flush marginal ASR, then
# return non-zero while entering the later joint-reconstruction routine; that
# case, and the deliberate early stop above, are both handled the same way.
set +e
python "$SCRIPT_DIR/02_validate_paml_marginal.py" \
  --rst "$PAML_DIR/rst" \
  --phylip "$SEQFILE" \
  --status-out "$OUTDIR/02_${GENE}.paml_marginal_validation.tsv" \
  > "$PAML_DIR/marginal_validation.log" 2>&1
MARGINAL_STATUS=$?
set -e

if [[ "$MARGINAL_STATUS" -ne 0 ]]; then
  echo "[ERROR] codeml did not produce a complete, validated marginal ancestral reconstruction for $GENE." >&2
  echo "[ERROR] codeml exit status: $CODEML_STATUS" >&2
  cat "$PAML_DIR/marginal_validation.log" >&2 || true
  tail -80 "$PAML_DIR/codeml.log" >&2 || true
  printf 'gene\tcodeml_exit_code\tmarginal_asr_validated\tpipeline_action\n%s\t%s\tFalse\tFAIL\n' \
    "$GENE" "$CODEML_STATUS" > "$OUTDIR/02_${GENE}.codeml_run_status.tsv"
  exit 1
fi

if [[ "$STOPPED_BEFORE_JOINT" -eq 1 ]]; then
  echo "[WARN] codeml exit status $CODEML_STATUS reflects Pensieve stopping it deliberately after a complete marginal ASR was validated." >&2
  echo "[WARN] Pensieve does not require PAML joint reconstruction; continuing with the validated marginal ASR." >&2
  printf 'gene\tcodeml_exit_code\tmarginal_asr_validated\tpipeline_action\n%s\t%s\tTrue\tSTOPPED_BEFORE_JOINT_RECONSTRUCTION_WITH_VALIDATED_MARGINAL_ASR\n' \
    "$GENE" "$CODEML_STATUS" > "$OUTDIR/02_${GENE}.codeml_run_status.tsv"
elif [[ "$CODEML_STATUS" -ne 0 ]]; then
  echo "[WARN] codeml returned exit status $CODEML_STATUS for $GENE AFTER a complete marginal ASR was validated." >&2
  echo "[WARN] Pensieve does not require PAML joint reconstruction; continuing with the validated marginal ASR." >&2
  if grep -q "Joint reconstruction" "$PAML_DIR/codeml.log" 2>/dev/null; then
    echo "[WARN] codeml log reached the joint-reconstruction phase before the non-zero exit." >&2
  fi
  printf 'gene\tcodeml_exit_code\tmarginal_asr_validated\tpipeline_action\n%s\t%s\tTrue\tCONTINUE_WITH_VALIDATED_MARGINAL_ASR\n' \
    "$GENE" "$CODEML_STATUS" > "$OUTDIR/02_${GENE}.codeml_run_status.tsv"
else
  printf 'gene\tcodeml_exit_code\tmarginal_asr_validated\tpipeline_action\n%s\t0\tTrue\tCONTINUE\n' \
    "$GENE" > "$OUTDIR/02_${GENE}.codeml_run_status.tsv"
fi

require_file "$PAML_DIR/rst" paml_rst
require_file "$OUTDIR/02_${GENE}.paml_marginal_validation.tsv" paml_marginal_validation

echo "[$(date)] ASR backends complete for $GENE"
