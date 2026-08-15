#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROGRAM_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

GENE=""; FASTA=""; TREE=""; WORKDIR="$PWD"; THREADS=8; TIME_HOURS=96
ALIGNMENT="perform"; DATED="yes"; RUN_FROM_STEP="AUTO"; RUN_UP_TO="FULL"
CLEAN=0; SKIP_PLOT=0; PAML_MODE="local"; DAT_DIR="$PROGRAM_DIR/dat"
INDELMAP_DIR="$PROGRAM_DIR/external/indelMaP"; RUN_INDELMAP="no"; ENV_NAME="Pensieve"
TIE_BREAK="none"; BREAKPOINT_TOLERANCE="0"; ON_MISSING_ROOT_SEQUENCE="warn"
ORF_LOSS_COST="1.0"; ORF_RESTORATION_COST="2.0"; PSEUDOGENIC_BOUNDARY_CAP_VOTES="2.0"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --gene) GENE="$2"; shift 2;;
    --fasta) FASTA="$2"; shift 2;;
    --tree) TREE="$2"; shift 2;;
    --workdir) WORKDIR="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --time-hours|--time) TIME_HOURS="$2"; shift 2;;
    --alignment) ALIGNMENT="$2"; shift 2;;
    --dated) DATED="$2"; shift 2;;
    --run_from_step|--run-from-step) RUN_FROM_STEP="$2"; shift 2;;
    --run_up_to|--run-up-to) RUN_UP_TO="$2"; shift 2;;
    --clean) CLEAN=1; shift;;
    --skip-plot) SKIP_PLOT=1; shift;;
    --paml-mode) PAML_MODE="$2"; shift 2;;
    --dat-dir) DAT_DIR="$2"; shift 2;;
    --indelmap-dir) INDELMAP_DIR="$2"; shift 2;;
    --indelmap) RUN_INDELMAP="$2"; shift 2;;
    --env-name) ENV_NAME="$2"; shift 2;;
    --tie-break) TIE_BREAK="$2"; shift 2;;
    --breakpoint-tolerance) BREAKPOINT_TOLERANCE="$2"; shift 2;;
    --on-missing-root-sequence) ON_MISSING_ROOT_SEQUENCE="$2"; shift 2;;
    --orf-loss-cost) ORF_LOSS_COST="$2"; shift 2;;
    --orf-restoration-cost) ORF_RESTORATION_COST="$2"; shift 2;;
    --pseudogenic-boundary-cap-votes) PSEUDOGENIC_BOUNDARY_CAP_VOTES="$2"; shift 2;;
    *) echo "[ERROR] Unknown argument: $1" >&2; exit 1;;
  esac
done

[[ -n "$GENE" && -n "$FASTA" && -n "$TREE" ]] || { echo "[ERROR] Need --gene --fasta --tree" >&2; exit 1; }
[[ "$TIME_HOURS" =~ ^[1-9][0-9]*$ ]] || { echo "[ERROR] --time must be a positive integer" >&2; exit 1; }
case "$ALIGNMENT" in perform|defined) ;; *) echo "[ERROR] --alignment perform|defined" >&2; exit 1;; esac
case "$DATED" in yes|no) ;; *) echo "[ERROR] --dated yes|no" >&2; exit 1;; esac
case "$RUN_INDELMAP" in yes|no|auto) ;; *) echo "[ERROR] --indelmap yes|no|auto" >&2; exit 1;; esac
case "$TIE_BREAK" in none|ancestral|terminal) ;; *) echo "[ERROR] --tie-break none|ancestral|terminal" >&2; exit 1;; esac
case "$PAML_MODE" in
  local|same) ;;
  slurm) echo "[WARN] nested --paml-mode slurm is deprecated; running codeml in the current allocation" >&2; PAML_MODE="local";;
  *) echo "[ERROR] --paml-mode local|same|slurm" >&2; exit 1;;
esac

mkdir -p "$WORKDIR"; WORKDIR="$(cd "$WORKDIR" && pwd)"
FASTA="$(readlink -f "$FASTA")"; TREE="$(readlink -f "$TREE")"
cd "$WORKDIR"
mkdir -p logs

trap 'status=$?; if [[ $status -ne 0 ]]; then echo "[ERROR] Pensieve failed while processing ${GENE:-unknown} (status=$status)" >&2; fi; exit $status' EXIT
trap 'echo "[ERROR] Failed at line $LINENO while processing ${GENE:-unknown}" >&2' ERR

require_file(){ [[ -s "$1" ]] || { echo "[ERROR] Missing or empty required file ($2): $1" >&2; exit 1; }; }
copy_if_exists(){ local s="$1" d="$2"; [[ -s "$s" ]] && cp -f "$s" "$d" || true; }

normalise_step(){
  local x="$1" kind="$2"
  case "$x" in
    AUTO|auto) [[ "$kind" == from ]] && echo diagnostics || { [[ "$SKIP_PLOT" -eq 1 ]] && echo integrate || echo plot; };;
    FULL|full) [[ "$SKIP_PLOT" -eq 1 ]] && echo integrate || echo plot;;
    diagnostics|indel_discovery|indel-discovery|macse|MACSE) echo diagnostics;;
    alignment|MSA|msa) echo alignment;;
    events|binary|indel_events|event_reconstruction) echo events;;
    asr|ASR|omega|paml|PAML) echo asr;;
    integrate|integration) echo integrate;;
    plot|R|r) echo plot;;
    *) echo "[ERROR] Unknown step '$x'. Use diagnostics, alignment, events, asr, integrate, plot." >&2; exit 1;;
  esac
}
step_index(){ case "$1" in diagnostics) echo 1;; alignment) echo 2;; events) echo 3;; asr) echo 4;; integrate) echo 5;; plot) echo 6;; esac; }
RUN_FROM_STEP="$(normalise_step "$RUN_FROM_STEP" from)"; RUN_UP_TO="$(normalise_step "$RUN_UP_TO" to)"
RUN_FROM_IDX="$(step_index "$RUN_FROM_STEP")"; RUN_UP_TO_IDX="$(step_index "$RUN_UP_TO")"
[[ "$RUN_FROM_IDX" -le "$RUN_UP_TO_IDX" ]] || { echo "[ERROR] run_from is later than run_up_to" >&2; exit 1; }
step_in_range(){ local i; i="$(step_index "$1")"; [[ "$i" -ge "$RUN_FROM_IDX" && "$i" -le "$RUN_UP_TO_IDX" ]]; }
prereq_needed(){ local i; i="$(step_index "$1")"; [[ "$i" -lt "$RUN_FROM_IDX" && "$i" -le "$RUN_UP_TO_IDX" ]]; }

clean_from(){
  [[ "$CLEAN" -eq 1 ]] || return 0
  case "$RUN_FROM_STEP" in
    diagnostics) rm -rf "results_00/$GENE" "results_01/$GENE" "results_02/$GENE" "results_03/$GENE" "final_results/$GENE";;
    alignment) rm -rf "results_02/$GENE" "results_03/$GENE" "final_results/$GENE";;
    events) rm -rf "results_03/$GENE" "final_results/$GENE"; rm -rf "results_02/$GENE/paml_codon_asr" "results_02/$GENE/indelmap_asr";;
    asr) rm -rf "results_02/$GENE/paml_codon_asr" "results_02/$GENE/indelmap_asr"; find "results_03/$GENE" -maxdepth 1 -type f ! -name "03_${GENE}.alignment_*" ! -name "03_${GENE}.frame_arithmetic_by_node.tsv" ! -name "03_${GENE}.pensieve_labelled_dated_tree.nwk" ! -name "03_${GENE}.internode_label_map.tsv" -delete 2>/dev/null || true; rm -rf "final_results/$GENE";;
    integrate) find "results_03/$GENE" -maxdepth 1 -type f \( -name "03_${GENE}.paml*" -o -name "03_${GENE}.ancestral*" -o -name "03_${GENE}.internode_label_crosswalk.tsv" -o -name "03_${GENE}.tree_topology*" -o -name "03_${GENE}.indelmap_concordance*" -o -name "04_${GENE}.*" \) -delete 2>/dev/null || true; rm -rf "final_results/$GENE";;
    plot) rm -rf "final_results/$GENE";;
  esac
}

require_diagnostics(){
  require_file "results_00/$GENE/00_${GENE}.common_species.fasta" step00_common_fasta
  require_file "results_00/$GENE/00_${GENE}.orf_status.tsv" step00_orf_status
  require_file "results_01/$GENE/01_${GENE}.macse_NT.fasta" step01_macse_nt
}
require_alignment(){
  require_file "results_02/$GENE/02_${GENE}.primary_codon_alignment_native.fasta" native_alignment
  require_file "results_02/$GENE/02_${GENE}.primary_codon_alignment.fasta" paml_safe_alignment
  require_file "results_02/$GENE/02_${GENE}.codon_for_paml.phy" paml_phylip
  require_file "results_02/$GENE/02_${GENE}.tree_for_asr.nwk" asr_tree
}
require_events(){
  require_file "results_03/$GENE/03_${GENE}.alignment_characters.tsv" alignment_characters
  require_file "results_03/$GENE/03_${GENE}.alignment_character_node_states.tsv" character_states
  require_file "results_03/$GENE/03_${GENE}.pensieve_labelled_dated_tree.nwk" pensieve_tree
}
require_asr(){
  local rst="results_02/$GENE/paml_codon_asr/rst"
  local validation="results_02/$GENE/02_${GENE}.paml_marginal_validation.tsv"
  local run_status="results_02/$GENE/02_${GENE}.codeml_run_status.tsv"
  require_file "$rst" paml_rst
  echo "[INFO] Validating PAML marginal ASR before integration/reuse for $GENE"
  python "$SCRIPT_DIR/02_validate_paml_marginal.py" \
    --rst "$rst" \
    --phylip "results_02/$GENE/02_${GENE}.codon_for_paml.phy" \
    --status-out "$validation"
  awk -F '\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="marginal_asr_valid") c=i; next} NR==2{exit !($c=="True")}' "$validation" || {
    echo "[ERROR] PAML marginal ASR is not validated for $GENE: $validation" >&2
    exit 1
  }
  if [[ ! -s "$run_status" ]]; then
    printf 'gene\tcodeml_exit_code\tmarginal_asr_validated\tpipeline_action\n%s\tUNKNOWN_PRIOR_RUN\tTrue\tREUSE_VALIDATED_EXISTING_MARGINAL_ASR\n' \
      "$GENE" > "$run_status"
  fi
}
require_integrated(){
  require_file "results_03/$GENE/03_${GENE}.paml_marginal_asr.fa" paml_marginal_asr
  require_file "results_03/$GENE/03_${GENE}.ancestral_integrated_alignment.fa" integrated_ancestral_alignment
  require_file "results_03/$GENE/04_${GENE}.ancestral_orf_walk.tsv" ancestral_orf_walk
  require_file "results_03/$GENE/04_${GENE}.orf_transitions_by_branch.tsv" orf_transitions
}

copy_final_outputs(){
  local important="final_results/$GENE/important_output" support="final_results/$GENE/supporting_files"
  mkdir -p "$important" "$support"
  copy_if_exists "results_03/$GENE/03_${GENE}.alignment_events.tsv" "$important/${GENE}.alignment_events.tsv"
  copy_if_exists "results_03/$GENE/03_${GENE}.alignment_characters.tsv" "$important/${GENE}.alignment_characters.tsv"
  copy_if_exists "results_03/$GENE/04_${GENE}.orf_transitions_by_branch.tsv" "$important/${GENE}.orf_transitions_by_branch.tsv"
  copy_if_exists "results_03/$GENE/04_${GENE}.ancestral_orf_walk.tsv" "$important/${GENE}.ancestral_orf_walk.tsv"
  copy_if_exists "results_03/$GENE/03_${GENE}.ancestral_integrated_alignment.fa" "$important/${GENE}.ancestral_integrated_alignment.fa"
  copy_if_exists "results_03/$GENE/03_${GENE}.phylogenetic_msa.fa" "$important/${GENE}.phylogenetic_msa.fasta"
  copy_if_exists "results_03/$GENE/03_${GENE}.ancestral_native_cds.fa" "$important/${GENE}.ancestral_native_cds.fa"
  copy_if_exists "results_03/$GENE/03_${GENE}.phylogenetic_sequence_order.tsv" "$important/${GENE}.phylogenetic_sequence_order.tsv"
  copy_if_exists "results_00/$GENE/00_${GENE}.phylogenetic_tip_order.tsv" "$important/${GENE}.phylogenetic_tip_order.tsv"
  copy_if_exists "results_03/$GENE/03_${GENE}.pensieve_labelled_dated_tree.nwk" "$important/${GENE}.pensieve_tree.nwk"
  copy_if_exists "results_03/$GENE/03_${GENE}.internode_label_crosswalk.tsv" "$important/${GENE}.internode_label_crosswalk.tsv"
  copy_if_exists "results_03/$GENE/03_${GENE}.indelmap_concordance.tsv" "$important/${GENE}.indelmap_concordance.tsv"
  copy_if_exists "results_02/$GENE/02_${GENE}.primary_codon_alignment_native.fasta" "$important/${GENE}.canonical_alignment.fasta"
  copy_if_exists "results_02/$GENE/02_${GENE}.primary_codon_alignment.fasta" "$support/${GENE}.paml_safe_alignment.fasta"
  copy_if_exists "results_02/$GENE/02_${GENE}.alignment_provenance.tsv" "$support/${GENE}.alignment_provenance.tsv"
  copy_if_exists "results_02/$GENE/02_${GENE}.masked_inframe_premature_stops_after_macse_correction.tsv" "$support/${GENE}.premature_stop_registry.tsv"
  copy_if_exists "results_03/$GENE/03_${GENE}.alignment_character_node_states.tsv" "$support/${GENE}.alignment_character_node_states.tsv"
  copy_if_exists "results_03/$GENE/03_${GENE}.frame_arithmetic_by_node.tsv" "$support/${GENE}.frame_arithmetic_by_node.tsv"
  copy_if_exists "results_03/$GENE/03_${GENE}.paml_marginal_asr.fa" "$support/${GENE}.paml_marginal_asr.fa"
  copy_if_exists "results_03/$GENE/03_${GENE}.paml_marginal_asr_all_declared_nodes.fa" "$support/${GENE}.paml_marginal_asr_all_declared_nodes.fa"
  copy_if_exists "results_02/$GENE/02_${GENE}.paml_marginal_validation.tsv" "$support/${GENE}.paml_marginal_validation.tsv"
  copy_if_exists "results_02/$GENE/02_${GENE}.codeml_run_status.tsv" "$support/${GENE}.codeml_run_status.tsv"
  copy_if_exists "results_03/$GENE/03_${GENE}.tree_topology_and_node_count_audit.tsv" "$support/${GENE}.tree_topology_and_node_count_audit.tsv"
  copy_if_exists "results_03/$GENE/04_${GENE}.orf_frame_crosscheck_mismatches.tsv" "$support/${GENE}.orf_frame_crosscheck_mismatches.tsv"
}

purge_legacy_reference_artifacts(){
  # v3.30 remains reference-free. Remove exact legacy files from older Pensieve
  # versions so reused work directories cannot retain misleading reference-based
  # diagnostics after a new run/resume.
  rm -f \
    "results_00/$GENE/00_${GENE}.reference_info.tsv" \
    "results_00/$GENE/00_${GENE}.reference_sequence.fasta" \
    "results_01/$GENE/01_${GENE}.macse_indels_relative_to_reference.tsv" \
    "final_results/$GENE/important_output/${GENE}.reference_info.tsv" \
    "final_results/$GENE/important_output/${GENE}.reference_sequence.fasta" \
    "final_results/$GENE/supporting_files/${GENE}.macse_indels_relative_to_reference.tsv" \
    2>/dev/null || true
}

require_file "$FASTA" input_fasta; require_file "$TREE" input_tree
clean_from
purge_legacy_reference_artifacts
mkdir -p "results_00/$GENE" "results_01/$GENE" "results_02/$GENE" "results_03/$GENE" "final_results/$GENE/important_output" "final_results/$GENE/supporting_files"

echo "============================================================"
echo "Pensieve $(cat "$PROGRAM_DIR/VERSION") | GENE=$GENE | alignment=$ALIGNMENT | dated=$DATED"
echo "range=$RUN_FROM_STEP->$RUN_UP_TO | indelmap=$RUN_INDELMAP | tie_break=$TIE_BREAK"
echo "workdir=$WORKDIR | start=$(date)"
echo "============================================================"

if step_in_range diagnostics; then
  echo "[$(date)] STEP diagnostics: prune/tree reconcile + raw ORF audit"
  python "$SCRIPT_DIR/00_prune_and_check_orf.py" --gene "$GENE" --fasta "$FASTA" --tree "$TREE" --outdir "results_00/$GENE"
  echo "[$(date)] STEP diagnostics: MACSE frameshift/STOP-aware alignment"
  python "$SCRIPT_DIR/01_run_macse_and_extract_events.py" --gene "$GENE" --step00-dir "results_00/$GENE" --outdir "results_01/$GENE"
elif prereq_needed diagnostics; then require_diagnostics; fi

if step_in_range alignment; then
  require_diagnostics
  echo "[$(date)] STEP alignment: build one canonical alignment + synchronized PAML-safe view"
  python "$SCRIPT_DIR/02_prepare_asr_inputs.py" --gene "$GENE" --results01-dir "results_01/$GENE" --results00-dir "results_00/$GENE" --outdir "results_02/$GENE" --alignment-mode "$ALIGNMENT"
elif prereq_needed alignment; then require_alignment; fi

if step_in_range events; then
  require_alignment
  echo "[$(date)] STEP events: breakpoint-coded phylogenetic event reconstruction"
  python "$SCRIPT_DIR/03_alignment_events.py" --gene "$GENE" \
    --alignment "results_02/$GENE/02_${GENE}.primary_codon_alignment_native.fasta" \
    --tree "results_02/$GENE/02_${GENE}.tree_for_asr.nwk" \
    --masked-stops "results_02/$GENE/02_${GENE}.masked_inframe_premature_stops_after_macse_correction.tsv" \
    --orf-status "results_00/$GENE/00_${GENE}.orf_status.tsv" \
    --orf-loss-cost "$ORF_LOSS_COST" --orf-restoration-cost "$ORF_RESTORATION_COST" \
    --pseudogenic-boundary-cap-votes "$PSEUDOGENIC_BOUNDARY_CAP_VOTES" \
    --outdir "results_03/$GENE" --dated "$DATED" --tie-break "$TIE_BREAK" \
    --breakpoint-tolerance "$BREAKPOINT_TOLERANCE" --coordinate-system "canonical_codon_alignment"
elif prereq_needed events; then require_events; fi

if step_in_range asr; then
  require_alignment; require_events
  echo "[$(date)] STEP asr: codeml substitution ASR + optional IndelMaP concordance backend"
  bash "$SCRIPT_DIR/02_run_asr_backends.sh" --gene "$GENE" --workdir "$WORKDIR" --threads "$THREADS" --time-hours "$TIME_HOURS" \
    --paml-mode local --indelmap "$RUN_INDELMAP" --indelmap-dir "$INDELMAP_DIR" --env-name "$ENV_NAME"
elif prereq_needed asr; then require_asr; fi

if step_in_range integrate; then
  require_alignment; require_events; require_asr
  echo "[$(date)] STEP integrate: parse/map PAML; compare optional IndelMaP without overwriting Pensieve states"
  python "$SCRIPT_DIR/03_integrate_asr_evidence.py" --gene "$GENE" --results02-dir "results_02/$GENE" --results00-dir "results_00/$GENE" --outdir "results_03/$GENE" --dated "$DATED" --on-missing-root-sequence "$ON_MISSING_ROOT_SEQUENCE"
  echo "[$(date)] STEP integrate: build native ancestors and sticky pseudogenization history"
  python "$SCRIPT_DIR/04_ancestral_orf_walk.py" --gene "$GENE" \
    --paml-tree "results_03/$GENE/03_${GENE}.paml_labeled_reporting_tree.nwk" \
    --pensieve-tree "results_03/$GENE/03_${GENE}.pensieve_labelled_dated_tree.nwk" \
    --ancestral-fasta "results_03/$GENE/03_${GENE}.paml_marginal_asr.fa" \
    --tip-alignment "results_02/$GENE/02_${GENE}.primary_codon_alignment_native.fasta" \
    --character-states "results_03/$GENE/03_${GENE}.alignment_character_node_states.tsv" \
    --characters "results_03/$GENE/03_${GENE}.alignment_characters.tsv" \
    --events "results_03/$GENE/03_${GENE}.alignment_events.tsv" \
    --frame-arithmetic "results_03/$GENE/03_${GENE}.frame_arithmetic_by_node.tsv" \
    --outdir "results_03/$GENE"
  require_integrated
  copy_final_outputs
elif prereq_needed integrate; then require_integrated; copy_final_outputs; fi

if step_in_range plot; then
  require_events; require_integrated
  if [[ "$SKIP_PLOT" -eq 0 ]]; then
    echo "[$(date)] STEP plot: history-aware tree and event map"
    ALN_LEN="$(awk 'NR==1{print $2}' "results_02/$GENE/02_${GENE}.codon_for_paml.phy")"
    Rscript "$SCRIPT_DIR/05_plot_events.R" --gene "$GENE" \
      --tree "results_03/$GENE/03_${GENE}.pensieve_labelled_dated_tree.nwk" \
      --events "results_03/$GENE/03_${GENE}.alignment_events.tsv" \
      --orf-transitions "results_03/$GENE/04_${GENE}.orf_transitions_by_branch.tsv" \
      --alignment-length "$ALN_LEN" \
      --outdir "final_results/$GENE/important_output" --dated "$DATED"
    require_file "final_results/$GENE/important_output/${GENE}.pseudogenization_tree.pdf" plot_pdf
  fi
fi

trap - EXIT ERR
echo "DONE $GENE $(date)"
