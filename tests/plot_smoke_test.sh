#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if ! command -v Rscript >/dev/null 2>&1; then
  echo "SKIP plot smoke test: Rscript is not installed in this packaging environment."
  exit 0
fi

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
cat > "$tmp/tree.nwk" <<'EOF'
((A:1,B:1)NodeAB:1,(C:1,D:1)NodeCD:1)UserRoot;
EOF
cat > "$tmp/events.tsv" <<'EOF'
gene	event_id	character_class	event_type	biological_interpretation	alignment_start	alignment_end	event_length	length_mod_3	frame_effect	origin_node	origin_is_tip	parent_node	branch	shared_event	n_affected_tips	affected_tips	reversal_below_origin	secondary_changes_below_origin	root_state	parsimony_score	delta_parsimony_support	ambiguous_origin	direction_confident	parent_age	child_age	age_interval	n_observed_present	n_observed_absent	n_unknown	observed_present_tips	breakpoint_relationships	coordinate_system
T	E1	indel	gap_change	deletion	10	12	3	0	in_frame	NodeAB	False	UserRoot	UserRoot->NodeAB	True	2	A,B	False	0	residue	1	1	False	True	2	1	1-2	2	2	0	A,B	shared_core	canonical_codon_alignment
T	E2	indel	gap_change	deletion	40	40	1	1	frameshift	A	True	NodeAB	NodeAB->A	False	1	A	False	0	residue	1	1	False	True	1	0	0-1	1	3	0	A	terminal	canonical_codon_alignment
T	E3	stop_mask	stop_change	ambiguous_stop_change	70	72	3	0	nonsense	NodeCD	False	UserRoot	UserRoot->NodeCD	True	2	C,D	False	0	ambiguous	1	0	True	False	2	1	1-2	2	2	0	C,D	NA	canonical_codon_alignment
EOF
cat > "$tmp/orf.tsv" <<'EOF'
gene	branch	parent_node	child_node	child_is_tip	parent_coding_status	child_coding_status	orf_transition	transition_evidence	uncertainty_reason	known_pseudogenic_history	confident_disabling_events_on_branch	ambiguous_disabling_events_on_branch	all_events_on_branch	n_events_on_branch
T	UserRoot->NodeAB	UserRoot	NodeAB	False	unavailable	disrupted	root_adjacent_disruption_first_loss_unresolved	child_sequence_disrupted_without_confident_catalogued_branch_event	no_distinct_clock0_biological_root_sequence	True	NA	NA	E1	1
T	NodeAB->A	NodeAB	A	True	disrupted	intact	apparent_orf_restoration	inherited_pseudogenic_history_plus_current_orf_intact	NA	True	E2	NA	E2	1
T	NodeAB->B	NodeAB	B	True	disrupted	disrupted	already_pseudogenic	inherited_pseudogenic_history	NA	True	NA	NA	NA	0
T	UserRoot->NodeCD	UserRoot	NodeCD	False	unavailable	uncertain	ambiguous_disabling_event	ambiguous_disabling_event	event_origin_or_direction_ambiguous	unknown	NA	E3	E3	1
T	NodeCD->C	NodeCD	C	True	uncertain	intact	history_uncertain_no_confirmed_event	no_confident_disabling_event	entering_pseudogenic_history_unresolved	unknown	NA	NA	NA	0
T	NodeCD->D	NodeCD	D	True	uncertain	intact	history_uncertain_no_confirmed_event	no_confident_disabling_event	entering_pseudogenic_history_unresolved	unknown	NA	NA	NA	0
EOF
mkdir -p "$tmp/out"
Rscript "$ROOT/scripts/05_plot_events.R" \
  --gene T --tree "$tmp/tree.nwk" --events "$tmp/events.tsv" \
  --orf-transitions "$tmp/orf.tsv" --alignment-length 120 --outdir "$tmp/out" --dated yes >/dev/null
for f in T.pseudogenization_tree.pdf T.pseudogenization_tree.png T.event_map.pdf T.event_map.png T.plotted_events.tsv; do
  [[ -s "$tmp/out/$f" ]] || { echo "Missing/empty plot output: $f" >&2; exit 1; }
done
# Also execute the undated branch path.
rm -f "$tmp/out"/T.*
Rscript "$ROOT/scripts/05_plot_events.R" \
  --gene T --tree "$tmp/tree.nwk" --events "$tmp/events.tsv" \
  --orf-transitions "$tmp/orf.tsv" --alignment-length 120 --outdir "$tmp/out" --dated no >/dev/null
[[ -s "$tmp/out/T.pseudogenization_tree.pdf" && -s "$tmp/out/T.event_map.pdf" ]] || {
  echo "Undated plot smoke output missing" >&2; exit 1;
}
echo "Plot smoke test passed (dated + undated)."
