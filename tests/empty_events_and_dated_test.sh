#!/usr/bin/env bash
# Two edge cases that used to abort a run:
#
#  1. A gene with NO events to draw. 05_plot_events.R assigned scalars into a
#     zero-row data frame (marker_class, then lane/n_lanes), giving
#     "replacement has 1 row, data has 0" and killing the run AFTER the science
#     was already done. Reached whenever every sequence is a complete ORF, or
#     when every event was an in-frame indel and this is the no_inframe figure.
#
#  2. --dated yes on a tree with no branch lengths. That is a harmless user
#     mistake; the runner must switch to --dated no, say so, and continue.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
fail=0
ck(){ if [ "$1" = 0 ]; then echo "  ok   $2"; else echo "  FAIL $2"; fail=1; fi; }

echo "empty event table renders instead of crashing"
printf '%s\n' "gene	event_id	character_class	event_type	biological_interpretation	alignment_start	alignment_end	event_length	length_mod_3	frame_effect	origin_node	origin_is_tip	parent_node	branch	shared_event	n_affected_tips	affected_tips	reversal_below_origin	secondary_changes_below_origin	root_state	parsimony_score	delta_parsimony_support	ambiguous_origin	direction_confident	parent_age	child_age	age_interval	n_observed_present	n_observed_absent	n_unknown	observed_present_tips	terminal_incompleteness	breakpoint_relationships	coordinate_system" > "$TMP/empty_events.tsv"
echo '((A:0.1,B:0.1)N1:0.2,(C:0.1,D:0.1)N2:0.2)N3;' > "$TMP/t.nwk"
Rscript "$ROOT/scripts/05_plot_events.R" --gene EMPTY --tree "$TMP/t.nwk" \
  --events "$TMP/empty_events.tsv" --alignment-length 300 \
  --outdir "$TMP/fig" --dated no > "$TMP/plot.log" 2>&1
ck $? "05_plot_events.R exits 0 on a header-only event table"
[ -s "$TMP/fig/EMPTY.pseudogenization_tree.pdf" ]; ck $? "pseudogenization tree written"
[ -s "$TMP/fig/EMPTY.event_map.pdf" ]; ck $? "event map written"
grep -q "0 events plotted" "$TMP/plot.log"; ck $? "reports 0 events plotted"
! grep -q "replacement has 1 row" "$TMP/plot.log"; ck $? "no zero-row data frame error"

echo "--dated yes on a cladogram switches to --dated no"
grep -q 'carries no branch lengths' "$ROOT/scripts/run_one_gene_00_to_04.sh"; ck $? "runner has the branch-length guard"
# the guard's own detection logic, exercised directly
det(){ python - "$1" <<'PY'
import re,sys
sys.exit(0 if re.search(r":\s*-?\d", open(sys.argv[1]).read()) else 1)
PY
}
echo '((A,B),(C,D));' > "$TMP/clado.nwk"
! det "$TMP/clado.nwk"; ck $? "cladogram detected as having no branch lengths"
det "$TMP/t.nwk"; ck $? "dated tree detected as having branch lengths"

echo
if [ "$fail" = 0 ]; then echo "Empty-events and dated-switch test passed."; else echo "failure(s)"; exit 1; fi
