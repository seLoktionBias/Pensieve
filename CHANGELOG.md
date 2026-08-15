# Changelog

## v4.2 - a real STOP-overfitting bug fixed, and a correction to v4.1's boundary-capping scope

Prompted by user-reported results on `Bat_genes_from_Song` (`song_pensieve_v4_output`) plus a detailed independent review.

### Real bug: raw premature STOPs downstream of a frameshift could still be reported as independent nonsense mutations (CNGA3/Phyllostomus discolor)
- Reported directly: CNGA3's real, well-documented biology is a single early 1bp deletion in *Phyllostomus discolor* causing pseudogenization -- correcting that one lesion alone restores the reading frame for the rest of the gene, so no downstream premature stop in the shifted frame should be reported as its own independent event. `song_pensieve_v4_output`'s CNGA3 run (present in v4.0/v4.1, absent from the v3.37 baseline in `final_results`) over-reported these as independent STOP substitutions.
- Root cause, precisely characterized: `scripts/03_alignment_events.py`'s `stop_mask_characters()` had a fallback (`elif seg and seg == codon: PRESENT`) that resurrected ANY tip as a carrier of an already-founded STOP character purely because its own canonical-alignment segment textually matched the allele -- even when that tip's own raw-STOP occurrence had been explicitly excluded from founding a character because it was a known consequence of an earlier, unrelated frameshift (`frame_shifted_at_stop=True`). A frame-dependent downstream "STOP" only needs to coincide, by chance, with some OTHER unrelated lineage's real, independently-validated nonsense mutation at the same homologous codon position for this fallback to wrongly attach it to that character.
- `scripts/02_prepare_asr_inputs.py`'s `build_stop_registry()` also only gated `pseudogenizing_event_candidate`/`independent_stop_candidate` on frame-phase (`frame_shifted_at_stop`) and a raw mapping-succeeded check (`valid`); it never required the mapped columns to be **contiguous** (no MACSE `!`/gap wedged inside the 3bp span) nor re-verified that the actual **corrected homologous codon**, read back from the final canonical alignment, is still exactly TAA/TAG/TGA.
- Fixed with three independent, additive gates in `build_stop_registry()` -- `independent` (frame not shifted), `contiguous` + exactly-3bp span, and `codon_confirmed` (corrected homologous codon re-derived from the canonical alignment is still a STOP) -- all three must pass for `pseudogenizing_event_candidate`. New audit fields `corrected_homologous_codon`, `codon_confirmed_stop`; every rejection records a specific `reason`.
- Fixed the resurrection fallback in `stop_mask_characters()` more surgically than a blanket removal: a tip is now blocked from being resurrected as a carrier ONLY when its own occurrence at that exact span/allele was specifically classified `frame_shifted_at_stop=True` (an informative, specific piece of negative evidence). A tip with NO row at all, or one that was simply never independently classified for an unrelated/unregistered reason, still falls through to the original segment-text-match check -- this preserves a real, previously-fixed PDE6H regression (`test_shared_stop_not_missed_when_only_one_tip_is_registered`) where a genuinely shared mutation's second occurrence legitimately never got its own independent classification, and forcing it to ABSENT would have reintroduced that older bug.
- Verified end-to-end on real CNGA3/Phyllostomus discolor data (`CNGA3_Phyllostomus_stop_before_after.tsv`) and with new synthetic regression coverage.

### v4.1 correction: pseudogenic-component boundary capping was also overwriting the component-entry node's own state (ChatGPT review, ORF-aware-parsimony spec section 10)
- v4.1's `sankoff()` applied `clamp_boundary_evidence()` in the traceback/up-pass to decide a component-entry node's OWN C=ABSENT-vs-PRESENT state, not just the message it hands to its parent. This is stronger than the sampling-density fix requires and could overwrite real internal evidence AT the entry node itself.
- Fixed: capping is now applied ONLY in the down-pass, to what an entry node hands to its PARENT's own cost computation. The entry node's own reported state and confidence (`delta_parsimony_support`) always use full, unbounded evidence, exactly like every node strictly inside the component.
- Non-obvious mathematical finding made while fixing `tests/orf_aware_parsimony_test.py`'s `test_pseudogene_sampling_density_invariance` for this correction: ordinary binary Sankoff already bounds a SINGLE child's own pull on its DIRECT parent's cost difference at the transition cost (gain_cost/loss_cost), regardless of how large that child's own internal cost difference grows -- an inherent min-plus property, true with or without Stage C. This means the entry node's OWN reported confidence genuinely does (and per spec 10, should) grow unbounded with clade size -- confirmed with the same 2->320 tip sweep as v4.1, now correctly identical whether or not `--orf-status`/capping is supplied (0.95 to 96.35 either way) -- while its PARENT was already sampling-density invariant on its own for this topology, independent of Stage C. Stage C's `--pseudogenic-boundary-cap-votes` only makes an observable difference at the parent level when set BELOW that natural per-mutation transition-cost ceiling; the test now demonstrates this explicitly with a cap deliberately set under the ceiling (0.3 votes) alongside the loose/default one (2.0 votes). This also explains why v4.1's restoration-penalty sensitivity grid (ratios 2/4/8/16, all above the natural ceiling for indel gain/loss costs) gave byte-identical real GUCY2F output at every ratio: all four were already no-ops relative to uncapped for that gene's actual component topology.

### `--pseudogenic-boundary-cap-votes` defaults unchanged (2.0); documentation corrected
- README's "Pseudogenized clades do not get extra votes just for being densely sampled" section corrected to describe the v4.2 scope precisely (capping affects only the message crossing the component boundary upward, never the entry node's own state).

### VERSION bumped to 4.2.

## v4.1 - ORF-history-aware structural reconstruction (pseudogene-clade sampling bias), a real cladogram tip-alignment fix

### Real bug: a densely sampled pseudogenized clade could make its own pre-loss ancestor look implausibly damaged, purely from tip count
- Reported via a detailed review (real Bat_genes_from_Song/GUCY2F data): structural characters (indels, premature-STOP alleles) were reconstructed independently by binary Sankoff parsimony in `scripts/03_alignment_events.py`, with every tip fed in as if it were an independent observation. Tips descending from the same pseudogenization event are NOT independent evidence about the state of the ancestor that predates that event -- they share one correlated post-loss history. Confirmed directly against real `song_pensieve_output` GUCY2F data: `Node15->Molossus_molossus` (a genuinely fresh, independent loss off an intact parent) was reported `already_pseudogenic` (pale/inherited) instead of its own confident loss, and two branches (`Node5->Pteropus_alecto`, `Node5->Rousettus_aegyptiacus`) were reported `apparent_orf_restoration` (implying a prior loss that never happened) instead of plain `intact`.
- Root cause, precisely characterized with a real, minimal synthetic reproduction: for a single indel/STOP character, a pseudogenic clade's own entry-node Sankoff cost pair is NOT distorted merely by clade size when the clade is internally homogeneous (ordinary Sankoff already handles that correctly, at any size) -- the real sampling-density artifact only appears for a *heterogeneous* within-clade character (present in some, not all, descendants -- exactly how independently-accumulated passenger mutations look in real post-pseudogenization lineages). With the SAME 50/50 carrier:non-carrier composition held fixed and only clade size varied (2, 5, 10, ... 320 tips), the entry node's reported confidence (`delta_parsimony_support`) grew from 0.95 to 96.35 -- unbounded, purely from sample count, with the underlying biology genuinely unchanged.
- Fixed with a new three-stage, tree- and ORF-history-aware reconstruction layer in `scripts/03_alignment_events.py`:
  - **Stage A** (`pseudogenic_components()`): a provisional, tip-evidence-only Functional/Pseudogenic history is reconstructed once per gene, from `00_<gene>.orf_status.tsv`'s own `complete_orf` calls ONLY -- never from any structural character's own placement, to avoid using a potentially-biased indel/STOP reconstruction to define pseudogenization and then using that pseudogenization to "fix" the same reconstruction. Reuses the existing, already-tested `sankoff()` engine with an asymmetric F->P (`--orf-loss-cost`, default 1.0) vs P->F (`--orf-restoration-cost`, default 2.0) cost, exactly like the existing STOP/frameshift gain/loss asymmetry. The root is never forced: a genuine tie stays reported as ambiguous, exactly as `sankoff()` already does for every other character.
  - **Stage B**: every branch entering a contiguous pseudogenic (P) component -- confidently-F parent, confidently-P child -- is recorded as a component boundary. An uncertain entering history never invents a confident first-loss boundary.
  - **Stage C** (`clamp_boundary_evidence()`, threaded through `sankoff()` as `component_entry_idx`/`boundary_cap`): a pseudogenic component's own C=ABSENT-vs-PRESENT cost difference -- both what it hands to its parent's own down-pass sum, and its own up-pass decision -- is capped at `--pseudogenic-boundary-cap-votes` (default 2.0) "equivalent independent tip votes" worth of evidence, regardless of real descendant count. Everything strictly INSIDE a component keeps its full, unbounded internal resolution (post-loss mutations are still fully reconstructed and reported); only the component's influence on ancestors ABOVE its own boundary is bounded.
  - Simply re-running the character reconstruction on a "corrected" input does not fix this kind of bias by itself (unlike the v4.0 MACSE block-corruption fix) -- the bound has to be applied inside the Sankoff cost propagation itself, which is what Stage C does.
  - Two new outputs per gene: `03_<gene>.provisional_orf_history.tsv` (Stage A's own per-node F/P/ambiguous call) and `03_<gene>.pseudogenic_components.tsv` (Stage B's identified entry branches), for full auditability.
  - `scripts/run_one_gene_00_to_04.sh` threads `--orf-loss-cost`/`--orf-restoration-cost`/`--pseudogenic-boundary-cap-votes` through as overridable flags with the defaults above; `04_ancestral_orf_walk.py`'s existing sticky-pseudogenic-history semantics (fixed separately in v4.0) are unchanged and complementary -- Stage A's provisional history only bounds character-level evidence upstream of a component, it is not a second, competing authority on final per-node coding status.
- Restoration:loss cost ratio chosen from an explicit sensitivity grid (2, 4, 8, 16): real GUCY2F data (Bat_genes_from_Song, `--dated no`) produced byte-identical `alignment_events.tsv`/`provisional_orf_history.tsv`/`pseudogenic_components.tsv` at every ratio tested, and two synthetic stress cases (an all-pseudogenic tree; two independent losses vs one ancestral loss plus repeated restoration) gave identical, correct results throughout. Per the smallest-stable-ratio rule, the default is 2.0 (see `restoration_penalty_sensitivity.tsv`).
- Verified against real GUCY2F data end-to-end (fresh run vs the existing, untouched `song_pensieve_output` baseline): 6 of 36 branches changed classification, all in the expected direction -- `Node15->Molossus_molossus` now reports `pseudogenization` (a fresh, independent loss) instead of `already_pseudogenic`; `Node5->Pteropus_alecto`/`Node5->Rousettus_aegyptiacus` now report plain `intact` instead of `apparent_orf_restoration`; two branches previously miscategorized `already_pseudogenic` now correctly report `confirmed_disabling_event_first_loss_unresolved` (still a confident, solid-red disabling call, just no longer wrongly downgraded to "inherited"); one branch now correctly reports `sequence_state_uncertain` instead of a false confident pale call. The literal nucleotide content at the originally-reported problem region (columns 2790-2810) was itself already correct in both the old and new run -- the real defect was in the branch-level ORF-history classification, not that specific sequence window.
- New test suite `tests/orf_aware_parsimony_test.py`: `test_pseudogene_sampling_density_invariance` (the critical mandated test -- demonstrates the previously-unbounded confidence growth and confirms it is now capped), `test_multiple_losses_preferred_over_repeated_restoration_when_supported`, `test_root_not_forced_functional` (including a genuinely-ambiguous-root case), `test_apparent_restoration_retains_pseudogenic_history`, `test_all_functional_matches_neutral_behavior` (byte-identical output to the pre-v4.1 behavior when no pseudogenic component exists at all -- guards against accidental global weighting), `test_unknown_boundary_remains_ambiguous`.

### Real bug: with `--dated no`, tree tips did not line up at a common x position in the plots
- Reported directly: undated cladogram-style plots showed tips at different depths depending on topology instead of a clean, aligned right edge.
- Root cause, in `scripts/05_plot_events.R`: undated mode set every edge length to 1 and let `ape::plot.phylo` lay the tree out as an ordinary phylogram, which still spaces tips by their own edge COUNT from the root -- a deeper lineage lands further right than a shallow one, which is not a real cladogram. Verified directly with a synthetic unbalanced tree: tip x-coordinates came out `(3, 3, 2, 1)` the old way.
- Fixed by passing `use.edge.length = dated` to `plot()` instead -- ape's own actual cladogram layout, which ignores edge lengths entirely and aligns every tip at the same x position regardless of topology depth. Same synthetic tree: tip x-coordinates are now `(3, 3, 3, 3)`.

### Smoke test hardened against a false positive from Python bytecode cache
- `tests/smoke_test.sh`'s "no hard-coded user paths" check recursively grepped `bin scripts templates`, which also matched `__pycache__/*.pyc` files -- compiled bytecode always embeds the absolute source file path as a debugging artifact, so simply having run any test beforehand could trip this check on a real user's own machine. Now excludes `__pycache__`.

### VERSION bumped to 4.1.

## v4.0 - removed terminal-stop reconstruction and fixed two real ancestral-history bugs

### Terminal stop codons are stripped from every input sequence; Pensieve no longer reconstructs them
- Requested directly: after three successive fix attempts across v3.35-v3.37, terminal-stop reconstruction at ancestral nodes kept surfacing new real failure modes (PAML's codon model structurally cannot represent a STOP as a state, so every parsimony reconstruction of one from the tip alignment was working around that limitation, not solving it).
- `strip_terminal_stop()` (`scripts/00_prune_and_check_orf.py`) now removes a single trailing TAA/TAG/TGA from every input sequence's own real (gap-stripped) ending, unconditionally, before anything else in the pipeline sees it. There is nothing left for any downstream step to reconstruct, mis-scatter, or disagree about at the end of the CDS.
- Real bug found while verifying this against real PDE6C data: the first implementation still required the sequence's overall gapless length to be a multiple of 3 before stripping, so frameshifted/pseudogenized sequences (whose length is *not* a multiple of 3 precisely because they're frameshifted) kept an untouched trailing stop -- e.g. `Desmodus_rotundus`, `Diaemus_youngii` and 6 other real PDE6C species still ended `...TTGTAA` after the "fix." The mod-3 gate is removed entirely: if a sequence's own last 3 real characters spell a stop codon, they are stripped, regardless of what the rest of its length is.
- All ancestral terminal-stop machinery removed from `scripts/04_ancestral_orf_walk.py`: `own_terminal_codon()`, `fitch_reconstruct_codon()`, and the `terminal_fitch_columns` computation/application block.
- 3 obsolete tests removed from `tests/backend_consistency_test.py`; new test `test_terminal_stop_stripped_from_every_input_sequence()` covers both the simple case and the frameshifted-length case above.

### Real bug: MACSE could corrupt the conserved gene-start block even for species with a completely normal ATG start, because their real lesion is elsewhere in the gene
- Reported directly, with real PDE6C data pasted side by side: `Desmodus_rotundus`, `Diaemus_youngii` and other vampire-bat-clade species have a completely ordinary `ATG` start -- their real pseudogenizing lesion is elsewhere in the gene -- yet MACSE still inserted a spurious frame-correction marker right at their conserved 5' start, corrupting a region that never needed correcting, and (via the shared alignment columns) corrupting *every other* species' start codon too.
- Root cause: MACSE can place its frame-restoration marker anywhere its own global DP alignment finds cheapest, not necessarily near the lineage's actual lesion. Simply re-running MACSE on a "corrected" version of the full input does *not* fix this: MACSE sees the identical conserved region again and reliably makes the identical, wrong placement decision again (verified directly -- an earlier fix attempt that did exactly this produced byte-identical corrupted output).
- Fixed with a new preprocessing step in `scripts/01_run_macse_and_extract_events.py`: `detect_conserved_block()` walks codon-by-codon from position 1 across all complete-ORF raw sequences (no alignment needed -- complete ORFs have no internal frameshift by definition, so they're directly comparable), extending the block as long as modal amino-acid identity stays above a low floor (0.3 -- calibrated against real PDE6C data, where an *ordinary*, non-indel polymorphic site can legitimately dip to 0.29-0.60 between otherwise 90-100% conserved neighbours; a real indel boundary was never observed above that). For every species where MACSE modified this block, `resolve_conserved_block()` scores the species' own raw block content against MACSE's proposed edit -- nucleotide AND amino-acid identity, both -- against a complete-ORF consensus, and keeps MACSE's edit only if it's a strict improvement on at least one axis without being worse on the other.
- Critically, the block region is then *never shown to MACSE at all*: `main()` strips the decided block content off the front of every sequence and runs MACSE only on the remainder, splicing the (quality-gated) block back onto the front afterward. This is what actually prevents the corruption, rather than merely correcting it after the fact.
- Decisions are logged to the new `01_<gene>.conserved_block_correction.tsv`.
- Verified against real PDE6C data: with the block extended to the data-driven 30 codons, all 103 species' first 90nt are byte-identical to their own real raw sequence in the canonical MACSE alignment output -- including `Mormoops_megalophylla` (`ATT...`) and `Thyroptera_tricolor` (`GTG...`), whose own genuine non-ATG starts (simple substitutions, not indels) are preserved untouched rather than being replaced by a MACSE guess that scored 0% amino-acid identity against the conserved consensus.

### Real bug: two independently pseudogenized sister branches with a genuinely functional common ancestor rendered pale pink, not red -- and in two real genes, no branch ever rendered red at all
- Reported directly, comparing `GUCY2F`/`CNGB3` plots against others: no red-colored branches anywhere, only pale pink, even where two sister lineages clearly lost the gene independently.
- Root cause #1 (the dominant one, explaining both genes having *zero* `pseudogenization` branches): `04_ancestral_orf_walk.py`'s `root_disabling` check trusted *any* `stop_mask` character reconstructing `stop_present` at the tree root, without checking how well-supported that reconstruction actually was. Real GUCY2F/CNGB3 data: several such characters were observed in as few as 1-3 tips out of 100+, with `n_observed_absent == 0` everywhere else (every other tip simply lacks alignment coverage there, not confidently "absent"). With zero contradicting evidence anywhere in the tree, Fitch parsimony reconstructs "present" at every node including the root at zero cost -- a missing-data artifact, not real evidence the gene was broken before the sampled tree began. That single unsupported character flipped `root_disabling` True for the *entire* gene, making every subsequent independent loss anywhere in the tree render as inherited (`already_pseudogenic`, pale) instead of its own fresh loss (`pseudogenization`, solid red). Fixed: a stop_mask character now only sets `root_disabling` if `n_observed_absent > 0` -- i.e. at least one tip is confidently observed *without* the stop, giving the character real phylogenetic contrast.
- Root cause #2 (a rarer, secondary case): sticky `known_pseudogenic_history` (correctly kept True for a node whose current sequence is intact but whose more distant ancestor was pseudogenized -- `apparent_orf_restoration`) was blindly propagated to that node's *children* too. A child with its own brand-new, independent disabling event off such a node rendered `already_pseudogenic` (pale) even though the immediate parent is displayed grey/intact -- visually and semantically contradictory. Fixed: a parent whose own current sequence is confidently intact (`parent_coding_status == "intact"`) always resets history for its children, regardless of what its own more distant ancestor did.
- Two new regression tests in `tests/backend_consistency_test.py`: `test_sparse_root_stop_evidence_does_not_poison_whole_tree` and `test_intact_parent_resets_history_for_independently_disabled_children`.

### VERSION bumped to 4.0.

### Real bug: a single, contiguous indel could be reported as several separate events
- Reported directly, from real CNGA3/GUCA1C data: e.g. Nyctalus_aviator's real, contiguous 6bp deletion (columns 400-405) was reported as two events -- 1bp at 400, 5bp at 401-405 -- and CNGA3's Dasypterus_ega's real 24bp insertion likewise split into 3bp + 21bp pieces.
- Root cause: the breakpoint decomposition in `scripts/03_alignment_events.py` correctly splits a run into separate characters whenever some *other*, unrelated tip's own indel only overlaps part of the span (a real, necessary behaviour for cross-species comparison) -- but reporting each resulting character as its own event fragments what is, from a single branch's own point of view, one real indel. This could even manufacture fake frameshift signal: two fragments individually not a multiple of 3 in length, even when their true combined span is in-frame.
- Fixed with a post-processing merge: on one branch, two events of the same type (`deletion`, `insertion_or_restoration`, or `ambiguous_indel_change`) whose spans are exactly contiguous are now reported as one merged event. Runs before the frame-arithmetic step, so a previously-fragmented in-frame indel is now correctly counted as in-frame rather than as two frameshifts.
- Extended to shared (ancestral-origin) events too, guarded by requiring an *identical* affected-tip set between the two fragments -- real CNGA3 data (`Node95->Node94`) showed an ancestral 18bp+3bp deletion pair fragment exactly the same way a single tip's own event does, but a separate pair of genuinely different characters (6 tips vs. a disjoint 20) that happen to start at the same column must never merge just because they share a branch label.
- A second real bug found while verifying the first: a *nested* interior event on the same branch (e.g. a 3bp event strictly inside a larger 18bp one) sorted in between two genuinely adjacent outer pieces by start column, breaking the adjacency scan entirely. Fixed by excluding nested/contained intervals from the merge scan -- they stay their own separate event, and the outer pieces correctly merge around them.
- Five new regression tests in `tests/backend_consistency_test.py` covering: contiguous same-type merging, adjacent-different-type non-merging, contiguous shared-event merging with matching tips, same-branch-different-tips non-merging, and the nested-interior-event case.
- Verified against all eight genes: event counts dropped substantially wherever fragmentation was occurring (e.g. GUCA1C 411->290, CNGB3 2566->1661, GRK7 1123->805, GUCY2F 4331->2890), and a small number of pseudogenization calls that had depended on fake fragment-level frameshift signal correctly flipped to intact once merged.

### Real bug: multi-FASTA species order was the reverse of the rendered tree's top-to-bottom order
- Reported directly: comparing `ancestral_integrated_alignment.fa`/`phylogenetic_msa.fasta` against the tree figure side by side, the FASTA's species order read bottom-to-top instead of top-to-bottom.
- Root cause: `ape::plot.phylo` (`scripts/05_plot_events.R`, `direction="rightwards"`) draws the first tip in plain left-to-right preorder at the *bottom* of the figure and the last at the *top* -- every FASTA/order-table writer in the pipeline used that same plain preorder, so files read bottom-to-top next to a plot a person reads top-to-bottom.
- Fixed at all four order-construction sites (`scripts/00_prune_and_check_orf.py`, `01_run_macse_and_extract_events.py`, `02_prepare_asr_inputs.py`, `04_ancestral_orf_walk.py`) with a reverse-sibling preorder (still parent-before-children, so still topologically valid; only the left-right order of sibling subtrees at each node is flipped) that matches the plot's own top-to-bottom order exactly. Applies to `phylogenetic_msa.fasta`, `ancestral_integrated_alignment.fa`, `canonical_alignment.fasta`, and both tip/sequence order tables.
- Verified directly: the species drawn at the very top of every rendered tree this session (`Corynorhinus_mexicanus`) is now the first record in the exported FASTA; the one drawn at the very bottom (`Cynopterus_sphinx`) is last.

### Real bug: a branch could be labelled "pseudogenization" even though its own reconstructed sequence was a complete ORF
- Reported directly, from real PDE6C data: `Node59` was reported as pseudogenized purely because two indel events on its branch were each individually flagged frameshifting (length not a multiple of 3) -- despite Node59's own reconstructed sequence, gap-stripped and checked directly, being a genuinely complete, in-frame ORF, consistent with every descendant also being intact. `Node35` showed the identical pattern.
- Root cause, in `scripts/04_ancestral_orf_walk.py`'s branch-history walk: a branch was classified `pseudogenization` (or `confirmed_disabling_event_first_loss_unresolved`) as soon as *any* event was catalogued as "confident disabling" on it, without ever checking whether the resulting child sequence was actually disrupted. The two events' combined effect on the real sequence was a net multiple of 3 (not a single clean 3bp event, but together removing exactly 3bp) -- most likely reflecting the alignment being DNA-level rather than strictly codon-level, not a real frameshift ever occurring.
- Fixed by checking the child's own reconstructed `coding_status` first, with top priority, in both places this pattern occurred: if it is genuinely `intact`, the branch is reported `intact` regardless of what events are catalogued on it. The events are not silently discarded -- `confident_disabling_events_on_branch` still lists them, and `transition_evidence` records `catalogued_disabling_events_but_child_orf_intact` so the discrepancy stays auditable.
- New regression test (`tests/backend_consistency_test.py::test_confident_disabling_events_never_override_a_genuinely_intact_child_orf`) reproduces the exact PDE6C shape: two individually-frameshifting indel events whose combined effect on the real sequence is in-frame.
- Verified against real data: PDE6C's false-positive pseudogenization-branch count dropped from 16 to 11; `Node59`/`Node35` (and their equivalents in other genes) now correctly render grey (intact) rather than red.

### Terminology: "Pseudogenizing indel (frameshift)" renamed to "Frameshift indel"
- Requested directly: labelling an event "pseudogenizing" on a branch/tip whose own CDS is complete (per the fix above) reads as a direct contradiction. The legend entry for this marker class is now simply "Frameshift indel" -- the colour, and every other label, is unchanged.

## v3.36 - fixed a real verification-run mistake that produced uniformly grey plots, and hardened the plotting script against it

### Real bug: every branch and species name rendered plain grey/black, with no pseudogenization colour anywhere, on freshly regenerated plots
- Reported directly: after the v3.35 terminal-stop fix, regenerated `pseudogenization_tree.pdf`/`event_map.pdf` for all four genes showed no colour distinction at all -- every branch grey, every tip label black, contradicting real, known pseudogenization (e.g. GUCA1C's `Hipposideros` clade, CNGA3's several confirmed premature-stop lineages).
- Root cause: this was a mistake in the ad hoc command used to regenerate the plots for verification (not a bug reachable through the normal pipeline, `scripts/run_one_gene_00_to_04.sh`, which was and remains correct) -- it passed `--transitions`/`--orf-walk` instead of `scripts/05_plot_events.R`'s actual flag, `--orf-transitions`. `get_arg()` silently returns its default (`NA`, "no ORF history file") for any flag it doesn't recognise, so the script ran to completion with no error, no warning, and every branch defaulted to the `intact` (grey) category -- a uniformly grey, uniformly wrong plot with nothing in the output to indicate anything had gone wrong.
- Fixed the immediate output: all four genes' `pseudogenization_tree`/`event_map` PDFs and PNGs regenerated with the correct flag and now show the expected colour-coded history (confirmed disabling in saturated red, inherited in pale red, uncertain in amber, intact in grey), matching the underlying `orf_transitions_by_branch.tsv` data, which was correct throughout -- only the plot rendering was affected.
- Hardened `scripts/05_plot_events.R` itself against this exact class of mistake happening again, silently, to anyone (not just this session's manual invocation): it now collects every `--xxx`-shaped token in its arguments and errors out immediately, before doing any work, if any of them isn't one of its known flags -- `Unrecognised flag(s): --transitions. Known flags: --gene, --tree, --events, --orf-transitions, ...` -- instead of quietly falling back to a default that produces a wrong-but-plausible-looking plot. Also now errors immediately if `--orf-transitions` is given but the file doesn't exist, rather than silently treating a bad path the same as "no history data at all."
- Verified the hardening: re-running the exact bad command that caused this (`--transitions` instead of `--orf-transitions`) now fails fast with the error above instead of producing a bad plot; `tests/plot_smoke_test.sh` and the full `tests/smoke_test.sh` suite still pass with the correct flag.

### Real bug: PNG images silently went stale while PDFs were refreshed, and prefixed duplicate files accumulated in `final_results`
- Reported directly: delivered PDFs reflected the latest fix, but the paired PNGs in the same `final_results/<GENE>/important_output/` folder were still the previous, stale render -- `05_plot_events.R` always writes both formats together, but the ad hoc verification workflow this session had only copied the `.pdf` into `final_results`, never the `.png` written alongside it in the working directory.
- A second, related mistake in the same ad hoc workflow: rather than writing `05_plot_events.R`'s output directly to `final_results/<GENE>/important_output` (as `run_one_gene_00_to_04.sh` always does) and reusing the pipeline's own file-copy naming, results were manually `cp -f`'d in with their original `03_`/`04_`-prefixed working-directory filenames, alongside the already-present, correctly-named unprefixed files `copy_final_outputs()` had written on the actual pipeline run -- leaving two versions of the same file (e.g. both `CNGA3.ancestral_orf_walk.tsv` and `04_CNGA3.ancestral_orf_walk.tsv`) sitting side by side with no indication which one was current.
- Fixed the immediate output: removed every stray `03_<GENE>.*`/`04_<GENE>.*` file from all four genes' `final_results/<GENE>/important_output`, and re-copied/re-rendered everything using the same unprefixed naming `copy_final_outputs()` uses, with `05_plot_events.R` writing PDF+PNG straight into `final_results/<GENE>/important_output` in one step so the two formats cannot drift apart again.
- No script change was needed for this half of the bug -- `scripts/run_one_gene_00_to_04.sh`'s own `copy_final_outputs()` and its `--outdir "final_results/$GENE/important_output"` plot invocation were correct throughout; this was purely a consequence of an ad hoc, partial re-run bypassing them instead of reusing them.

## v3.35 - fixed a real shared-ancestry miscall in STOP character reconstruction

### Real bug: the terminal-stop Fitch fix above was still wrong when an unrelated lineage's insertion scattered a tip's real codon across non-adjacent alignment columns
- Reported directly with three independent real examples: CNGA3's `Carollia_perspicillata`/`Glyphonycteris_daviesi`/`Trinycteris_nicefori` (and their ancestors `Node35`/`Node36`) showed an isolated, split-looking stop; GUCA1B's Miniopterus ancestor showed `TAG` where the tips are mostly `TAA`; GUCA1C's `Node9`/`Node10` (ancestors of `Hipposideros_jonesi`/`abae`/`caffer`, all three genuinely `TGA`) reconstructed as `TTA`. The user asked directly how the alignment/reconstruction actually works and gave the fix in one sentence: *"For parsimony, use the tip to internode approach, i.e. look at the state of the shared tips to assign states to their common ancestor and keep doing it for each internode."*
- Root cause: the previous fix (the `fitch_reconstruct_nucleotide()` entry directly below) still read each tip's terminal codon off the shared alignment's fixed final 3 *columns*. A tip's real biological sequence has no gaps in it at all -- gaps only exist in the shared multiple-alignment coordinate system, inserted wherever *other*, unrelated lineages have extra length nearby. So a tip's own real, contiguous terminal codon can land split across non-adjacent columns (or short of the alignment's widest point entirely) purely because of insertions in completely unrelated taxa, with nothing wrong in the tip's own sequence. Reading a fixed column window -- even treating gap as a legitimate state, as the previous fix did -- reads whatever fragments happen to land in that window, not the tip's real ending.
- Fixed by replacing the column-based reconstruction with two tip-sourced functions in `scripts/04_ancestral_orf_walk.py`: `own_terminal_codon()` derives each tip's true terminal codon directly from *its own* sequence with all gaps stripped first (ignoring alignment column position entirely), and `fitch_reconstruct_codon()` runs standard two-pass Fitch parsimony (postorder down-pass building candidate states by intersection/union, preorder up-pass resolving ties toward the parent's state) treating a whole codon as one atomic state, bottom-up from tips to root -- exactly the tip-to-internode approach requested.
- The earlier "gap is a legitimate parsimony state" design (previous entry) is superseded, not just patched: a clade not reaching the alignment's widest column no longer means "unknown/gap" for that clade's ancestor, it means that clade's own real CDS is shorter, and its ancestor now gets *that clade's own* real reconstructed ending rather than a blank stand-in borrowed from an unrelated longer-tailed lineage elsewhere in the tree.
- `tests/backend_consistency_test.py::test_terminal_stop_gap_majority_not_overwritten_with_fake_stop` renamed to `test_terminal_stop_gap_majority_uses_each_lineages_own_true_ending` and its expectation corrected (a genuinely shorter-CDS clade's ancestor now reconstructs its own true ending, not a placeholder gap). New regression test `test_terminal_stop_scattered_by_unrelated_insertion_still_reads_correctly` reproduces the exact real-world shape of this bug (an unrelated tip's extra base splits two tips' real terminal codon across non-adjacent columns) and asserts the shared ancestor still reconstructs the correct, real codon.
- Verified against real data: re-ran `04_ancestral_orf_walk.py` for PDE6H, CNGA3, GUCA1B and GUCA1C. CNGA3's `Node35`/`Node36` now reconstruct `TAG`, matching all three tips' own true (gap-stripped) terminal codon. GUCA1B's real Miniopterus ancestors (`Node67`, MRCA of `Miniopterus_schreibersii`/`natalensis`, both `TAA`; `Node68`, their ancestor together with `Miniopterus_australis`, `TAG`) now both resolve to `TAA`. GUCA1C's `Node9`/`Node10` now resolve to `TGA`, matching all three `Hipposideros` tips' own true terminal codon, even though those same tips' raw alignment rows still show trailing gaps at the shared alignment's final columns (a longer-tailed, unrelated lineage elsewhere pushes the shared coordinate system further right -- expected and correct, since the fix no longer depends on that column position at all).
- Separately investigated and confirmed real, but out of scope for this fix: GUCA1C's `Hipposideros_abae`/`caffer` carry a single-column MACSE placeholder (`!`, rendered as `-` downstream) a few bases after the start codon that `Hipposideros_larvatus`/`armiger`/`swinhoii`/`jonesi` don't. Traced back to MACSE's own raw output and further back to the pre-alignment input sequences: all six taxa's real, pre-alignment sequences are identical for the first 14 bases (differing only by one synonymous substitution at base 15, no real indel anywhere nearby) -- so this is a genuine, unforced MACSE local-alignment artifact in a low-complexity stretch, not a real biological difference between these lineages. Unlike the terminal-stop case, this lives inside MACSE's own alignment output, upstream of everything Pensieve's own reconstruction controls, so it was not corrected as part of this fix -- see the open question below.

### Event marker colour palette unified between the two figures
- Reported directly: `<GENE>.event_map` (Figure 2) coloured its event bars by raw `biological_interpretation` via a completely separate, independently-computed palette (`interp_colours`) from `<GENE>.pseudogenization_tree` (Figure 1)'s `marker_class`/`marker_colours` -- the exact same event could render in a different colour depending on which of the two figures it was viewed in (e.g. any "deletion", pseudogenizing or not, shared or not, was always plain red in Figure 2, but could be black/gold/khaki/grey in Figure 1 depending on what kind of deletion it actually was). Figure 2's own old palette also had two categories, `premature_stop_gained` and `ambiguous_indel_change`, sharing the identical colour.
- Fixed by computing the marker classification and colour palette exactly once, immediately after events are loaded, shared by both figures -- Figure 2 now uses the identical `marker_class`/`marker_colours`/`marker_labels` Figure 1 does, so the same event is always the same colour in both places. `interp_colours` and the separate, ad hoc classification block that used to live inside Figure 1's section are both removed.
- Two more issues fixed as part of unifying this: (1) the marker classifier's "is this a confirmed premature-stop gain" check used `event_type == "stop_mask_gain"`, not `biological_interpretation == "premature_stop_gained"` as `04_ancestral_orf_walk.py` itself uses to decide a branch's own colour -- an *ambiguous* stop change could keep `event_type == "stop_mask_gain"` and be miscoloured as confidently disabling; now matched exactly. (2) Added two categories the unified palette was missing relative to the old Figure-2-only one: `ambiguous_disabling_candidate` (an unresolved STOP or frameshift-length change -- scoped the same way `04_ancestral_orf_walk.py` scopes "ORF-relevant" ambiguity, so an ambiguous *in-frame* indel's polarity tie still stays in the ordinary shared/other categories) and `stop_lost_reversion` (a premature stop being lost, i.e. a reversion -- meaningfully different from a stop being gained).
- Internal-node labels in Figure 2 (added earlier this release) were reported too small to read on a real rendered figure; doubled in size. Reported again after doubling: on a real rendered 103-tip figure the doubled labels now overlap each other; reduced by 30% (`size = max(1.4, tip_font * 0.7)`, from `max(2.0, tip_font * 1.0)`) -- smaller than the doubled size but still larger than the original.
- Verified against real data: PDE6H, CNGA3, GUCA1B and GUCA1C regenerated end to end and visually reviewed; the same six-category legend now appears in both figures for a given gene, with distinguishable colours throughout, and internal-node labels are legible without overlapping.

### Event map (Figure 2): every internal node is now labelled
- Requested directly: on a 100+ tip tree, `<GENE>.event_map`'s tree panel had no internal-node labels at all, making it impossible to tell which branch corresponds to which node in `alignment_events.tsv`/`orf_transitions_by_branch.tsv`, especially for deep/backbone nodes far from any tip.
- `<GENE>.pseudogenization_tree` (Figure 1) deliberately labels only event-carrying nodes -- it already has event markers and their position/length labels competing for the same vertical space per branch, and labelling all ~100 nodes there would reintroduce the fanning/clutter this release otherwise fixed. Figure 2's tree panel has no such competition, so every internal node now gets a small boxed label at its own branch position.
- Verified against real data: regenerated `event_map` for PDE6H, CNGA3, GUCA1B and GUCA1C (up to 411 events across 91 branches on a 103-tip tree) and visually reviewed; every backbone node is legible and distinguishable, including in the densest cases.

### Real bug: PAML cannot reconstruct the ordinary terminal stop codon at any ancestral node
- Reported directly, found by inspecting real ancestral FASTAs: `ancestral_integrated_alignment.fa` showed a non-stop codon at every internal node's true C-terminus (`AAG` in a real PDE6H run, `TTA` in a real GUCA1C run) where every intact tip has a real stop codon.
- Root cause: PAML's codon-substitution model has a 61-state, sense-codon-only state space by construction -- it cannot represent a stop codon as a state anywhere, ever. The terminal stop column is masked to `N` for *every* tip (not just ancestors) before codeml ever sees it, so PAML's own reconstruction there isn't just wrong, it's uninformed: nothing in its model or its input constrains it. This is a structural property of the codon model, not a bug in any particular codeml run, and it is distinct from premature stops -- those are already correctly handled by Pensieve's own stop-mask character/parsimony system, precisely because that system was built to reconstruct disabling lesions PAML can't represent. The ordinary, universal terminal stop was never covered by that system (correctly -- it isn't a disabling lesion), so nothing corrected it for internal nodes.
- Considered borrowing this from IndelMaP (an independent, non-codon-model reconstruction) but chose internal parsimony instead: Pensieve's design deliberately never lets IndelMaP overwrite the PAML nucleotide scaffold (concordance evidence only, to avoid mixing two ASR methods' assumptions on the coding backbone), and the same tip-observed data IndelMaP would use is already directly available here.
- Fixed in `scripts/04_ancestral_orf_walk.py` with a small, self-contained two-pass Fitch parsimony reconstruction (`fitch_reconstruct_nucleotide()`) applied to just the alignment's final 3 columns (the same span the rest of the pipeline already treats as "terminal"), sourced directly from the real observed tip alignment -- bypassing PAML's incapable codon model for only this span, leaving every other column exactly as PAML reconstructed it.
- A second issue found while verifying the fix against real GUCA1C data: at those literal final columns, most real alignments have only a handful of taxa with the longest aligned CDS carrying a real base; the rest are a trailing gap because their own true CDS ends earlier, at a column specific to their lineage. An initial version of the fix treated gap as "unknown data", which let those few informative tips over-confidently dictate a stop codon for every ancestor in the tree, including ones nowhere near the long-tailed lineage -- exactly the kind of overreach this release has otherwise been about removing. Fixed by treating gap as its own legitimate parsimony state alongside A/C/G/T, so most ancestors now correctly reconstruct "gap" (matching most descendants) and only nodes genuinely ancestral to the long-tailed taxa get a real reconstructed stop codon.
- Two new regression tests in `tests/backend_consistency_test.py`: `test_terminal_stop_reconstructed_at_ancestral_nodes` (the simple, uniform case) and `test_terminal_stop_gap_majority_not_overwritten_with_fake_stop` (the gap-majority case, reproducing the real GUCA1C shape).
- Verified against real data: re-ran `04_ancestral_orf_walk.py` for PDE6H, GUCA1C, CNGA3 and GUCA1B. PDE6H (every tip shares a real terminal stop) now shows `TGA` at all 101 internal nodes, matching every tip. GUCA1C (only 3 of 103 tips -- a small Emballonuridae clade -- reach the final alignment columns) now shows gap at 99 internal nodes and the real stop (`TAG`) only at the 2 nodes genuinely ancestral to that clade. Branch-level pseudogenization calls (`orf_transitions_by_branch.tsv`, `alignment_events.tsv`) are unchanged for all four genes -- confirmed identical `pseudogenization_branches`/`apparent_restorations` counts before and after -- since a wrong-but-still-non-stop codon never happened to trip the existing "internal stop" disruption check either way; this fix corrects the exported sequence data, not any classification.

### Real bug, found by reviewing a genuine full HPC run (PDE6H, 103 species)
- Reported directly: `Desmodus_rotundus` and `Diaemus_youngii` share the exact same premature in-frame stop codon (`TGA` at canonical alignment columns 160-162) and are sister taxa (their MRCA is `Node50`), but Pensieve reported it as two *independent* pseudogenization events on their two terminal branches, left `Node50` unlabeled, and the plot showed no shared origin at all.
- Root cause, in `scripts/03_alignment_events.py`'s `stop_mask_characters()`: a tip was only ever marked `PRESENT` for a stop character if *that tip's own* upstream raw-sequence classification pass (`02_...premature_stop_masking`) had independently flagged its occurrence as an `independent_stop_candidate`. `Desmodus_rotundus`'s occurrence was flagged that way; `Diaemus_youngii`'s identical occurrence, at the identical canonical-alignment columns, was instead caught by a later, generic "mask any remaining stop before building the PAML-safe alignment" catch-all pass (recorded as `terminal_or_unregistered_stop_masked_for_paml_only` in `PDE6H.premature_stop_registry.tsv`) and was therefore never added to the character's carrier set -- and the fallback branch that decided every other tip's state checked only "is this a clean ACGT triplet", never "does it actually equal the character's own stop codon". So `Diaemus_youngii` was silently recorded `stop_absent` for a character it demonstrably carries, corrupting the Sankoff parsimony reconstruction for that character (single-tip presence in each of two lineages, rather than the correct two-tip presence in one shared clade) and everything downstream: `Node50`'s reconstructed state, `alignment_events.tsv`'s `origin_node`/`shared_event`/`affected_tips`, `orf_transitions_by_branch.tsv`'s per-branch calls, and the plotted tree (which only ever labels/highlights an internal node when an event's `origin_is_tip` is `False`).
- Fixed by checking the canonical alignment directly, not just each tip's own upstream classification: once a stop character's (position, codon) is established by at least one independently-classified occurrence, every other tip whose own canonical-alignment segment at those exact columns matches that exact codon is now also marked present. A tip's own classification pass still decides whether its occurrence is trustworthy enough to *found* a character (avoiding spurious characters from noise); it no longer gates whether that tip carries an *already-established* one.
- New regression test (`tests/backend_consistency_test.py::test_shared_stop_not_missed_when_only_one_tip_is_registered`) reproduces the exact real-world shape of the bug (one registered carrier, one carrier caught only by the generic catch-all, both sharing the identical codon/position) and asserts the character, node states, and event table are all now correct.
- Verified against the real data: re-ran the full, un-subsetted 103-species PDE6H gene locally end to end. Before the fix, `Node50` reconstructed `stop_absent` and the event table showed two independent tip-origin `stop_mask_gain` events. After the fix: `origin_node=Node50`, `origin_is_tip=False`, `branch=Node51->Node50`, `shared_event=True`, `affected_tips=Desmodus_rotundus,Diaemus_youngii`; `orf_transitions_by_branch.tsv` now shows `Node51->Node50: confident_disabling_event` (the real origin) and both `Node50->Desmodus_rotundus` / `Node50->Diaemus_youngii` as `inherited_pseudogenic_history` (not new independent events); the rendered plot now labels `Node50`, draws the shared branch as the saturated-red first-disabling branch with its position/length annotated, and draws both descendant branches pale red for inherited history. `Diaemus_youngii`'s own additional lineage-specific indel events (unrelated frameshift, found independently) are still correctly plotted as separate secondary markers on its own branch -- this gene has both a shared ancestral lesion and an independent, lineage-specific one, and the corrected output distinguishes the two.
- No change to indel-character reconstruction, PAML/ASR handling, or plotting logic itself -- only the stop-mask character's per-tip presence/absence detection.

### Real bug: root-adjacent uncertainty silently poisoned entire, otherwise-clean trees
- Reported directly after reviewing GUCA1B, GUCA1C and CNGA3 plots: huge stretches of each tree rendered amber ("history uncertain") even on branches where both the parent and child node were separately, confidently reported `intact`.
- Root cause, in `scripts/04_ancestral_orf_walk.py`'s branch walk: the biological root never gets a fabricated sequence, so its own pseudogenic-history state starts `None` (unresolved) for every gene by design. Once a node's history was `None`, the walk had no path back to a confident "clean" state -- only to "confirmed disabling" -- so `history_uncertain_no_confirmed_event` (amber) silently propagated down *every* descendant branch indefinitely, regardless of how many generations of confidently intact ORFs came after it. On a gene with any ambiguous/frame-shifting character near the root (routine on a 100+ tip tree with hundreds of indel characters), this could amber-out the vast majority of an otherwise clean tree.
- Fixed: once a node's own reconstructed ORF is confidently `intact` and no disabling event (confident or ambiguous) is catalogued on the branch leading to it, its history now resolves to "clean" from that point forward, regardless of what remained unresolved above it. Remaining amber is now reserved for nodes whose own status is genuinely unresolved (e.g. `uncertain`), not nodes that are themselves confidently intact.
- A second, compounding bug in the same file: a node's *own* ORF-completeness call (`coding_status`) was marked `uncertain` if **any** character had an ambiguous/undefined reconstructed state at that node -- including ordinary in-frame indels, whose gap/residue polarity cannot change the reading frame or stop-codon content of the resulting sequence either way. Genes with many ordinary in-frame indels (not itself unusual) had this trigger constantly. Fixed to only treat ambiguity as ORF-relevant when the character is a STOP allele or a frame-shifting-length indel -- the only kinds of ambiguity that could actually flip an intact/disrupted call.
- Two new regression tests in `tests/backend_consistency_test.py`: `test_uncertain_root_resolves_to_intact_once_descendants_are_confidently_intact` and `test_ambiguous_in_frame_indel_does_not_force_node_uncertain`.
- Verified against real data: CNGA3 (118 events/65 branches) and GUCA1B (106 events/22 branches) re-run locally end to end. Before the fix, both plots were dominated by amber/orange branches wall to wall. After the fix, both trees are overwhelmingly grey (intact) with amber surviving only in the handful of places the underlying reconstruction is genuinely ambiguous (confirmed by cross-checking `orf_transitions_by_branch.tsv`: remaining amber branches have `child_coding_status=uncertain`, not `intact`).

### Plot readability and colour-coding overhaul (`scripts/05_plot_events.R`)
- Branch colours: six different `orf_transition` categories -- including several that are in fact *confirmed* disabling events, just with an unresolved first-loss footnote -- all rendered in the same amber, making a confirmed lesion visually indistinguishable from genuine uncertainty. Collapsed to four display categories (grey = intact, bold red = confirmed disabling, pale red = inherited, amber = genuinely uncertain/no signal either way); the legend now has 4-6 clear rows instead of 8-11 near-duplicates.
- Event markers (the vertical bars) were a single fixed colour (`grey15`/black) for every event type, with no legend, so "black bar" conveyed nothing about what kind of change it was. Markers are now coloured and legended by type, classification reusing the exact same fields `04_ancestral_orf_walk.py` uses to decide a branch's own colour (so a marker's colour is never in tension with its branch's): shared/ancestral in-frame insertion (non-disrupting) is blue; shared/ancestral in-frame deletion (non-disrupting) is khaki; an ordinary lineage-specific in-frame indel (the majority of events on most genes, and not itself notable) is a pale, deliberately receding grey rather than solid black, so it no longer reads as important. Pseudogenizing substitution (premature stop) and pseudogenizing indel (frameshift) initially reused the branch reds, but markers sit directly on top of branches that are *already* red (bold for confirmed-disabling, pale for inherited), so a same-family red marker on a red branch was nearly invisible; they are now near-black and gold respectively, chosen specifically to read clearly against both red branch shades.
- Event position/length text labels were anchored at a fixed offset that scaled with how many events shared a branch, so branches with many events (a real, common case -- not a bug) fanned labels vertically into neighbouring branches' rows, making them illegible and impossible to attribute. Labels now sit at a small, constant offset from their own branch's row (events on one branch are already kept apart along the x/time axis, so this is sufficient), with explicit `vjust` so the text grows away from its marker bar instead of centering on top of it. Marker bars themselves were also shrunk (roughly half the previous width/height) to leave more room.
- Removed the redundant `" bp"` suffix from every individual event-length label (unit is implied and was repeated hundreds of times per figure); Figure 2's axis title (a single, one-time label) still reads "Alignment position (bp)".
- The combined branch-colour + event-marker legend routinely ran past the figure width and got visually cut off, especially once the marker-type legend was added. Both legends now wrap to 2 rows (`guide_legend(nrow = 2, byrow = TRUE)`) with roughly half the previous key/swatch size and smaller legend text, so the full legend reliably fits within the figure.
- Verified against real data: CNGA3, GUCA1B, GUCA1C and PDE6H plots regenerated end to end and visually reviewed after each round of changes; labels are legible, the full legend is visible (not cut off), markers are distinctly coloured and clearly readable against both grey and red branches, and colours match their underlying event data (cross-checked directly against `alignment_events.tsv`, e.g. `Nycteris_thebaica`'s 18-event branch in GUCA1B has exactly 2 frame-shifting/pseudogenizing markers among 16 ordinary in-frame ones).

## v3.34 - automatic staged-install fallback for OOM-killed solves

### Real HPC failure: `bash install.sh --backend=mamba` killed during "Solving environment"
- Reported on a real Slurm cluster after `ml miniforge`: `mamba env create -n Pensieve -f environment.yml` was killed mid-solve (`environment: line 2: <pid> Killed`) with no package-conflict message -- the classic signature of the dependency solver being OOM-killed, not a real dependency error. This environment solves ~20 conda-forge/bioconda packages simultaneously, including a full R/tidyverse stack, which is a heavy single solve; HPC login nodes commonly cap per-process memory well below what that needs, even for mamba's faster libmamba solver.
- `install.sh` already shipped a lower-memory `--backend=staged` path (creates the env with minimal deps, then installs numpy/scipy/pandas/biopython/ete3, then the R stack, then macse/paml/muscle/iqtree/emboss as separate, much smaller solves) but nothing tried it automatically, so a user hitting this had to already know it existed.
- `auto`, `mamba`, and `conda` backends now automatically retry with the staged installer whenever the single-shot `env create`/`env update` fails or is killed (any non-zero, non-"command missing" exit status). `--no-staged-fallback` (or `PENSIEVE_NO_STAGED_FALLBACK=1`) disables this for scripting contexts that want strict single-attempt semantics. The final `auto`-backend error message, reached only if even the staged fallback fails on both mamba and conda, now explicitly explains the OOM-kill signature and suggests requesting an interactive compute allocation (`srun --mem=8G ... --pty bash`) and `conda config --set channel_priority strict`.
- Verified with a mock `mamba` that fails only the single-shot `env create` call (returns 137, as a real OOM kill would) and succeeds on every staged call: `--backend=mamba` now recovers and exits 0; `--backend=mamba --no-staged-fallback` correctly still fails without retrying. Added as a permanent regression in `tests/install_portability_test.sh`, not just verified ad hoc.

## v3.33 - local install hardening (no Slurm) and dynamic PAML dat-dir

### Installation fixes (found by actually running install.sh on a real Mac/conda setup)
- Fixed a real bug in `install.sh`'s final summary message: unquoted backticks around `` `--env-mode conda` `` / `` `run -n` `` inside an unquoted heredoc were interpreted by bash as command substitution, printing spurious `command not found: --env-mode` / `command not found: -n` errors on every successful install. Replaced with plain quotes.
- The closing banner no longer hard-codes a version string (`Pensieve v3.31 retains...`); it now reads `VERSION` at install time.
- Verified for real (not mocked) against this machine's actual `conda` (no standalone `mamba` binary present; `conda`'s libmamba solver backend is what's actually installed) using `install.sh --backend=auto`: environment update, IndelMaP clone via `git clone https://github.com/acg-team/indelMaP.git` into `external/indelMaP`, and re-running is idempotent (`git pull` instead of re-cloning). `paml` (4.10.10, bioconda) was already declared in `environment.yml` and installs/updates correctly.

### PAML dat/ directory is now auto-detected, never hard-coded
- `bin/pensieve` previously always passed `--dat-dir <package>/dat` (an empty placeholder; see `dat/README.md`) to the runner, regardless of where the actual environment's PAML installation was.
- Added `resolve_paml_dat_dir()`: for `--env-mode conda|mamba` it asks that exact environment (`<manager> run -n <env> python3 -c "...shutil.which('codeml')..."`) where `codeml` actually resolves and uses `<that prefix>/dat` if present; for `venv`/`inherit` it resolves from the venv or the calling shell's PATH. Falls back to the package placeholder only if `codeml` cannot be resolved. `--dat-dir` remains available to force a path explicitly. Dry runs never shell out to probe this.
- Confirmed empirically on this machine: bioconda's `paml` package does install a real `dat/` directory (`dayhoff.dat`, `wag.dat`, `grantham.dat`, etc., 18 files) directly under `<conda env prefix>/dat`, e.g. `/Users/.../envs/Pensieve/dat`. `bin/pensieve` now finds it there automatically with zero hard-coded paths.
- For transparency: those files are PAML's *empirical amino-acid* substitution matrices. Pensieve's own `templates/dummy_codon_asr.ctl` runs a codon model (`seqtype=1`, no `aaRatefile`) and does not currently read anything from this directory — confirmed by running the real `codeml` binary against Pensieve's exact control file with no `dat/` directory present anywhere nearby, which completed the same marginal ASR either way. `--dat-dir` is wired correctly for forward/CLI compatibility, not because today's run needs it.

### Two more real bugs found only by running a genuine end-to-end gene (PDE6H, 18 species, real codeml/MACSE/IndelMaP on this Mac)
- **rst parser (`scripts/03_integrate_asr_evidence.py`, shared by `02_validate_paml_marginal.py`): false "duplicate marginal node records" on every real multi-node run.** Real codeml output between the marginal DNA records and `(2) Joint reconstruction` includes an `Overall accuracy of the N ancestral sequences:` block and an `Amino acid sequences inferred by codonml.` block containing its own `Node #N  <one-letter AA string>` lines. The parser's stop-guard checked for the literal string `"Overall accuracy of the reconstruction"`, which never actually occurs in real PAML output (it's always `"Overall accuracy of the <N> ancestral sequences:"`), so parsing continued into the amino-acid block, where each `Node #N` line matched the same node-record pattern as the real DNA records and was recorded as a second, "duplicate" entry for every node -- hard-failing every real run at the ASR step (`Duplicate marginal node records in PAML rst: PAML_Node19, PAML_Node20, ...`). Fixed by matching on the real prefix (`"Overall accuracy of the"`) and adding `"Amino acid sequences inferred by codonml"` as an explicit second stop marker. Every previous mock/synthetic rst fixture in the test suite omitted this block, so this was invisible to mocked testing; a new regression test (`test_realistic_rst_amino_acid_section_not_treated_as_duplicate`) reproduces the real block structure so this can't silently return.
- **IndelMaP invocation (`scripts/02_run_asr_backends.sh` and `scripts/03_integrate_asr_evidence.py`) used the wrong CLI entirely.** `02_run_asr_backends.sh` called `indelMaP_ASR.py` with `--alignment`/`--tree`/`--output`, but the real script (confirmed via `--help` against a fresh `git clone` of `acg-team/indelMaP`) requires `--msa_file/-m`, `--tree_file/-t`, `--alphabet/-a` (`DNA` or `Protein`), and `--ancestral_reconstruction` to actually emit reconstructed sequences; every real run failed immediately with an argparse usage error and IndelMaP concordance was silently unavailable for every gene (harmless, since it's optional evidence, but never actually functioning). Separately, `03_integrate_asr_evidence.py` looked for output files named `<gene>.indelmap_ASR_tree.nwk` / `<gene>.indelmap_ASR_internal_ancestral_reconstruction.fas`, but the real tool appends its suffixes directly to `--output_file` with no `_ASR` infix (`<gene>.indelmap_tree.nwk` / `<gene>.indelmap_internal_ancestral_reconstruction.fas`), so concordance would have been reported "not_available" even after the CLI flags were fixed. Both fixed and verified end-to-end: PDE6H's real IndelMaP run now completes and reports 68/68 comparable states in agreement with Pensieve's own parsimony reconstruction (0 disagreements).

### Real end-to-end validation
- Ran the full pipeline (`diagnostics` through `plot`) locally against a genuine 18-species subset of PDE6H (from `bat1k_t2.nwk` + `PDE6H_bat1k.fa`), through the real `conda`-installed `Pensieve` environment (real MACSE, codeml, IndelMaP, R/ggplot2) -- no mocks. `codeml` reproducibly crashed during joint reconstruction on this exact machine/build too (not just large HPC trees), and was correctly recovered from via the validated marginal ASR from v3.32. The subset deliberately includes three species with genuine premature stop codons/frameshifts (`Desmodus_rotundus`, `Diaemus_youngii`, `Triaenops_persicus`) alongside their intact relatives (including `Diphylla_ecaudata`, an intact sister vampire-bat species); the pipeline correctly and independently flagged all three as `pseudogenization`/`confident_disabling_event` branches while every clean species (including `Diphylla_ecaudata`) remained `intact`, and IndelMaP's independent reconstruction agreed on all 68 comparable structural states.

## v3.32 - stop codeml before joint reconstruction

### PAML ASR reliability
- `02_run_asr_backends.sh` now launches `codeml` in the background and polls its `rst` file instead of blocking on the whole process. As soon as `rst` shows PAML has finished writing `(1) Marginal reconstruction of ancestral sequences` and moved on to the `(2) Joint reconstruction of ancestral sequences` marker, Pensieve stops `codeml` itself (SIGTERM, then SIGKILL after a short grace period) rather than waiting on — and risking an OOM kill/crash during — the far more memory-hungry joint reconstruction that Pensieve never reads or needs.
- This directly fixes real-world failures observed on large gene trees (e.g. `GUCA1B`, `GUCA1C`, `PDE6H` in a bat-vision-gene run) where `codeml` was killed mid-`Joint reconstruction` and the run failed with `[ERROR] Failed at line 218` even though the marginal ASR Pensieve actually needs had already completed.
- The stop is provably safe: a stream is written to disk in order, so if the `Joint reconstruction of ancestral sequences` marker is observable in `rst`, everything written before it (the entire marginal section) is guaranteed to already be flushed.
- The new `02_T.codeml_run_status.tsv` value `STOPPED_BEFORE_JOINT_RECONSTRUCTION_WITH_VALIDATED_MARGINAL_ASR` distinguishes a deliberate early stop from the pre-existing `CONTINUE_WITH_VALIDATED_MARGINAL_ASR` recovery (codeml exited non-zero entirely on its own). Both are non-fatal and both require the same validated marginal ASR.
- Poll interval defaults to 10s and is configurable via `PENSIEVE_ASR_POLL_SECONDS` for testing or site tuning.
- Verified against a new regression test (`tests/paml_exit_and_reference_free_test.py::test_slow_codeml_stopped_before_joint_reconstruction`) using a mock `codeml` that would otherwise hang for 60s in "joint reconstruction"; Pensieve now stops it in ~2s once the marginal section is complete, including correct parsing of codon-spaced sequences (e.g. `node #7 ATG AAA CCC GGG`).
- No change to the `n_tips - 2` fix, the rst-declared-node-range parser, or the whitespace-tolerant codon-sequence parser introduced in v3.30 — those were re-verified, remain correct, and are unaffected by this release.

### Known issue outside this package
- Failing jobs in the reported cluster run were submitted via a bare `conda run <script> ...` wrapper (no `-n <env>`), which does not match anything `bin/pensieve` itself generates for any `--env-mode`. That wrapper is external to this package (not present in this repository) and should be replaced by calling `bin/pensieve --mode slurm --env-mode inherit` (if `mamba activate Pensieve` is already run before submission) or `--env-mode mamba --env-name Pensieve` (to let Pensieve run `mamba run -n Pensieve ...` itself). See the accompanying report PDF for details.

## v3.31 - portable installation and runtime environment handling

### Installation portability
- Removed any assumption that mamba or Miniforge exists. `install.sh` supports `auto`, `mamba`, `conda`, `staged`, `venv`, and `current` backends.
- No installer path executes `ml miniforge`, `module load`, `mamba activate`, or `conda activate`. Site-specific module loading is left to the user/HPC site.
- Added `requirements-pip.txt` for Python-only installations. `--backend=venv` creates a standard Python virtual environment; `--backend=current` installs Python dependencies into the current Python without creating any environment.
- venv/current explicitly report that MACSE, codeml/PAML and R/Rscript are external scientific runtimes that must be supplied separately.

### Local and Slurm execution
- Added `--env-mode inherit|conda|mamba|venv` (default `inherit`). Local execution therefore works in any already prepared environment without an activation command.
- Slurm scripts no longer hard-code `ml miniforge` or `mamba activate Pensieve`. `inherit` uses the exported submission environment; explicit conda/mamba modes use `<manager> run -n ENV`; venv mode sources the requested virtual environment.
- Added optional `--slurm-module NAME` for clusters that genuinely require a module. Pensieve never invents or assumes the module name.
- Preserved the single-allocation Slurm design and all v3.30 scientific inference code unchanged.

### Validation
- Added portability regression tests for installer help/dry-runs and local/Slurm environment-mode generation.
- Verified all files under `scripts/` are byte-identical to v3.30, so this release changes installation/execution plumbing only.

## v3.30 - rst-driven PAML validation + compensated-STOP and audit cleanup

### PAML safe recovery
- Replaced the erroneous `n_tips - 2` ASR-completeness assumption with PAML's own `rst` declaration: `Nodes X to Y are ancestral`; the required marginal count is exactly `Y-X+1`.
- Validation reads only `(1) Marginal reconstruction of ancestral sequences`, requires every declared `node #X..#Y` exactly once at the expected alignment length, and ignores joint-ASR records for core completeness.
- A non-zero codeml exit after complete marginal ASR is therefore a warning/continue condition; a missing or malformed declared marginal sequence remains fatal.
- The validator and downstream integrator now use the same parser, eliminating competing interpretations of real `rst` files.
- Added an exact 103-tip regression fixture matching `Nodes 104 to 205 are ancestral`: 102 marginal records + 102 joint records = 204 global `node #` lines, while only the 102 marginal records define Pensieve success.
- All PAML-declared marginal sequences are retained in an audit FASTA. Any TreeView serialization/root-only PAML vertex is separated during topology mapping and never used to fabricate a biological-root sequence.

### STOP/frame logic
- Replaced the old “any upstream MACSE marker” rule with frame phase at the STOP position. Upstream MACSE partial-codon correction lengths are summed modulo 3.
- A STOP in non-zero phase remains a likely frameshift consequence; a STOP after compensated/restored phase (`0 mod 3`) remains eligible as an independent allele-specific nonsense event.

### Auditability and ORF-history interpretation
- Removed the unreachable legacy integration functions identified during code review, including obsolete IndelMaP/PAML projection and internode-registry scaffolding.
- Replaced generic branch-level `uncertain` output with explicit categories and added `transition_evidence` / `uncertainty_reason`.
- Root-adjacent sequence-only disruption is now explicitly reported as `root_adjacent_disruption_first_loss_unresolved`; Pensieve remains conservative and does not fabricate a root sequence or claim an unsupported first-loss branch.

### Plot validation
- Added `tests/plot_smoke_test.sh`, which executes the R plotting script for dated and undated synthetic cases whenever `Rscript` is available and validates non-empty PDF/PNG/event outputs. The packaging environment lacks Rscript, so this test reports an explicit SKIP here and will execute in the normal Pensieve R-enabled environment.

## v3.29 - phylogenetic ordering for MSA/manual inspection

### Sequence ordering
- Removed alphabetical FASTA ordering from Step 00. Pruned tip FASTAs now follow the left-to-right terminal order of the rooted user tree.
- MACSE NT and AA outputs are explicitly rewritten into that phylogenetic tip order, including when existing MACSE outputs are reused.
- Canonical native, PAML-safe, IndelMaP-input, AA and PHYLIP alignments inherit the same tree order.
- PAML marginal ancestor FASTA is written in rooted phylogenetic preorder after topology-based node mapping rather than by PAML node number/alphabetical label.
- Combined observed-tip + reconstructed-internode MSA is written in deterministic depth-first preorder, placing each available ancestor immediately before its descendant clade.
- `--alignment defined` still preserves every alignment column exactly; only FASTA record order is changed for inspection.

### New inspection outputs
- Added `GENE.phylogenetic_msa.fasta` to `important_output/` for direct manual inspection of tips and reconstructed internodes together.
- Added `GENE.phylogenetic_tip_order.tsv` and `GENE.phylogenetic_sequence_order.tsv` to make row order explicit and reproducible.
- The order table records node type, parent, depth, sequence availability and descendant tips. The clock=0 biological root is listed but never given a fabricated sequence.

### Validation
- Added synthetic tests using deliberately non-alphabetical tree order (`C,A,D,B`) and a combined internode/tip preorder case.
- All prior reference-free, breakpoint, PAML-recovery, root-policy and pseudogenic-history smoke tests remain in the suite.


## v3.28 - reference-free cleanup + safe marginal-PAML recovery

### Complete removal of reference dependency
- Pensieve remains entirely reference-free: no biological reference species/sequence is selected or consulted anywhere in event discovery, coordinates, ASR, plotting, or diagnostics.
- All lesion coordinates are 1-based/inclusive canonical-alignment coordinates; `GENE.canonical_alignment.fasta` is exported to the important final output.
- The runner now actively deletes the known legacy `reference_info`, `reference_sequence`, and `macse_indels_relative_to_reference` artifacts when an older work directory is reused, preventing stale files from implying a reference-based analysis.

### Safe codeml non-zero handling
- Added `scripts/02_validate_paml_marginal.py`, which strictly validates the actual Pensieve dependency: one expected-length marginal ancestral sequence for every clock=0 PAML internal vertex (`n_tips - 2`), together with a parseable PAML-labelled tree.
- The backend removes stale PAML products before every new codeml attempt, so validation can never accidentally accept an `rst` from an earlier run.
- A non-zero codeml exit is now **non-fatal only when the current-run marginal ASR passes strict validation**. This safely handles PAML runs that finish marginal reconstruction and then fail during the optional joint-reconstruction phase.
- If marginal ASR is incomplete, missing, wrong-length, or unmappable, Pensieve still fails immediately.
- Added `02_GENE.codeml_run_status.tsv` and `02_GENE.paml_marginal_validation.tsv`; both are copied to final supporting output.
- Resume/integration paths also validate an existing `rst` before accepting it.

### Validation
- Added smoke tests for codeml exit 1 after a complete marginal ASR (must continue), codeml exit 1 with an incomplete marginal ASR (must fail), reference-free Step 00/01 contracts, and deletion of stale legacy reference artifacts.

## v3.27 - fully reference-free coordinates and diagnostics

### Reference dependency removed
- Removed `--reference` from Step 00 and stopped selecting any automatic or user-specified reference species.
- Removed `00_GENE.reference_info.tsv` and `00_GENE.reference_sequence.fasta`; they are no longer created.
- Removed `--reference-species` from Step 01.
- Removed `01_GENE.macse_indels_relative_to_reference.tsv` and all `reference_species`, `reference_nt_start`, `reference_nt_end`, `deletion_relative_to_reference`, and `insertion_relative_to_reference` fields/labels.
- MACSE `!` diagnostics now contain only species-local provenance plus MACSE alignment coordinates and explicitly assign no insertion/deletion polarity.
- All evolutionary indel direction is inferred only by the rooted-tree event engine from canonical alignment states.
- The canonical native alignment is copied to `final_results/GENE/important_output/GENE.canonical_alignment.fasta` and is the sole coordinate authority for event/STOP tables and plots.
- Added smoke-test guards that fail if reference-selection/reference-relative output contracts reappear.

## v3.26 - coherent breakpoint-event / PAML integration

### Pipeline architecture
- Rewired the executable order to `diagnostics -> alignment -> events -> asr -> integrate -> plot`.
- Fixed the v3.25 production failure in which `02_prepare_asr_inputs.py` was never called before `02_GENE.codon_for_paml.phy` was required.
- Removed legacy authoritative dependencies on `candidate_indel_frameshift_events.tsv` and `paml_indelmap_asr_combined.fa`.
- Event reconstruction now occurs before ancestral-sequence integration and actually supplies the structural states used in the final ancestors.

### Canonical alignment
- `--alignment perform`: MACSE NT alignment is now the single canonical coordinate system; the second MUSCLE alignment/projection layer is removed.
- `--alignment defined`: the user alignment is authoritative; Pensieve never inserts/reorders columns. Invalid codon-alignment length now fails rather than being silently padded.
- Two synchronized views are written: native structural (`! -> -`) and PAML-safe (`! -> N`, STOP -> `NNN`).
- Corrected the v3.25 documentation/logic claim that every MACSE `!` is a one-base deletion. v3.26 treats `!` as a partial-codon placeholder; insertion/deletion direction comes from aligned occupancy + phylogeny.
- Removed the premature-STOP-only gate that could exclude frameshift-only taxa from MACSE-guided representation.

### Event reconstruction
- Retained the v3.25 breakpoint decomposition: shared core + lineage-specific extension + independent interior events.
- Fixed state coding so a smaller interior gap is ABSENT for a larger event when known residues occur inside the larger interval, while a different larger deletion spanning an entire smaller character remains UNKNOWN. This is the GUCA1B Miniopterus/Nycteris distinction.
- Exact parsimony ties remain biologically ambiguous. `--tie-break` only chooses a representative history; it can no longer turn a tied root polarity into a confident deletion/insertion.
- STOP characters are allele-specific (`TAA`, `TAG`, `TGA`).
- Signed frame arithmetic separates deletion (negative) from residue gain/insertion (positive), and premature STOPs no longer count as frame changes.

### PAML / IndelMaP / nodes
- Kept `clock=0`, `fix_blength=0`, `RateAncestor=1`, `cleandata=0`. IQ-TREE is not part of the core ASR path.
- PAML node labels are mapped to deterministic Pensieve `Node<i>` labels; the biological root is explicitly `NA` for PAML marginal sequence under clock=0 and is never fabricated.
- IndelMaP is optional concordance only and never supplies authoritative ancestral gaps. Its failure warns and core inference continues.
- Removed nested PAML Slurm child jobs; in top-level Slurm mode codeml runs inside the parent allocation, eliminating the orphan/unfinished-child failure mode.

### ORF/history logic
- Final native ancestors now start from the PAML substitution scaffold and receive Pensieve structural/STOP states directly.
- ORF status is `intact`, `disrupted`, `uncertain` or `unavailable`; missing terminal STOP remains acceptable.
- Pseudogenic history is sticky. A compensatory indel can yield `apparent_orf_restoration` but does not turn the lineage grey/functional or create a second first-loss event.
- A frameshifting GAP state at the sampled root is not automatically called a pre-root deletion because ancestral absence + later insertion is an alternative history.

### CLI, resume and testing
- `--indelmap`, `--tie-break` and `--breakpoint-tolerance` are now actually forwarded by the top-level CLI.
- New canonical stages are documented while historical step names remain aliases.
- Resume checks prerequisites rather than requiring files to be newer than the current invocation.
- Added synthetic runner, GUCA1B, root-tie, insertion, STOP-allele, frameshift-only, defined-alignment, mock-PAML/root-policy and compensatory-history tests.
- `bash tests/smoke_test.sh` passes in the packaging environment. Real MACSE/codeml/IndelMaP/R execution still requires acceptance testing in the Pensieve conda environment.

> The v3.25 changelog below is historical. Where it conflicts with v3.26 (especially the interpretation of MACSE `!`), v3.26 supersedes it.

## v3.25 - breakpoint-coded events, no IQ-TREE, optional IndelMaP

### Event reconstruction (rewritten)
- Indels are coded by BREAKPOINTS, not by alignment columns
  (`scripts/03_alignment_events.py`). v3.24 cut a new block whenever the exact
  carrier set changed between adjacent columns, so any unrelated lineage with an
  overlapping indel split someone else's event. On GUCA1B that fragmented the
  42 bp Miniopterus deletion into three blocks, and in general it can turn one
  in-frame deletion into several apparent frameshifts.
- A run that EXTENDS a shared deletion is decomposed into the shared core plus
  its lineage-specific extension. GUCA1B 646-696 in M. australis becomes
  646-687 (shared by all three Miniopterus) plus 688-696 (australis only).
- A contained run sharing NEITHER breakpoint is an independent interior event
  and splits nothing: Nycteris thebaica 661-669 stays on its own branch.
- A MACSE '!' is now treated as what it is - a missing nucleotide. The native
  alignment renders it as a gap, so a shared frameshift marker is an ordinary
  shared 1 bp deletion. The whole `aligned_macse_frameshift_marker` event class,
  its carrier bookkeeping and the destructive mask-projection step are gone.
- delta-parsimony support is reported per node. delta = 0 is an exact tie and is
  flagged; `--tie-break` chooses the display resolution but never hides the tie.
- Frame arithmetic per node is computed without any ancestral sequence and is
  cross-checked against the sequence-based ORF call.

### Ancestral sequences and the ORF walk (new)
- `scripts/04_ancestral_orf_walk.py` relabels codeml node numbers to Pensieve
  `Node<i>` labels by rooted descendant-tip set, so tree, FASTA and every table
  share one namespace. `03_<gene>.paml_to_pensieve_node_map.tsv` records it.
- Two ancestral alignments: `ancestral_corrected.fa` (frame-restored, what
  codeml saw) and `ancestral_native.fa` (lesions re-applied). codeml emits a
  base in every column including ones the ancestor never had; those are
  overwritten by the reconstructed gap state.
- Root-to-tip ORF walk. The branch where intact -> broken is the
  pseudogenization event; descendants are `already_pseudogenic`.

### Removed
- IQ-TREE, and with it `02_root_iqtree_for_indelmap.py`,
  `02_prepare_fixed_iqtree_tree.py`, `02_validate_fixed_iqtree_tree.py`,
  `02_prepare_indelmap_runtime.py`, `indelmap_calculateC_full_precision.py`
  (1,023 lines). codeml estimates its own branch lengths.
- `02_build_alignment_event_candidates.py` and
  `03_project_alignment_mask_events.py`, superseded by the new engine.

### Changed
- `templates/dummy_codon_asr.ctl`: `fix_blength = 0` (codeml estimates branch
  lengths), `cleandata = 0` retained, `method = 1` for speed on large trees.
- IndelMaP is optional (`--indelmap yes|no`, default no) and is a concordance
  check only. If it is missing the run warns and continues.
- Figures are named after the gene and only event-bearing internodes are
  labelled. Event labels are the alignment start above the branch and the length
  below. Branches are coloured by ORF status: saturated red where the gene
  broke, pale red for descendants, grey for intact.
- Figure height, tip font and the tip-label strip scale with the tip count and
  the longest name. v3.24 hard-coded `--height 10`, which at 103 tips gives
  6.99 pt per row while tip labels are 9.10 pt - overlap was guaranteed.
- Help: `-h`, `-help`, `--help` give the summary; `--help --long` (or `-long`)
  gives the full manual.

### Not included
- Branch-model omega and Meredith-style dating of pseudogenization are deferred.

## v3.24 — Alignment-block Sankoff parsimony as the authoritative event-placement layer

- Adopted the strongest part of the reviewed Claude implementation: shared-event origins are reconstructed from vertical aligned tip states on the dated tree, independently of PAML/IndelMaP ancestral sequences.
- Candidate events are reconstructed at the **block level**, so a multi-base gap or masked STOP is one evolutionary character rather than several independent sites.
- Added equal-cost binary Sankoff parsimony with explicit gain, loss/reversal, root-state and ambiguity reporting.
- Reversals no longer fragment a shared ancestral event into multiple tip events.
- If a gap state is ancestral at the root, Pensieve reports descendant gap losses as candidate insertions/restorations rather than inventing a giant root-level deletion.
- For gap blocks, other non-ACGT states are unknown rather than evidence of gap absence.
- For `N`, only correction-supported masked STOP/MACSE frameshift carriers are treated as event presence; natural/unrecorded `N` remains unknown.
- PAML/IndelMaP are now cross-check/evidence layers for event placement rather than the authority that discovers shared events. Gap transitions report IndelMaP concordance.
- Correction-supported N events are projected according to the parsimony reconstruction, including the correct shared ancestral nodes.
- Added `03_GENE.alignment_block_parsimony_summary.tsv`, `...events.tsv`, `...node_states.tsv`, and warnings/projection audits.
- The final event table replaces legacy ASR-driven placements for alignment candidates with the vertical alignment-parsimony placements.
- Verified against the uploaded PDE6H dataset: the Desmodus/Diaemus masked STOP at columns 160–162 is one gain on `PAML_Node142 -> PAML_Node143`; Triaenops column 178 is one terminal gain; the three majority-gap blocks are correctly polarized as single-species gap losses/restorations.

## v3.23 — Alignment-column event reconstruction and complete event reporting

- Event identity is now based on primary-alignment columns, never species-specific ungapped codon positions.
- Shared masked premature STOPs are merged across species even when upstream insertions shift raw codon numbering.
- Correction-supported MACSE `! -> N` frameshift markers are promoted from diagnostics to event candidates.
- A vertical reference-free audit records every variable gap block and every aligned `N` block.
- Recorded STOP/frameshift masks are projected as binary events onto maximal uniform carrier clades and overlaid onto matching PAML-labelled ancestral sequences with explicit provenance.
- Every gap event is reportable and plottable, including singleton, in-frame, nondeleterious and one-base events.
- Gap events lacking an unambiguous ASR transition receive an explicitly labelled observed-carrier fallback instead of disappearing.
- Plot labels display `NodeN` rather than `PAML_NodeN`; machine-readable labels remain unchanged.
- Raw IndelMaP internal ancestral reconstruction FASTA and ASR tree are collected into `supporting_files`.
- IQ-TREE fixed-topology runs include `--keep-ident`.
- Added actual PDE6H regression coverage: the Desmodus/Diaemus STOP at aligned columns 160–162 maps once to `PAML_Node142 -> PAML_Node143`, and the Triaenops one-base marker is retained.


## v3.22 — Direct fixed-topology IQ-TREE tree and correct rooted/unrooted node model

### Tree topology and branch lengths

- Replaced the post-IQ-TREE branch-length-transfer workflow.
- The rooted user tree remains the biological reporting authority.
- The dated degree-2 root is suppressed before IQ-TREE and CODEML `clock = 0`.
- Added a complete user-tree registry containing one dated-root edge row and `n - 2` non-root tripartition rows.
- IQ-TREE 3 is run with `-t FIXED_TREE --tree-fix`.
- The raw IQ-TREE `.treefile` is copied byte-for-byte for PAML.
- Pensieve no longer serializes IQ-TREE lengths through Bio.Phylo or rounds small positive values to zero.
- Renamed the validator to `02_validate_fixed_iqtree_tree.py` because it validates and copies; it does not transfer lengths.

### PAML

- Explicitly set `clock = 0` and use the unrooted degree-3 IQ-TREE tree.
- The expected PAML internal-node count is `n - 2`, not `n - 1`.
- The dated root is no longer treated as a PAML node and receives no invented PAML sequence.
- Set `fix_blength = 1`, using IQ-TREE lengths as CODEML starting values rather than fixing GTR distances under the codon model.
- Clarified that the outermost PAML trifurcation is a Newick serialization vertex, not the dated root.

### IndelMaP

- Added `02_root_iqtree_for_indelmap.py` to insert the dated root on the exact corresponding IQ-TREE edge.
- Every non-root IQ-TREE branch-length token is preserved verbatim.
- The root edge is divided with Decimal arithmetic and audited; the two output lengths must sum exactly to the original edge.
- Added an isolated full-precision IndelMaP runtime that removes upstream four-decimal percentile rounding, preventing valid tiny branches from becoming zero and causing `log(0)`.
- The raw IQ-TREE tree itself is never floored or modified to work around IndelMaP.

### Internode mapping

- Non-root internodes are matched across user tree, IQ-TREE, PAML and IndelMaP by exact root-independent three-way tip partitions.
- Rooted descendant species are reported after mapping back onto the dated user tree.
- The dated root is represented as `UserRoot`, with PAML label `NA` and its actual IndelMaP root label.
- Added explicit root-state provenance and separate `03_GENE.indelmap_dated_root_asr.fa`.
- Removed unused legacy rooted-descendant and rerooting functions from the integration script.

### Shell, output and documentation

- Retained strict-shell-safe `final_root` initialization.
- Retained plot-only final-output recollection.
- Updated CLI help, README, validation tests and resume guidance.
- Removed obsolete empirical amino-acid `.dat` installation instructions.

## Earlier versions

Versions 3.15–3.21 are superseded for tree-dependent Results 02–04. Their historical notes are intentionally omitted from this clean changelog because several described node/root models are no longer valid.
