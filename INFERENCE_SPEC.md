# Pensieve v3.35 inference specification

This file defines the scientific contract implemented by the code. If code and documentation disagree, this specification is the intended behaviour and the discrepancy should be treated as a bug.

## 1. Authorities

1. The **rooted user tree** is the biological topology/reporting authority.
2. The **canonical alignment** is the sole coordinate authority. Pensieve never selects or uses a reference species/sequence for event discovery, coordinates, polarity, diagnostics, ASR, or plotting.
3. Pensieve's **breakpoint/STOP character reconstruction** is the authoritative structural-event history.
4. **PAML** is the nucleotide-substitution ASR scaffold, not an indel-history estimator.
5. **IndelMaP** is independent concordance evidence and must never silently overwrite authoritative structural states.
6. A codeml process exit code does not by itself define ASR success. The `rst` line `Nodes X to Y are ancestral` is the sole source of truth for how many PAML ancestral sequences that run declares: `Y-X+1`. Pensieve validates only section `(1) Marginal reconstruction of ancestral sequences`, requiring every declared `node #X ... node #Y` exactly once at the expected alignment length. A later joint-reconstruction failure is non-fatal when this marginal section is complete. Pensieve never guesses ASR completeness from tip count or `n-2`.
7. `Node<i>` labels belong to Pensieve. PAML/IndelMaP labels are crosswalk metadata.

## 2. Canonical alignment

- `perform`: MACSE NT alignment is canonical.
- `defined`: user alignment is canonical and must not be altered.
- There is never a second global MSA in v3.31.
- Every reported event/STOP coordinate is 1-based and inclusive in this canonical alignment.
- The exact canonical native alignment is copied into the important final output directory so coordinates can be inspected directly.

Two synchronized views share identical columns:

- native: biological residues/gaps; MACSE `!` placeholder rendered `-`;
- PAML-safe: MACSE `!` rendered `N`; exact STOP codons masked `NNN`.

`! -> -` is a representation change only. It is not an insertion/deletion polarity call.

## 3. Indel characters

Per-tip maximal gap runs are decomposed by breakpoint relationships.

- matching breakpoints -> same character;
- one run extending another -> shared core + extension;
- contained run sharing neither breakpoint -> independent interior event.

For a character interval in a tip:

- matching decomposed gap run -> `GAP`/present;
- at least one known residue inside the interval -> `RESIDUE`/absent;
- whole interval unavailable because a different larger gap spans it -> `UNKNOWN`;
- unresolved N/no informative residue -> `UNKNOWN`.

This rule intentionally makes a small interior deletion ABSENT for a larger event, while a larger encompassing deletion is UNKNOWN for the small interior character.

## 3a. Phylogenetic row order

Sequence-record order is a visualization/manual-inspection property, not an inference prior. Pensieve therefore derives it only from the rooted pruned user tree and never from species names. Tip-only FASTAs follow left-to-right terminal order. Combined tip+internode FASTAs use depth-first preorder, placing each reconstructed ancestor immediately before its descendant clade. Reordering rows must never alter alignment columns, coordinates, event states, or likelihood calculations.

For `--alignment defined`, the user's alignment columns are immutable; only the order of FASTA records may be rearranged to follow the tree.

## 4. STOP characters

Only mapped raw premature in-frame stops are candidates. Terminal STOPs are never pseudogenizing events merely because they are terminal.

STOP event identity includes the exact allele. `TAA`, `TAG`, and `TGA` are separate mutational characters even at the same aligned position.

STOP classification uses the reading-frame phase at the STOP position, not merely the presence of any earlier MACSE marker. Pensieve sums the lengths of upstream MACSE partial-codon `!` runs modulo 3:

- non-zero phase -> the STOP is retained diagnostically as a likely frameshift consequence and is not promoted to an independent nonsense character;
- phase restored to zero (for example, compensating upstream frame corrections whose lengths sum to a multiple of three) -> the STOP remains eligible as an independent allele-specific nonsense character;
- if phase classification is unavailable/uncertain, the STOP must be retained with uncertainty rather than silently discarded.

## 5. Parsimony

Default transition costs are equal.

Events are transitions on parent->child branches; carrier MRCA alone is not an event origin.

Exact equal-cost state assignments are uncertainty. `--tie-break` can choose a representative assignment but cannot alter:

- `ambiguous_origin=True`;
- `direction_confident=False`;
- ambiguous biological interpretation.

## 6. Direction

For a confident indel character:

- `RESIDUE -> GAP`: deletion;
- `GAP -> RESIDUE`: insertion/residue restoration.

The second category remains intentionally cautious because simple gap-state parsimony cannot always distinguish a true insertion from regain/restoration after a prior loss outside the sampled history.

## 7. Frame arithmetic

Only confident indel transitions contribute.

- deletion = negative event length;
- residue gain/insertion/restoration = positive event length;
- current frame offset = signed cumulative change mod 3.

Premature STOPs never alter frame offset.

A current offset of zero does not erase the number/history of earlier frameshifting events.

## 8. PAML

Required core settings:

```text
clock = 0
fix_blength = 0
RateAncestor = 1
cleandata = 0
```

PAML estimates its own branch lengths. PAML posterior probabilities apply to its substitution scaffold, not to the final lesion-overlaid native ancestor.

PAML ASR completeness and biological node mapping are separate questions. Pensieve first validates every ancestral sequence declared by `Nodes X to Y are ancestral`. During topology mapping, a PAML TreeView serialization/root vertex may be present in addition to the degree-3 biological internodes; such a sequence is retained in an all-declared audit FASTA but is not used to fabricate a distinct biological-root sequence.

The biological root has no distinct marginal PAML sequence under the unrooted `clock=0` interpretation and must not be fabricated.

## 9. Native ancestral CDS

For each reconstructable internal node:

1. start from PAML marginal nucleotide sequence;
2. overlay Pensieve `GAP` states as `-`;
3. keep PAML bases for Pensieve `RESIDUE` states;
4. restore exact reconstructed STOP alleles;
5. flag ambiguous structural states rather than guessing;
6. remove alignment gaps to obtain the native CDS.

## 10. ORF state

Missing terminal STOP is acceptable.

Definite disruption is based on:

- missing start `ATG`;
- CDS length not divisible by 3;
- internal exact in-frame STOP.

Unknown bases/ambiguous structural states produce `uncertain` rather than automatically intact/disrupted when they prevent a definitive call.

## 11. Pseudogenic history

ORF state and history are not interchangeable.

- first confidently acquired disabling lesion on an otherwise no-loss history -> `pseudogenization`;
- descendants -> `already_pseudogenic`;
- apparent intact ORF after known earlier loss -> `apparent_orf_restoration`, while inherited history stays true;
- uncertain entering history -> do not claim a first-loss branch.

Branch-history output must distinguish different uncertainty mechanisms instead of collapsing them into one generic `uncertain` label. In particular:

- `confirmed_disabling_event_first_loss_unresolved`: a real disabling event is present, but entering pseudogenic history is unresolved;
- `root_adjacent_disruption_first_loss_unresolved`: a sequence-only disruption is seen immediately below the biological root, but the absent distinct clock=0 root sequence prevents assigning that branch as the first loss;
- `sequence_disruption_first_loss_unresolved`: a disrupted child is observed while entering history is unresolved;
- `ambiguous_disabling_event`: event reconstruction itself is ambiguous;
- `history_uncertain_no_confirmed_event`: entering history is unresolved and there is no confirmed disabling event;
- `sequence_state_uncertain` / `sequence_state_unavailable`: sequence evidence itself prevents a definitive call.

`transition_evidence` and `uncertainty_reason` must accompany these states so the reason is explicit.

A compensatory indel is therefore never interpreted as proof of biological resurrection.

A frameshifting GAP state reconstructed at the sampled root does **not** by itself prove a pre-root deletion: it may instead represent ancestral absence before a residue insertion. Pre-root indel-based pseudogenic history is therefore not asserted from root gap occupancy alone.

## 12. IndelMaP

IndelMaP may agree or disagree. Both are informative.

Its role is to write concordance evidence. A missing/crashed IndelMaP run must not prevent the core PAML+Pensieve reconstruction from completing.

## 13. Plotting

- strong red: first confident disabling branch;
- pale red: inherited pseudogenic history;
- amber: unresolved first-loss/history state;
- apparent restoration: pale-red inherited history plus explicit restoration annotation, never a grey functional reset.

Only event-bearing internodes need event labels.
