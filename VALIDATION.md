# Pensieve v3.35 validation

## Reference-free invariant

The smoke suite verifies that Step 00/01 expose no reference-selection arguments or reference-relative output contracts, that no `reference_info`, `reference_sequence`, or `relative_to_reference` files are created, and that the canonical alignment is copied to `final_results/.../important_output/GENE.canonical_alignment.fasta`.


## Safe codeml exit-code invariant

ASR completeness is determined from PAML's own declaration in `rst`, not from tip count. If the file says `Nodes X to Y are ancestral`, Pensieve requires exactly `Y-X+1` marginal sequence records, one each for nodes X..Y, from section `(1) Marginal reconstruction of ancestral sequences`. Joint-ASR records are ignored for core validation.

The smoke suite contains three explicit codeml/PAML fixtures:

1. A small `rst` says `Nodes 5 to 7 are ancestral`, contains all three marginal records, enters joint reconstruction, and codeml exits with status 1. Pensieve must validate the marginal section, emit a warning, record `CONTINUE_WITH_VALIDATED_MARGINAL_ASR`, and continue.
2. The same declared range is missing one marginal sequence and exits with status 1. Pensieve must reject it and fail.
3. A regression fixture matching the reported 103-species production case says `Nodes 104 to 205 are ancestral`. Pensieve must calculate 102 expected marginal sequences, accept 102 marginal node records, and ignore the additional 102 joint `node #` records (204 global `node #` lines total) for completeness decisions. No `n_tips-2` ASR-count guess is allowed.

The backend deletes stale `rst`/PAML output files before each new codeml attempt, and the resume path revalidates an existing `rst` before accepting it. The exact same marginal parser is used for validation and downstream integration.

## Test command

```bash
bash tests/smoke_test.sh
```

The biological suites from v3.30 plus the v3.31 portability suite passed in the build environment used to package v3.31.

## What the smoke suite actually tests

### CLI and static contracts

- `-h`, `-help`, `--help`, and `--help -long` all work.
- long help reports v3.31 and the true stage order.
- Python and shell syntax checks pass.
- `clock=0`, `fix_blength=0`, `cleandata=0`, `RateAncestor=1` are present in the codeml template.
- the authoritative runner contains no legacy `candidate_indel_frameshift_events.tsv` or `paml_indelmap_asr_combined` dependency.
- no nested codeml Slurm child-job registry remains.

### The v3.25 orchestration failure

A mock-MACSE/mock-codeml runner test executes the core workflow through integration:

```text
diagnostics -> alignment -> events -> asr -> integrate
```

and verifies that:

```text
results_02/T/02_T.codon_for_paml.phy
```

is created before codeml is invoked; the event layer is completed before ancestral integration; the PAML marginal scaffold is parsed; Pensieve structural states are overlaid; the lesion-aware ORF walk completes; and the integrated ancestral alignment is copied to `final_results`. This specifically covers the production error observed in v3.25:

```text
Missing or empty required file (step02_codon_for_paml)
```

### GUCA1B event structure

A 720-column synthetic alignment encodes:

- Miniopterus australis: 646-696 gap;
- M. natalensis: 646-687;
- M. schreibersii: 646-687;
- Nycteris thebaica: 661-669;
- residue-bearing outgroups.

The test requires:

- one confident shared deletion 646-687 affecting exactly the three Miniopterus;
- one australis-specific 688-696 extension;
- one independent Nycteris 661-669 event, not merged into/splitting the shared core.

### Exact parsimony ties

A synthetic 2-of-4 monophyletic gap pattern has equally parsimonious deletion and insertion histories. v3.31 must report representative event rows but requires:

```text
direction_confident = False
ambiguous_origin = True
biological_interpretation = ambiguous_indel_change
```

### Insertion/residue gain

A separate topology makes ancestral GAP + one residue gain strictly more parsimonious than multiple deletions. The test requires a confident `insertion_or_restoration` call, showing that native `! -> -` representation does not force deletion polarity.

### STOP identity and compensated frame restoration

TGA and TAA at the same aligned span are required to become two separate STOP characters. A separate regression case places MACSE `!` runs of lengths 1 and 2 upstream of a raw premature STOP; because the cumulative frame correction is `0 mod 3`, the STOP must remain an independent nonsense candidate rather than being discarded merely because upstream markers exist.

### Alignment preparation

The test verifies:

- a frameshift-only raw taxon (length not divisible by three, no premature STOP gate) still uses its MACSE `!` evidence;
- native `! -> -` and PAML-safe `! -> N` are synchronized;
- a raw premature TAA remains TAA in native coordinates and becomes NNN only in PAML view;
- the premature-stop registry contains canonical alignment coordinates;
- `02_T.codon_for_paml.phy` is produced;
- `--alignment defined` preserves every supplied sequence column and records that no realignment/column insertion occurred.

### PAML declaration, mapping and root policy

Mock PAML `rst` fixtures carry an explicit `Nodes X to Y are ancestral` declaration and all declared marginal sequences. The tests require:

- declared-range completeness is checked before topology mapping;
- all declared PAML marginal sequences are retained in an audit FASTA;
- degree-3 biological PAML internodes map to Pensieve nodes by topology;
- a TreeView serialization/root-only PAML vertex, when present, is kept for audit but excluded from biological-root inference;
- the biological root maps to PAML `NA` / `NO_PAML_SEQUENCE_EXPECTED`;
- the ORF walk reports the biological-root sequence status as `unavailable`, not fabricated.

### Compensatory restoration

A synthetic deep-tree case contains:

1. a confident 1-bp deletion on an internal branch;
2. a later confident 1-bp residue restoration on a descendant branch;
3. an apparently intact descendant ORF.

The required history is:

```text
first branch:      pseudogenization
descendant branch: apparent_orf_restoration
known history:     True
```

The compensatory event must not become a second pseudogenization and must not reset inherited history.


## Phylogenetic MSA ordering

A deliberately non-alphabetical rooted tree with terminal order `C, A, D, B` is used to test ordering. The suite requires Step 00, reused/alphabetized MACSE NT and AA outputs, and the final canonical alignment to be rewritten as `C, A, D, B` without changing alignment columns.

A second synthetic tree verifies the combined tip+internode MSA. With `((C,A)NodeCA,(D,B)NodeDB)UserRoot`, the expected FASTA order is `NodeCA, C, A, NodeDB, D, B`; the clock=0 root is present in `phylogenetic_sequence_order.tsv` but has no fabricated sequence.

## What was not executed in this build environment

Real external-program biological runs were not available here:

- real MACSE executable;
- real codeml executable;
- R/ggplot2 plotting runtime.

Accordingly:

- external command wrappers were syntax-checked;
- MACSE orchestration was exercised with a deterministic mock executable;
- PAML parsing/internode mapping was exercised with deterministic `rst` fixtures including the exact `Nodes 104 to 205` / 102-marginal / 204-total-`node #` regression pattern;
- an execution-level plotting smoke test (`tests/plot_smoke_test.sh`) is included and runs automatically when `Rscript` is available; in this packaging environment it is explicitly skipped because `Rscript` is not installed.

This smoke suite is deliberately stronger than v3.25's because it exercises the actual runner beyond the stage where v3.25 failed, through mock codeml ASR, PAML parsing, structural overlay and ORF-history integration, rather than testing only isolated event functions. It is still **not a substitute for a real production acceptance run** in the Pensieve conda environment.

## Recommended first production acceptance run

Before launching all genes, run at least:

1. GUCA1B, manually checking 646-687, 688-696 and Nycteris 661-669;
2. PDE6H, manually checking the known shared premature STOP pattern;
3. one `--alignment defined` gene to verify the curated alignment is retained exactly.

Any discrepancy should be treated as a regression and added to `tests/` before another package revision.

## v3.31 installation/runtime portability validation

`tests/install_portability_test.sh` checks that:

- installer help exposes `auto`, `mamba`, `conda`, `staged`, `venv`, and `current`;
- `venv` and `current` dry-runs use standard Python/pip and do not require conda/mamba;
- executable installation/runtime code contains no forced `ml miniforge`, `mamba activate`, or `conda activate`;
- local `--env-mode inherit` runs directly in the caller environment;
- Slurm `inherit` scripts contain no automatic activation/module command;
- explicit conda/mamba Slurm modes use non-interactive `<manager> run -n ENV`;
- venv mode sources only the user-selected virtual environment;
- `--slurm-module NAME` is opt-in and preserves the user-supplied site module name.

All files under `scripts/` are byte-identical to v3.30; v3.31 changes only installation/execution plumbing, version/help/docs, and tests.
