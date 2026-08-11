This directory is a fallback only. As of v3.33, `bin/pensieve` auto-detects
PAML's own dat/ directory (empirical amino-acid substitution matrices such as
dayhoff.dat/wag.dat/grantham.dat, installed by the bioconda `paml` package
directly under `<conda/mamba env prefix>/dat`) from wherever `codeml` actually
resolves for the chosen --env-mode, and passes that real path via --dat-dir.
This placeholder is used only if codeml cannot be resolved that way, or if
--dat-dir is not overridden explicitly.

Pensieve's own codon-model control file (templates/dummy_codon_asr.ctl) does
not currently read anything from this directory -- it has no aaRatefile /
empirical amino-acid model setting -- so nothing here is required for a run to
succeed today. --dat-dir exists for forward/CLI compatibility.
