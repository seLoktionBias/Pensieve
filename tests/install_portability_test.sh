#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

TMP="tests/tmp_install_portability"
rm -rf "$TMP"
mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT

cat > "$TMP/T.fa" <<'FA'
>A
ATGAAATAG
>B
ATGAAATAG
FA
cat > "$TMP/tree.nwk" <<'TR'
(A:1,B:1);
TR

# Installer exposes portable backends and does not execute a site-specific module.
bash install.sh --help > "$TMP/install_help.txt"
grep -q -- '--backend=auto|mamba|conda|staged|venv|current' "$TMP/install_help.txt"
grep -q -- '--venv-path=PATH' "$TMP/install_help.txt"

bash install.sh --backend=venv --env-only --dry-run > "$TMP/venv_dry.txt"
grep -q ' -m venv ' "$TMP/venv_dry.txt"
grep -q 'requirements-pip.txt' "$TMP/venv_dry.txt"

bash install.sh --backend=current --env-only --dry-run > "$TMP/current_dry.txt"
grep -q ' -m pip install ' "$TMP/current_dry.txt"

# Executable installation/runtime code must not hard-code a specific HPC module
# or force shell activation of a conda/mamba environment.
if grep -R -n -E '^[[:space:]]*(ml[[:space:]]+miniforge|module[[:space:]]+load[[:space:]]+miniforge|mamba[[:space:]]+activate|conda[[:space:]]+activate)' \
    install.sh bin scripts Makefile 2>/dev/null; then
  echo "Hard-coded Miniforge/activation command remains in executable code" >&2
  exit 1
fi

COMMON=(python bin/pensieve --gene T --fasta "$TMP/T.fa" --tree "$TMP/tree.nwk" --workdir "$TMP/run" --dry-run)

# Local inherit: direct execution, no activation wrapper.
"${COMMON[@]}" --mode local --env-mode inherit > "$TMP/local_inherit.txt"
if grep -Eq 'mamba activate|conda activate|ml miniforge' "$TMP/local_inherit.txt"; then
  echo "Local inherit mode contains forced activation" >&2; exit 1
fi

# Slurm inherit: no module or environment activation at all.
"${COMMON[@]}" --mode slurm --env-mode inherit > "$TMP/slurm_inherit.txt"
grep -q 'Environment mode: inherit' "$TMP/slurm_inherit.txt"
if grep -Eq 'mamba activate|conda activate|ml miniforge|module load miniforge' "$TMP/slurm_inherit.txt"; then
  echo "Slurm inherit mode contains forced site-specific activation" >&2; exit 1
fi

# Explicit managers use non-interactive `run -n`, not activate.
"${COMMON[@]}" --mode slurm --env-mode conda --env-name Pensieve > "$TMP/slurm_conda.txt"
grep -q 'conda run -n Pensieve' "$TMP/slurm_conda.txt"
! grep -q 'conda activate' "$TMP/slurm_conda.txt"

"${COMMON[@]}" --mode slurm --env-mode mamba --env-name Pensieve > "$TMP/slurm_mamba.txt"
grep -q 'mamba run -n Pensieve' "$TMP/slurm_mamba.txt"
! grep -q 'mamba activate' "$TMP/slurm_mamba.txt"

# Venv mode uses only the user-selected venv path.
"${COMMON[@]}" --mode slurm --env-mode venv --venv-path "$TMP/myvenv" > "$TMP/slurm_venv.txt"
grep -q "source .*myvenv/bin/activate" "$TMP/slurm_venv.txt"

# Site module support is opt-in and preserves the user-supplied module name.
"${COMMON[@]}" --mode slurm --env-mode inherit --slurm-module custom_stack > "$TMP/slurm_module.txt"
grep -q 'module load custom_stack' "$TMP/slurm_module.txt"

# Long help documents all runtime modes.
python bin/pensieve --help -long > "$TMP/long_help.txt"
grep -q -- '--env-mode MODE' "$TMP/long_help.txt"
grep -q 'inherit|conda|mamba|venv' "$TMP/long_help.txt"
grep -q -- '--slurm-module NAME' "$TMP/long_help.txt"

# A single-shot `env create -f environment.yml` solve can be OOM-killed by
# the dependency solver, especially on HPC login nodes -- a real failure
# reported against this exact package, distinct from a real package
# conflict. auto/mamba/conda backends must retry with the staged installer
# automatically; --no-staged-fallback must disable that retry.
MOCKBIN="$TMP/mockbin"
mkdir -p "$MOCKBIN"
cat > "$MOCKBIN/mamba" <<'MOCKEOF'
#!/usr/bin/env bash
if [[ "$1" == "env" && "$2" == "list" ]]; then
  echo "# conda environments:"; echo "base   /fake/base"; exit 0
fi
if [[ "$1" == "env" && "$2" == "create" ]]; then
  echo "Solving environment: - environment: line 2: 99999 Killed" >&2
  exit 137
fi
if [[ "$1" == "create" || "$1" == "install" ]]; then
  echo "staged call ok: $*"; exit 0
fi
echo "UNHANDLED: $*" >&2; exit 1
MOCKEOF
chmod +x "$MOCKBIN/mamba"

PATH="$MOCKBIN:$PATH" bash install.sh --backend=mamba --env-name=T --env-only \
  > "$TMP/staged_fallback.txt" 2>&1
grep -q 'out-of-memory kill' "$TMP/staged_fallback.txt"
grep -q 'Using staged mamba installation' "$TMP/staged_fallback.txt"
grep -q 'Installation finished using backend: mamba' "$TMP/staged_fallback.txt"

if PATH="$MOCKBIN:$PATH" bash install.sh --backend=mamba --env-name=T --env-only --no-staged-fallback \
    > "$TMP/no_staged_fallback.txt" 2>&1; then
  echo "--no-staged-fallback should not silently recover from a killed solve" >&2
  exit 1
fi
if grep -q 'Using staged mamba installation' "$TMP/no_staged_fallback.txt"; then
  echo "--no-staged-fallback must not retry with the staged installer" >&2
  exit 1
fi

echo "Installation/runtime portability tests passed."
