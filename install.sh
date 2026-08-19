#!/usr/bin/env bash
set -uo pipefail

# Pensieve portable installer.
#
# No HPC module command is assumed. In particular, Pensieve never executes
# `ml miniforge` or `module load ...` during installation. If a site requires a
# module to expose conda/mamba, the user should load that module before running
# this installer.
#
# Installation backends:
#   auto     try mamba, then conda; if neither exists, explain venv/current
#   mamba    create/update named conda-style environment with mamba
#   conda    create/update named conda-style environment with conda
#   staged   lower-solve-pressure staged install; mamba first, then conda
#   venv     create a Python virtual environment and install Python dependencies
#   current  install Python dependencies into the current Python environment
#
# venv/current deliberately do NOT pretend to install non-Python tools. MACSE,
# codeml/PAML and R must already be available on PATH (or be provided by the
# user's system/module/container) for a full run.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

ENV_NAME="${PENSIEVE_ENV_NAME:-Pensieve}"
ENV_FILE="$ROOT_DIR/environment.yml"
PIP_FILE="$ROOT_DIR/requirements-pip.txt"
VENV_PATH="${PENSIEVE_VENV_PATH:-$ROOT_DIR/.venv}"
BACKEND="${PENSIEVE_INSTALL_BACKEND:-auto}"
ENV_ONLY=0
DRY_RUN=0
NO_STAGED_FALLBACK="${PENSIEVE_NO_STAGED_FALLBACK:-0}"

usage(){
  cat <<'USAGE'
Usage: bash install.sh [options]

Options:
  --backend=auto|mamba|conda|staged|venv|current
                              Installation method (default: auto)
  --env-name=NAME             Conda/mamba environment name (default: Pensieve)
  --venv-path=PATH            Python venv path (default: PACKAGE/.venv)
  --env-only                  Deprecated no-op (retained for compatibility)
  --no-staged-fallback        Do not automatically retry a failed/killed
                              mamba or conda solve with the lower-memory
                              staged installer (default: retry automatically;
                              see HPC note below)
  --dry-run                   Print planned commands without executing them
  -h, --help                  Show this help

Environment variables:
  PENSIEVE_INSTALL_BACKEND    Same values as --backend
  PENSIEVE_ENV_NAME           Conda/mamba environment name
  PENSIEVE_VENV_PATH          Python virtual-environment path
  PENSIEVE_NO_STAGED_FALLBACK Same as --no-staged-fallback

Examples:
  # mamba or conda, whichever is available
  bash install.sh

  # explicitly use conda
  bash install.sh --backend=conda

  # Python virtual environment; external tools must be available separately
  bash install.sh --backend=venv
  source .venv/bin/activate

  # no new environment: install Python packages into the active/current Python
  bash install.sh --backend=current

HPC note:
  If your cluster requires a module to expose conda/mamba, load it yourself
  before installation, e.g. `module load <your-site-module>`. Pensieve does not
  hard-code any site-specific module such as Miniforge.

  A single-shot `env create -f environment.yml` solve (~20 conda-forge/
  bioconda packages including a full R stack) can need more memory than an
  HPC login node allows, and fails as a plain "Killed" process, not a
  package-conflict error. auto/mamba/conda backends automatically retry with
  the staged installer (one dependency group per solve, much lower peak
  memory) when that happens; pass --backend=staged directly to skip the
  first attempt, or --no-staged-fallback to disable the retry entirely.
USAGE
}

for arg in "$@"; do
  case "$arg" in
    --env-only) ENV_ONLY=1 ;;
    --dry-run) DRY_RUN=1 ;;
    --no-staged-fallback) NO_STAGED_FALLBACK=1 ;;
    --backend=*) BACKEND="${arg#--backend=}" ;;
    --env-name=*) ENV_NAME="${arg#--env-name=}" ;;
    --venv-path=*) VENV_PATH="${arg#--venv-path=}" ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $arg" >&2; usage >&2; exit 1 ;;
  esac
done

have_cmd(){ command -v "$1" >/dev/null 2>&1; }

run(){
  printf '[Pensieve]'
  printf ' %q' "$@"
  printf '\n'
  [[ "$DRY_RUN" -eq 1 ]] || "$@"
}

env_exists_with(){
  local manager="$1"
  "$manager" env list 2>/dev/null | awk '{print $1}' | grep -qx "$ENV_NAME"
}

print_external_tool_note(){
  cat <<'NOTE'

[Pensieve] Python environment setup is complete.
For a full Pensieve run, these non-Python programs must also be available on PATH:
  macse     frameshift-aware coding alignment
  codeml    PAML ancestral substitution reconstruction
  Rscript   plotting (with ape, ggplot2, dplyr, readr, stringr, tidyr, optparse)
NOTE
}

try_mamba_env(){
  have_cmd mamba || return 127
  echo "[Pensieve] Using mamba environment backend"
  if env_exists_with mamba; then
    run mamba env update -n "$ENV_NAME" -f "$ENV_FILE"
  else
    run mamba env create -n "$ENV_NAME" -f "$ENV_FILE"
  fi
}

try_conda_env(){
  have_cmd conda || return 127
  echo "[Pensieve] Using conda environment backend"
  if env_exists_with conda; then
    run conda env update -n "$ENV_NAME" -f "$ENV_FILE"
  else
    run conda env create -n "$ENV_NAME" -f "$ENV_FILE"
  fi
}

try_staged_mamba(){
  have_cmd mamba || return 127
  echo "[Pensieve] Using staged mamba installation"
  if ! env_exists_with mamba; then
    run mamba create -y -n "$ENV_NAME" -c conda-forge python=3.11.15 pip git
  fi
  run mamba install -y -n "$ENV_NAME" -c conda-forge numpy scipy pandas biopython ete3
  run mamba install -y -n "$ENV_NAME" -c conda-forge r-base r-ape r-ggplot2 r-dplyr r-readr r-stringr r-tidyr r-optparse
  run mamba install -y -n "$ENV_NAME" -c conda-forge -c bioconda iqtree paml macse emboss
}

try_staged_conda(){
  have_cmd conda || return 127
  echo "[Pensieve] Using staged conda installation"
  if ! env_exists_with conda; then
    run conda create -y -n "$ENV_NAME" -c conda-forge python=3.11.15 pip git
  fi
  run conda install -y -n "$ENV_NAME" -c conda-forge numpy scipy pandas biopython ete3
  run conda install -y -n "$ENV_NAME" -c conda-forge r-base r-ape r-ggplot2 r-dplyr r-readr r-stringr r-tidyr r-optparse
  run conda install -y -n "$ENV_NAME" -c conda-forge -c bioconda iqtree paml macse emboss
}

try_staged(){
  try_staged_mamba && return 0
  local status=$?
  [[ "$status" -eq 127 ]] || echo "[Pensieve] staged mamba failed with status $status; trying conda" >&2
  try_staged_conda
}

# A single `env create -f environment.yml` solve (~20 conda-forge/bioconda
# packages including a full R stack) can need more memory than the solver has
# to work with, especially on HPC login nodes with tight per-process memory
# caps. That failure mode looks like "Killed" with no package-conflict
# message, not a real dependency error, so it is worth automatically retrying
# with the staged/incremental installer (each solve only adds a small group
# of packages to an already-mostly-solved environment, so peak memory is much
# lower) rather than making the user discover --backend=staged themselves.
mamba_with_staged_fallback(){
  try_mamba_env; local status=$?
  [[ "$status" -eq 0 || "$status" -eq 127 || "$NO_STAGED_FALLBACK" -eq 1 ]] && return "$status"
  echo "[Pensieve] mamba env create/update failed or was killed (status=$status)." >&2
  echo "[Pensieve] This is commonly an out-of-memory kill from the dependency solver, not a real package conflict; retrying with a lower-memory staged install." >&2
  try_staged_mamba
}

conda_with_staged_fallback(){
  try_conda_env; local status=$?
  [[ "$status" -eq 0 || "$status" -eq 127 || "$NO_STAGED_FALLBACK" -eq 1 ]] && return "$status"
  echo "[Pensieve] conda env create/update failed or was killed (status=$status)." >&2
  echo "[Pensieve] This is commonly an out-of-memory kill from the dependency solver, not a real package conflict; retrying with a lower-memory staged install." >&2
  try_staged_conda
}

try_venv(){
  have_cmd python3 || have_cmd python || { echo "[ERROR] python3/python not found" >&2; return 127; }
  local py="python3"
  have_cmd python3 || py="python"
  echo "[Pensieve] Creating/updating Python venv: $VENV_PATH"
  if [[ ! -x "$VENV_PATH/bin/python" ]]; then
    run "$py" -m venv "$VENV_PATH" || return $?
  fi
  run "$VENV_PATH/bin/python" -m pip install --upgrade pip || return $?
  run "$VENV_PATH/bin/python" -m pip install -r "$PIP_FILE" || return $?
  print_external_tool_note
}

try_current(){
  local py=""
  if have_cmd python3; then py="python3"; elif have_cmd python; then py="python"; else
    echo "[ERROR] python3/python not found" >&2; return 127
  fi
  echo "[Pensieve] Installing Python dependencies into current Python: $(command -v "$py")"
  run "$py" -m pip install -r "$PIP_FILE" || return $?
  print_external_tool_note
}

install_env(){
  local status=0
  case "$BACKEND" in
    auto)
      if have_cmd mamba; then
        mamba_with_staged_fallback; status=$?
        [[ "$status" -eq 0 ]] && return 0
        echo "[Pensieve] mamba (including the staged fallback) failed with status $status; trying conda" >&2
      fi
      if have_cmd conda; then
        conda_with_staged_fallback; status=$?
        [[ "$status" -eq 0 ]] && return 0
        echo "[Pensieve] conda (including the staged fallback) failed with status $status" >&2
      fi
      cat >&2 <<EOF2
[ERROR] No usable mamba/conda installation succeeded, including the
lower-memory staged fallback.
A solve that just gets "Killed" with no package-conflict message is almost
always the dependency solver being killed for using too much memory -- common
on HPC login nodes, which are often capped well below what solving this
environment's ~20 conda-forge/bioconda packages (including a full R stack)
can need. If you are on a login node, try an interactive compute allocation
instead, e.g.:
  srun --mem=8G --time=00:30:00 --pty bash
  ml <your-site-module>
  bash install.sh --backend=staged
'conda config --set channel_priority strict' often reduces solver memory/time
substantially too. Otherwise choose one of:
  bash install.sh --backend=venv
  bash install.sh --backend=current
or load your site's conda/mamba module yourself and rerun this installer.
EOF2
      return 127
      ;;
    mamba) mamba_with_staged_fallback ;;
    conda) conda_with_staged_fallback ;;
    staged) try_staged ;;
    venv) try_venv ;;
    current) try_current ;;
    *) echo "[ERROR] Invalid backend '$BACKEND'. Use auto, mamba, conda, staged, venv, or current." >&2; return 2 ;;
  esac
}

install_env || exit $?

cat <<EOF2

[Pensieve] Installation finished using backend: $BACKEND
EOF2
case "$BACKEND" in
  mamba|conda|staged|auto)
    cat <<EOF2
Environment name: $ENV_NAME
For interactive use, activate that environment using your normal conda/mamba
workflow if desired. Pensieve's generated Slurm scripts do not require shell
activation when '--env-mode conda' or '--env-mode mamba' is selected; those
modes use the manager's non-interactive 'run -n' command.
EOF2
    ;;
  venv)
    echo "Activate with: source $VENV_PATH/bin/activate"
    ;;
  current)
    echo "Pensieve will use the current environment/PATH. No environment was created."
    ;;
esac

echo "Pensieve v$(cat "$ROOT_DIR/VERSION" 2>/dev/null || echo '?') installed from $ROOT_DIR."
