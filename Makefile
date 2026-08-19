SHELL := /bin/bash
ENV_NAME := Pensieve

.PHONY: install env env-mamba env-conda env-staged venv current test clean help

help:
	@echo "Pensieve Makefile"
	@echo "  make install            Portable auto install (mamba, then conda)"
	@echo "  make env                Create/update env only with auto backend"
	@echo "  make env-mamba          Create/update named environment with mamba"
	@echo "  make env-conda          Create/update named environment with conda"
	@echo "  make env-staged         Lower-solve-pressure staged mamba/conda installation"
	@echo "  make venv               Create PACKAGE/.venv and install Python dependencies"
	@echo "  make current            Install Python dependencies into current Python; no env"
	@echo "  make test               Run Pensieve smoke tests"
	@echo "  make clean              Remove temporary test files"

env:
	bash install.sh --env-only --backend=auto

env-mamba:
	bash install.sh --env-only --backend=mamba

env-conda:
	bash install.sh --env-only --backend=conda

env-staged:
	bash install.sh --env-only --backend=staged

venv:
	bash install.sh --env-only --backend=venv

current:
	bash install.sh --env-only --backend=current


install:
	bash install.sh

test:
	bash tests/smoke_test.sh

clean:
	rm -rf tests/tmp_smoke
