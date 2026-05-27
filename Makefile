# OPMD MSI biomarker pipeline — common targets.

PY ?= python
PIP ?= $(PY) -m pip
PYTEST ?= $(PY) -m pytest

.PHONY: help install install-fetch demo demo-clean test test-cov fetch harmonize contrasts filter train rank report lint typecheck clean

help:
	@echo "Targets:"
	@echo "  install        Install package (core deps only)."
	@echo "  install-fetch  Install + heavy GEO/TCGA fetcher extras."
	@echo "  demo           Run end-to-end synthetic demo (no network)."
	@echo "  test           Run unit + integration tests."
	@echo "  test-cov       Run tests with coverage report."
	@echo "  fetch          Fetch real GEO/TCGA datasets per manifests/datasets.yaml."
	@echo "  harmonize      Cross-platform normalization."
	@echo "  contrasts      Compute 3 harmonized contrasts."
	@echo "  filter         Apply surface-localization + baseline-normal filters."
	@echo "  train          Fit L1 Cox + L1 logistic with leave-one-cohort-out CV."
	@echo "  rank           Compute SHAP-based composite ranking."
	@echo "  report         Render reviewer-ready ranking + per-marker reports."
	@echo "  clean          Remove caches and demo outputs."

install:
	$(PIP) install -e .

install-fetch:
	$(PIP) install -e ".[fetch,dev]"

demo:
	$(PY) -m opmd_msi.cli demo --out outputs/demo

demo-clean:
	rm -rf outputs/demo

test:
	$(PYTEST)

test-cov:
	$(PYTEST) --cov=opmd_msi --cov-report=term-missing

fetch:
	$(PY) scripts/01_fetch_data.py --manifest manifests/datasets.yaml --out data/raw

harmonize:
	$(PY) scripts/02_harmonize.py --raw data/raw --out data/harmonized

contrasts:
	$(PY) scripts/03_contrasts.py --harmonized data/harmonized --out outputs/contrasts

filter:
	$(PY) scripts/04_filter.py --contrasts outputs/contrasts --out outputs/filtered

train:
	$(PY) scripts/05_train_models.py --filtered outputs/filtered --harmonized data/harmonized --out outputs/models

rank:
	$(PY) scripts/06_rank.py --models outputs/models --contrasts outputs/contrasts --out outputs/ranking

report:
	$(PY) -m opmd_msi.cli report --ranking outputs/ranking --out outputs/report

lint:
	ruff check src tests

typecheck:
	mypy src

clean:
	rm -rf outputs/demo .pytest_cache .ruff_cache .mypy_cache build dist *.egg-info
	find . -type d -name __pycache__ -exec rm -rf {} +
