# Repository Guidelines

## Project Structure & Module Organization
Core Python code lives in `src/boltz`, split between `data/` for parsing, tokenization, MSA, and dataset pipelines, and `model/` for layers, losses, optimizers, and model assemblies. Tests live under `tests/`, with layer-specific coverage in `tests/model/layers/` and broader regression or utility checks in files such as `tests/test_regression.py` and `tests/test_utils.py`. Training, evaluation, and preprocessing entry points live in `scripts/train/`, `scripts/eval/`, and `scripts/process/`. User-facing docs and example inputs are in `docs/` and `examples/`.

## Build, Test, and Development Commands
Install a local dev environment with `pip install -e .[test,lint]` and add `[cuda]` on supported GPUs. Run inference locally with `boltz predict examples/prot.yaml --help` or `boltz predict <input_path> --use_msa_server`. Run tests with `pytest`; skip long-running cases with `pytest -m "not slow"`, or target a file, for example `pytest tests/test_utils.py`. Lint the codebase with `ruff check src tests scripts` and format with `ruff format`.

To serve the results viewer, start an HTTP server from the repository root, not from `results_viewer/`, because `results_viewer/database.csv` contains structure paths like `test_eve_experiments/...` and `runs/...` that must resolve from the repo root. Use `python3 -m http.server 8766 --bind 127.0.0.1 --directory .` and open `http://127.0.0.1:8766/results_viewer/`.

## Coding Style & Naming Conventions
Use 4-space indentation, type hints on public interfaces, and NumPy-style docstrings when documentation is needed. Follow existing Python naming: `snake_case` for functions and modules, `PascalCase` for classes, and descriptive config names such as `structure.yaml` or `confidence.yaml`. Ruff is configured as the primary style gate in `pyproject.toml`; keep imports clean and avoid introducing patterns already ignored by the project only when justified.

## Testing Guidelines
Pytest is the test runner. Name new tests `test_<behavior>.py`, and keep unit coverage close to the affected module when practical, for example `tests/model/layers/` for new layer logic. Use the existing markers `slow` and `regression` consistently so selective runs remain reliable. Add regression tests for bug fixes and run the narrowest relevant `pytest` target before opening a PR.

## Commit & Pull Request Guidelines
Recent history uses short, imperative subjects such as `bump version`, `Update prediction.md...`, and `Minor tweaks to the documentation`. Keep commit titles concise and action-oriented; separate unrelated changes into different commits. PRs should explain the user-visible or training/evaluation impact, link issues when applicable, and include command evidence for validation. Attach screenshots only for documentation or output-format changes.

## Configuration & Data Notes
Python `>=3.10,<3.13` is required. Large training and preprocessing workflows depend on external datasets and tools documented in `docs/training.md` and `scripts/process/README.md`; avoid hardcoding local paths or credentials in committed configs.
