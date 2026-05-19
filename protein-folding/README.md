# protein-folding

`protein-folding` is a lightweight orchestration repo for protein-folding experiments.
It keeps experiment preparation, queueing, run metadata, publishing, and the static
results viewer separate from the upstream `boltz` and `lmi4boltz` source trees.

## Layout

```text
protein-folding/
  experiments/
  queue/
  published/
  viewer/
  scripts/
  configs/
  state/
  docs/
```

## Core workflow

1. Prepare experiment-specific inputs and emit Boltz YAMLs plus queue job JSON files.
2. Run the stateless queue runner on the compute host.
3. Publish completed run bundles into `viewer/published/` or to a viewer host.
4. Merge published viewer metadata into `viewer/database.csv`.
5. Serve the `viewer/` directory as a standalone static site.

## Commands

From the `protein-folding` repo root:

```bash
python3 scripts/prepare_experiment.py --help
python3 scripts/run_queue.py --help
python3 scripts/publish_run.py --help
python3 scripts/merge_viewer_db.py --help
python3 scripts/serve_viewer.py --help
```

## Host model

- Compute host: runs `boltz` or `lmi4boltz` and writes experiment run directories.
- Viewer host: receives published bundles under `viewer/published/`, rebuilds
  `viewer/database.csv`, and serves the `viewer/` directory.

## Runtime layout

The runtime repos are intentionally not tracked here:

- `boltz/`
- `lmi4boltz/`

Point the runner at them with config values or environment overrides.
