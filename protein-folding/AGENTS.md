# protein-folding Working Notes

## Purpose

`protein-folding/` is a new orchestration and deployment workspace that separates:

- experiment preparation
- queueing and run execution
- viewer publishing
- standalone web deployment

from the upstream runtime repositories:

- `boltz/`
- `lmi4boltz/`

The long-term goal is:

1. prepare experiment-specific inputs
2. emit queueable JSON job specs and Boltz YAMLs
3. run jobs with a shared stateless runner
4. publish viewer-ready artifacts
5. deploy only the standalone `viewer/` folder to a web host

## Current State

### Implemented

- Filesystem-backed queue scaffold:
  - `queue/pending/`
  - `queue/running/`
  - `queue/succeeded/`
  - `queue/failed/`
  - `queue/by-hash/`
- Deterministic job hashing based on canonical JSON plus YAML content digest
- Shared queue runner:
  - `scripts/run_queue.py`
- Shared publish step:
  - `scripts/publish_run.py`
- Legacy importers:
  - `scripts/import_legacy_runs.py`
  - `scripts/import_results_database.py`
- Standalone viewer DB merge:
  - `scripts/merge_viewer_db.py`
- Standalone viewer server:
  - `scripts/serve_viewer.py`
- Standalone deploy helpers:
  - `scripts/rsync_viewer.sh`
  - `scripts/install_viewer_systemd.sh`
- Standalone web root:
  - `viewer/index.html`
  - `viewer/database.csv`
  - `viewer/published/...`
  - `viewer/experiments/<experiment_slug>/summary.md`

### Important design decision

The deployable website artifact is now only:

- `protein-folding/viewer/`

That folder contains everything the browser needs. The server does not need the rest of the repo to serve the site.

## Viewer Design

The current `viewer/index.html` intentionally preserves the old `results_viewer` look and feel.

Changes made relative to the original viewer:

- reads `./database.csv` instead of the old shared path
- serves from the viewer root directly, so `http://host:port/` opens the site
- supports deep links with query params:
  - `?e=<experiment_slug>`
  - `?e=<experiment_slug>&t=<target_slug>`
  - `?e=<experiment_slug>&t=<target_slug>&s=<config|sample_rank>`
- updates the URL as selection changes
- has an `Experiment Notes` button that opens a modal with markdown notes

Recent bug fixes already applied:

- experiment dropdown now resets target selection correctly instead of leaking a stale target from another experiment
- target dropdown now shows only the target label, not `folder / target`
- target dropdown explicitly selects the first valid option when the previous selection no longer applies
- `serve_viewer.py` uses a reusable TCP server so Ctrl-C restarts are less likely to fail with `Address already in use`

## Viewer Data Layout

Everything served by the standalone site lives under `viewer/`:

- `viewer/index.html`
- `viewer/database.csv`
- `viewer/experiments/<experiment_slug>/summary.md`
- `viewer/published/<experiment_slug>/<bundle_slug>/...`

Artifact paths in `viewer/database.csv` are relative to `viewer/`, not the repo root.

This is important: it allows rsyncing only `viewer/` to the remote host.

## Notes / Markdown Export

Experiment descriptions are exported during `merge_viewer_db.py`.

Rules:

- if a curated source markdown exists, copy it through
- otherwise generate a compact markdown summary automatically from merged metadata

Curated sources currently wired:

- `flu_trimer_ligand_experiments/summary.md`
- `flu_trimer_single_ligand_affinity/summary.md`
- `test_eve_experiments/summary.md`
- `test_eve_experiments/multimers/summary.md`

Generated notes are produced for experiments that do not have a curated source file.

The exporter lives in:

- `scripts/utils/experiment_notes.py`

## Imported Data Status

The standalone viewer currently contains imported legacy data from:

- Emma EVE gag runs
- EVE structures
- EVE multimers
- flu trimer ligand screen
- flu trimer single-ligand affinity
- influenza H5 affinity matrix
- influenza MSA ladder
- influenza SIA/NGC ladder
- influenza glycan linkage
- `trix_s40_n6_p2_r5`

At the time of writing, `viewer/database.csv` contains 1222 rows after merge.

The actual generated `viewer/database.csv` is ignored and should not be committed.

## Ignore / Commit Policy

Generated deployment data should stay out of git:

- `viewer/database.csv`
- everything under `viewer/published/`

Tracked:

- code
- configs
- tests
- templates
- docs
- `.gitkeep` placeholders

## Deployment Workflow

### Local refresh

Run:

```bash
python3 protein-folding/scripts/merge_viewer_db.py
```

This rebuilds:

- `viewer/database.csv`
- `viewer/experiments/<experiment_slug>/summary.md`

### Ship to remote host

Run:

```bash
./protein-folding/scripts/rsync_viewer.sh user@host:/srv/protein-folding-viewer
```

This syncs only the standalone viewer web root.

### Install server on remote host

Run:

```bash
sudo ./protein-folding/scripts/install_viewer_systemd.sh /srv/protein-folding-viewer
```

Installed files:

- `/etc/systemd/system/protein-folding-viewer.service`
- `/etc/default/protein-folding-viewer`

The service serves the viewer root directly via:

- `python3 -m http.server`

Default URL:

- `http://host:8766/`

## Queue / Runner Model

Queue jobs are JSON files with:

- experiment metadata
- target metadata
- `yaml_path`
- program name
- structured args
- optional env overrides
- deterministic `job_hash`

Job hashing is implemented in:

- `scripts/utils/job_spec.py`

The stateless runner:

- claims jobs from `queue/pending/`
- validates the hash
- checks `queue/by-hash/`
- moves work to `queue/running/`
- executes `boltz` or `lmi4boltz`
- writes run-local outputs
- moves job to `succeeded/` or `failed/`

Current support:

- `program = "boltz"`
- `program = "lmi4boltz"`

## Tests

Current focused tests live under:

- `tests/test_job_spec.py`
- `tests/test_run_queue.py`
- `tests/test_merge_viewer_db.py`
- `tests/test_import_legacy_runs.py`
- `tests/test_import_results_database.py`
- `tests/test_experiment_notes.py`

Expected check:

```bash
pytest protein-folding/tests -q
```

## Known Gaps / Next Work

These pieces are not fully finished yet:

1. Experiment preparation is still only scaffolded.
   - We do not yet have a full queue-emitting preparer for the Emma EVE gag workflows.

2. Legacy run import is implemented, but the long-term path should be:
   - prepare new jobs in `protein-folding`
   - run via `run_queue.py`
   - publish directly into `viewer/published/`
   - avoid depending on one-off legacy imports

3. The old top-level `published/` tree still exists as a historical artifact from an earlier phase.
   - The active deployment path is `viewer/published/`

4. The viewer markdown renderer is intentionally lightweight.
   - It supports headings, bullet lists, inline code, fenced code blocks, and simple pipe tables
   - It is not a full Markdown engine

5. Some experiment notes are auto-generated because no curated markdown exists yet.
   - If better summaries are wanted, add real source `summary.md` files and wire them in `SOURCE_SUMMARIES`

## Most Relevant Files

- `protein-folding/scripts/run_queue.py`
- `protein-folding/scripts/merge_viewer_db.py`
- `protein-folding/scripts/import_legacy_runs.py`
- `protein-folding/scripts/import_results_database.py`
- `protein-folding/scripts/utils/job_spec.py`
- `protein-folding/scripts/utils/viewer_db.py`
- `protein-folding/scripts/utils/experiment_notes.py`
- `protein-folding/viewer/index.html`
- `protein-folding/deploy/systemd/protein-folding-viewer.service`
- `protein-folding/scripts/rsync_viewer.sh`

## Recent Commits

Useful checkpoints in the parent repo history:

- `8d9f418` `Add protein-folding scaffold`
- `6879517` `Refine protein-folding viewer workflow`

## Working Assumptions

- The compute side may continue to live in the current `boltz` workspace for now.
- The remote viewer host only needs the standalone `viewer/` folder.
- Viewer data is regenerated locally and then synced, not edited manually on the server.
- We prefer preserving the old viewer UX and extending it carefully instead of redesigning it.
