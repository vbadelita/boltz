# Architecture

`protein-folding` uses filesystem-backed queues and published static artifacts.

- Queue input: one deterministic JSON job spec per fold request.
- Queue state: files move across `pending/`, `running/`, `succeeded/`, and `failed/`.
- Dedupe index: `queue/by-hash/`.
- Run output: self-contained run directories under `experiments/<slug>/runs/<run_id>/`.
- Published output: immutable bundles under `published/<experiment_slug>/<run_id>/`.
- Viewer index: a merged flat CSV at `viewer/database.csv`.
