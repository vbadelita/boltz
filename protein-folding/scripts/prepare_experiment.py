#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from utils.common import QUEUE_DIRS, REPO_ROOT, ensure_repo_layout
from utils.job_spec import dump_json, ensure_job_hash, load_json


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Enqueue prepared job specs into protein-folding/queue/pending."
    )
    parser.add_argument(
        "--job",
        action="append",
        default=[],
        help="Path to a job JSON file. Repeatable.",
    )
    parser.add_argument(
        "--jobs-dir",
        type=Path,
        help="Directory containing one or more job JSON files.",
    )
    parser.add_argument(
        "--queue-dir",
        type=Path,
        default=QUEUE_DIRS["pending"],
        help="Pending queue directory.",
    )
    return parser.parse_args()


def collect_job_paths(args: argparse.Namespace) -> list[Path]:
    paths = [Path(item).resolve() for item in args.job]
    if args.jobs_dir:
        paths.extend(sorted(args.jobs_dir.resolve().glob("*.json")))
    return paths


def main() -> int:
    args = parse_args()
    ensure_repo_layout()
    job_paths = collect_job_paths(args)
    if not job_paths:
        print("No job JSON files supplied.", file=sys.stderr)
        return 1

    queue_dir = args.queue_dir.resolve()
    queue_dir.mkdir(parents=True, exist_ok=True)

    enqueued = 0
    last_destination = queue_dir
    for path in job_paths:
        job = ensure_job_hash(load_json(path))
        destination = queue_dir / f"{job['job_hash'].split(':', 1)[1]}.json"
        last_destination = destination.parent
        if destination.exists():
            continue
        dump_json(destination, job)
        enqueued += 1

    print(f"Queued {enqueued} job(s) into {last_destination.relative_to(REPO_ROOT).as_posix()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
