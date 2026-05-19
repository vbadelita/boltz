#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from utils.common import HOSTS_CONFIG_PATH, PUBLISHED_ROOT, REPO_ROOT, repo_relative

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover
    import tomli as tomllib  # type: ignore[no-redef]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Publish a completed run bundle.")
    parser.add_argument("--run-dir", type=Path, required=True, help="Completed run directory.")
    parser.add_argument(
        "--dest-root",
        type=Path,
        help="Local publish root override. Defaults to the configured published root.",
    )
    parser.add_argument(
        "--ssh-target",
        help="Optional SSH target for rsync publish, for example user@viewer-host.",
    )
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def load_hosts_config() -> dict[str, Any]:
    if not HOSTS_CONFIG_PATH.exists():
        return {}
    return tomllib.loads(HOSTS_CONFIG_PATH.read_text())


def publish_locally(run_dir: Path, destination: Path) -> None:
    if destination.exists():
        shutil.rmtree(destination)
    shutil.copytree(run_dir, destination)


def publish_via_rsync(run_dir: Path, destination: Path, ssh_target: str, dry_run: bool) -> None:
    command = ["rsync", "-av", f"{run_dir}/", f"{ssh_target}:{destination.as_posix()}/"]
    if dry_run:
        print(" ".join(command))
        return
    subprocess.run(command, check=True)  # noqa: S603


def main() -> int:
    args = parse_args()
    hosts_config = load_hosts_config()
    run_dir = args.run_dir.resolve()
    experiment_slug = run_dir.parents[1].name
    run_id = run_dir.name
    dest_root = (
        args.dest_root.resolve()
        if args.dest_root
        else (REPO_ROOT / hosts_config.get("viewer", {}).get("publish_root", "published"))
    )
    destination = dest_root / experiment_slug / run_id
    destination.parent.mkdir(parents=True, exist_ok=True)

    if args.ssh_target:
        publish_via_rsync(run_dir, destination, args.ssh_target, args.dry_run)
    elif not args.dry_run:
        publish_locally(run_dir, destination)

    print(f"Published {repo_relative(run_dir)} -> {destination.as_posix()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
