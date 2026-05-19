from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
QUEUE_ROOT = REPO_ROOT / "queue"
PUBLISHED_ROOT = REPO_ROOT / "published"
VIEWER_DB_PATH = REPO_ROOT / "viewer" / "database.csv"
HOSTS_CONFIG_PATH = REPO_ROOT / "configs" / "hosts.toml"
RUNTIME_CONFIG_PATH = REPO_ROOT / "configs" / "runtime.toml"

QUEUE_DIRS = {
    "pending": QUEUE_ROOT / "pending",
    "running": QUEUE_ROOT / "running",
    "succeeded": QUEUE_ROOT / "succeeded",
    "failed": QUEUE_ROOT / "failed",
    "by_hash": QUEUE_ROOT / "by-hash",
}


def ensure_repo_layout() -> None:
    """Create the expected top-level directories if they do not exist."""

    for path in [*QUEUE_DIRS.values(), PUBLISHED_ROOT, VIEWER_DB_PATH.parent]:
        path.mkdir(parents=True, exist_ok=True)


def repo_relative(path: Path) -> str:
    """Return a portable repo-relative path."""

    return path.relative_to(REPO_ROOT).as_posix()
