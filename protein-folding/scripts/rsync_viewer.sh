#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 <user@host:/remote/viewer/root> [extra-rsync-arg]" >&2
  echo "Example: $0 ubuntu@viewer:/srv/protein-folding-viewer" >&2
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SOURCE_DIR="${REPO_ROOT}/viewer/"
DEST="$1"

if [[ ! -d "${SOURCE_DIR}" ]]; then
  echo "Viewer directory not found: ${SOURCE_DIR}" >&2
  exit 1
fi

EXTRA_ARGS=()
if [[ $# -eq 2 ]]; then
  EXTRA_ARGS+=("$2")
fi

rsync -avzP --delete "${EXTRA_ARGS[@]}" "${SOURCE_DIR}" "${DEST}/"
