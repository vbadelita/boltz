#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"

default_host="${RESULTS_VIEWER_HOST:-localhost}"
default_port="${RESULTS_VIEWER_PORT:-9000}"

host="${1:-${default_host}}"
port="${2:-${default_port}}"
viewer_url="http://${host}:${port}/results_viewer/"

cd "${repo_root}"

python3 scripts/process/combine_results_viewer_db.py

echo "Serving repository root from ${repo_root}"
echo "Viewer URL: ${viewer_url}"

python3 -m http.server "${port}" --bind "${host}" --directory "${repo_root}"
