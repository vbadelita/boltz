#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 <viewer-root> [service-name]" >&2
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
DEPLOY_DIR="${REPO_ROOT}/deploy/systemd"
VIEWER_ROOT="$1"
SERVICE_NAME="${2:-protein-folding-viewer}"
UNIT_PATH="/etc/systemd/system/${SERVICE_NAME}.service"
ENV_PATH="/etc/default/${SERVICE_NAME}"

install -d "${VIEWER_ROOT}"
install -m 0644 "${DEPLOY_DIR}/protein-folding-viewer.service" "${UNIT_PATH}"
install -m 0644 "${DEPLOY_DIR}/protein-folding-viewer.env" "${ENV_PATH}"

sed -i "s#^VIEWER_ROOT=.*#VIEWER_ROOT=${VIEWER_ROOT}#" "${ENV_PATH}"

systemctl daemon-reload
systemctl enable --now "${SERVICE_NAME}.service"

echo "Installed ${SERVICE_NAME}.service"
echo "Viewer root: ${VIEWER_ROOT}"
echo "Unit: ${UNIT_PATH}"
echo "Env: ${ENV_PATH}"
