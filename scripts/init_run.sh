#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <run-name> <command...>" >&2
  exit 1
fi

run_name="$1"
shift

run_dir="runs/${run_name}"
mkdir -p "${run_dir}"

cat > "${run_dir}/repro.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail

# Reproduce run '${run_name}' from the repository root.

$*
EOF

chmod +x "${run_dir}/repro.sh"
echo "Created ${run_dir}/repro.sh"
