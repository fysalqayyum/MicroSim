#!/usr/bin/env bash

set -euo pipefail

base_ref="${1:-upstream/main}"

if ! git rev-parse --verify "${base_ref}^{commit}" >/dev/null 2>&1; then
  echo "Unknown base ref: ${base_ref}" >&2
  echo "Usage: $0 [base-ref]" >&2
  exit 2
fi

patch_file="$(mktemp)"
trap 'rm -f "${patch_file}"' EXIT

git diff --binary "${base_ref}" -- . \
  ':!scripts/pre-publish-check.sh' > "${patch_file}"

if grep -Eiq \
  'BEGIN [A-Z ]*PRIVATE KEY|github_pat_[A-Za-z0-9_]+|ghp_[A-Za-z0-9]+|AKIA[0-9A-Z]{16}|(api[_-]?key|access[_-]?token|client[_-]?secret|password)[[:space:]]*[:=][[:space:]]*[^[:space:]]+' \
  "${patch_file}"; then
  echo "Potential credential detected in outgoing changes." >&2
  exit 1
fi

if grep -Eiq \
  '^\+.*(/Users/[^/[:space:]]+|/home/[^/[:space:]]+|/scratch/[^/[:space:]]+|tabuk-gpu|utfqayyum)' \
  "${patch_file}"; then
  echo "Machine-specific path, host, or account detected in outgoing changes." >&2
  exit 1
fi

bad_name=0
large_file=0
while IFS= read -r -d '' file; do
  [[ -e "${file}" ]] || continue
  case "${file}" in
    *.pem|*.key|*.p12|*.pfx|*.env|*/.env|*/.DS_Store|*.tdb)
      echo "Disallowed outgoing filename: ${file}" >&2
      bad_name=1
      ;;
  esac
  size="$(wc -c < "${file}" | tr -d ' ')"
  if [[ "${size}" -gt 5242880 ]]; then
    echo "Outgoing file exceeds 5 MiB: ${file} (${size} bytes)" >&2
    large_file=1
  fi
done < <(git diff --name-only --diff-filter=ACMR -z "${base_ref}")

if [[ "${bad_name}" -ne 0 || "${large_file}" -ne 0 ]]; then
  exit 1
fi

git diff --check "${base_ref}"

echo "Pre-publication checks passed against ${base_ref}."
