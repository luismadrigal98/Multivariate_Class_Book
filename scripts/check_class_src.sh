#!/usr/bin/env bash
# Report whether class_src/ is up to date with R/.
#
#   bash scripts/check_class_src.sh
#
# Files carrying a HAND-OWNED line in their header are maintained by hand; the
# generator skips them and so does this check. Fixes made in R/ do NOT reach
# them -- they are listed separately so they never go quietly stale.
set -e
cd "$(dirname "$0")/.."
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
cp -r class_src "$tmp/before"
python3 scripts/make_class_src.py > /dev/null

owned=$(grep -l 'HAND-OWNED' class_src/*.R 2>/dev/null | xargs -r -n1 basename || true)

if diff -rq "$tmp/before" class_src > /dev/null; then
  echo "class_src/ is in sync with R/"
else
  echo "class_src/ was STALE and has been regenerated:"
  diff -rq "$tmp/before" class_src | sed 's/^/  /'
fi

if [ -n "$owned" ]; then
  echo
  echo "hand-owned (generator skipped, R/ fixes do not reach these):"
  for f in $owned; do echo "  $f"; done
fi
