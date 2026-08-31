#!/usr/bin/env bash
# Fail if class_src/ is out of date with respect to R/.
# Run after editing anything in R/:   bash scripts/check_class_src.sh
set -e
cd "$(dirname "$0")/.."
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
cp -r class_src "$tmp/before"
python3 scripts/make_class_src.py > /dev/null
if diff -rq "$tmp/before" class_src > /dev/null; then
  echo "class_src/ is in sync with R/"
else
  echo "class_src/ was STALE and has been regenerated:"
  diff -rq "$tmp/before" class_src | sed 's/^/  /'
  exit 1
fi
