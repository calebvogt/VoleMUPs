#!/bin/bash
# Run the Linux steps of the VoleMUPs pipeline in order (inside WSL or on the Mac/Linux), one at a time so that the
# thread budget in config.yml holds. Steps 01 (downloads) and 06 --download need internet and can run on Windows.
#   bash run_all.sh                       # everything after the downloads
#   bash run_all.sh step04 step05         # only some steps
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R="$HERE/run_wsl.sh"
steps=("$@")
[ ${#steps[@]} -eq 0 ] && steps=(step02 step03 step04 step06 step05 step07)
for s in "${steps[@]}"; do
  echo "=== $s $(date '+%F %T')"
  case $s in
    step02) bash "$R" step02_annotated.py ;;
    step03) bash "$R" step03_regions.py --define --fetch --search ;;
    step03s) bash "$R" step03_regions.py --search ;;
    step04) bash "$R" step04_genomes.py ;;
    step04_*) bash "$R" step04_genomes.py --only "${s#step04_}" ;;
    step05) bash "$R" step05_tree.py ;;
    step05f) bash "$R" step05_tree.py --skip-focused ;;
    step06) bash "$R" step06_expression.py --count ;;
    step07) bash "$R" step07_collect.py ;;
    *) echo "unknown step $s" ;;
  esac
  rc=$?
  echo "=== $s done $(date '+%F %T') (exit $rc)"
done
