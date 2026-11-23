#!/bin/bash
# Run one VoleMUPs step inside WSL (or natively on the Mac/Linux) with the micromamba env "genomics";
# HMMER and seqkit come from the env "vmsearch". Usage (from Windows PowerShell):
#   wsl -d Ubuntu -- bash /mnt/c/GitHub/VoleMUPs/scripts/run_wsl.sh step02_annotated.py [args]
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENVS="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs"
export PATH="$ENVS/genomics/bin:$ENVS/vmsearch/bin:$PATH"
cd "$HERE"
exec nice -n 10 python -u "$@"
