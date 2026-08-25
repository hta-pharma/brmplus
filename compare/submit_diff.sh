#!/bin/sh
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
code_dir="${SLURM_SUBMIT_DIR:-$script_dir}"
if [ ! -f "$code_dir/run_rarediff.R" ]; then
  code_dir="$script_dir"
fi
result_dir="${RESULT_DIR:-$code_dir/results/rarediff}"

module load StdEnv/2023
module load gcc/12.3   r/4.3.1


n=${1:-50}
param=${2:-RR}
event=${3:-rare11}
hypothesis=${4:-alternative}

case "$param" in
  RR|RD) ;;
  *) echo "param must be RR or RD" >&2; exit 2 ;;
esac
case "$event" in
  rare1|rare11) ;;
  *) echo "event must be rare1 or rare11" >&2; exit 2 ;;
esac
case "$hypothesis" in
  null|alternative) ;;
  *) echo "hypothesis must be null or alternative" >&2; exit 2 ;;
esac

cd "$code_dir"
mkdir -p "$result_dir/Rout"
R --vanilla --max-connections=512 CMD BATCH --no-save --no-restore \
  "--args n=$n param='$param' event='$event' hypothesis='$hypothesis' result_dir='${result_dir}'" \
  "$code_dir/run_rarediff.R" \
  "$result_dir/Rout/rarediff_${param}_${event}_${hypothesis}_N_${n}_R_1000.Rout"
