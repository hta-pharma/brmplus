#!/bin/sh

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
code_dir="${SLURM_SUBMIT_DIR:-$script_dir}"
if [ ! -f "$code_dir/run_simulation.R" ]; then
  code_dir="$script_dir"
fi
result_dir="${RESULT_DIR:-$code_dir/results/simulation}"

module load StdEnv/2023
module load gcc/12.3 r/4.3.1

n_values=${1:-50,200,500}
param_values=${2:-RR,RD}
event_values=${3:-common,rare}
hypothesis_values=${4:-null,alternative}

cd "$code_dir"
mkdir -p "$result_dir/Rout"

R --vanilla --max-connections=512 CMD BATCH --no-save --no-restore \
  "--args n_values='$n_values' param_values='$param_values' event_values='$event_values' hypothesis_values='$hypothesis_values' result_dir='$result_dir'" \
  "$code_dir/run_simulation.R" \
  "$result_dir/Rout/simucpp_all_scenarios_R_1000.Rout"
