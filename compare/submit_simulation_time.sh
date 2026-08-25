#!/bin/bash
# Submit with, for example:
#   sbatch submit_simulation_time.sh 100 RR rare alternative
# Use 32 cores for the inner exact() bootstrap:
#   sbatch --cpus-per-task=32 submit_simulation_time.sh 100 RR rare alternative /path/to/results 32
# Arguments:
#   1 n, 2 param, 3 event, 4 hypothesis, 5 result_dir, 6 exact_ncores

#SBATCH --job-name=method_time
#SBATCH --account=def-liteep
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=104
#SBATCH --mem=60G
#SBATCH --output=method_time_%j.out
#SBATCH --error=method_time_%j.err

set -euo pipefail

module load StdEnv/2023
module load gcc/12.3 r/4.3.1

# Each inner exact() worker must remain single-threaded.  Otherwise every
# worker can start additional BLAS/OpenMP threads and oversubscribe the node.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
code_dir="${SLURM_SUBMIT_DIR:-$script_dir}"
if [[ ! -f "$code_dir/run_simulation_time.R" ]]; then
  code_dir="$script_dir"
fi
n="${1:-100}"
param="${2:-RR}"
event="${3:-rare}"
hypothesis="${4:-alternative}"
result_dir="${5:-$code_dir/results/time}"
exact_ncores="${6:-100}"

if ! [[ "${exact_ncores}" =~ ^[1-9][0-9]*$ ]]; then
  echo "exact_ncores must be a positive integer" >&2
  exit 2
fi

# Always use outer-serial/inner-parallel mode in run_simulation_time.R.
ncores=1
exact_parallel=TRUE
allocated_cpus="${SLURM_CPUS_PER_TASK:-104}"
if (( exact_ncores > allocated_cpus )); then
  echo "exact_ncores (${exact_ncores}) exceeds allocated CPUs (${allocated_cpus})" >&2
  exit 2
fi

cd "${code_dir}"
mkdir -p "${result_dir}"

echo "Starting timing simulation"
echo "param=${param}, event=${event}, hypothesis=${hypothesis}"
echo "n=${n}, seeds=1:1000"
echo "outer_ncores=${ncores}, exact_parallel=${exact_parallel}, exact_ncores=${exact_ncores}"
echo "result_dir=${result_dir}"

Rscript --vanilla "${code_dir}/run_simulation_time.R" \
  "param='${param}'" \
  "event='${event}'" \
  "hypothesis='${hypothesis}'" \
  "n=${n}" \
  "ncores=${ncores}" \
  "exact_parallel=${exact_parallel}" \
  "exact_ncores=${exact_ncores}" \
  "result_dir='${result_dir}'"

echo "Timing simulation completed"
