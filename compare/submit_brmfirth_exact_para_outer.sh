#!/bin/bash
# Outer-parallel simulation using run_brmfirth_exact_para.R.
# The exact-CI bootstrap inside each simulation replicate is forced to serial.
#
# Usage:
#   sbatch submit_brmfirth_exact_para_outer.sh \
#     [n] [param] [event] [hypothesis] [result_dir] [outer_ncores]
#
# Example:
#   sbatch submit_brmfirth_exact_para_outer.sh \
#     100 RR rare alternative /path/to/results 100

#SBATCH --job-name=brm_exact_outer
#SBATCH --account=def-liteep
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=104
#SBATCH --mem=60G
#SBATCH --output=brm_exact_outer_%j.out
#SBATCH --error=brm_exact_outer_%j.err

set -euo pipefail

module load StdEnv/2023
module load gcc/12.3 r/4.3.1

# Prevent nested BLAS/OpenMP parallelism in every outer worker.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
code_dir="${SLURM_SUBMIT_DIR:-$script_dir}"
if [[ ! -f "$code_dir/run_brmfirth_exact_para.R" ]]; then
  code_dir="$script_dir"
fi
run_file="${code_dir}/run_brmfirth_exact_para.R"
ci_file="${code_dir}/CI_exact_adapt_para.R"

n="${1:-100}"
param="${2:-RR}"
event="${3:-rare}"
hypothesis="${4:-alternative}"
result_dir="${5:-$code_dir/results/exact_para}"
requested_outer_ncores="${6:-100}"

for integer_arg in "${n}" "${requested_outer_ncores}"; do
  if ! [[ "${integer_arg}" =~ ^[1-9][0-9]*$ ]]; then
    echo "n and outer_ncores must be positive integers" >&2
    exit 2
  fi
done

case "${param}" in
  RR|RD) ;;
  *) echo "param must be RR or RD" >&2; exit 2 ;;
esac
case "${event}" in
  common|rare) ;;
  *) echo "event must be common or rare" >&2; exit 2 ;;
esac
case "${hypothesis}" in
  null|alternative) ;;
  *) echo "hypothesis must be null or alternative" >&2; exit 2 ;;
esac

if [[ ! -f "${run_file}" ]]; then
  echo "Cannot find ${run_file}" >&2
  exit 2
fi
if [[ ! -f "${ci_file}" ]]; then
  echo "Cannot find ${ci_file}" >&2
  exit 2
fi

# Leave CPUs available for the master R process and OS. The default allocation
# is 104 CPUs, so this normally starts 100 single-threaded outer workers.
allocated_cpus="${SLURM_CPUS_PER_TASK:-104}"
worker_limit=$((allocated_cpus - 1))
if (( worker_limit < 1 )); then
  worker_limit=1
fi

outer_ncores="${requested_outer_ncores}"
if (( outer_ncores > worker_limit )); then
  outer_ncores="${worker_limit}"
fi
if (( outer_ncores > 1000 )); then
  outer_ncores=1000
fi

# These two values deliberately disable the inner parallel layer.
exact_parallel=FALSE
exact_ncores=1

cd "${code_dir}"
mkdir -p "${result_dir}"

echo "Starting BRM-Firth exact simulation"
echo "param=${param}, event=${event}, hypothesis=${hypothesis}"
echo "n=${n}, seeds=1:1000"
echo "outer_ncores=${outer_ncores}"
echo "exact_parallel=${exact_parallel}, exact_ncores=${exact_ncores}"
echo "result_dir=${result_dir}"

Rscript --vanilla --max-connections=512 "${run_file}" \
  "param='${param}'" \
  "event='${event}'" \
  "hypothesis='${hypothesis}'" \
  "n=${n}" \
  "ncores=${outer_ncores}" \
  "exact_parallel=${exact_parallel}" \
  "exact_ncores=${exact_ncores}" \
  "result_dir='${result_dir}'"

echo "BRM-Firth exact simulation completed"
