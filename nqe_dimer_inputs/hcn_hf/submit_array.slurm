#!/bin/bash -l
#
# SLURM submission for the HCN...HF CNEO-DFT re-optimization + Hessian
# inputs (hcn_hf/inputs/HCNHF.cneo-dft-*.py) on the lab partition. One
# array task per input, one node per task.
#
# Dual-purpose script (same pattern as the mc-development submit scripts):
#   * Run it directly (./submit_array.slurm) on the login node and it
#     submits itself as three arrays, one per electronic-basis tier.
#   * Under SLURM (sbatch ... submit_array.slurm <dz|tz|qz>) it is the
#     array job script.
#
# Tier resources, sized for a 2-day queue limit. A 5-atom CNEO-DFT job is
# cheap at DZ/TZ. At QZ the February 2026 runs (single-threaded) needed
# ~200 s per SCF and ~340 s per gradient, and their Hessians crashed in the
# coupled-perturbed solver after 100 Krylov iterations rather than timing
# out; the inputs now allow 300 iterations with a level-shift retry, and
# the QZ tier gets a full node's threads. PYSCF_MAX_MEMORY (MB) is exported
# per tier and read by the inputs.
#
#   tier  files                            mem    cpus  walltime  PYSCF_MAX_MEMORY
#   dz    *.aug-cc-pvdz.*.py                16G     8   1 day      12000
#   tz    *.aug-cc-pvtz.*.py                64G    16   2 days     48000
#   qz    *.aug-cc-pvqz.*.py               128G    32   2 days    100000
#
# Python runs under `timeout` ~2 h short of the walltime so partial results
# (trajectory, optimized geometry, log) still copy back to OUTPUT_DIR.
#
# Per-input output/error files (named like the input, run date appended):
#   <OUTPUT_DIR>/<input-stem>.<dd-mm-yyyy>.out / .err
# plus <stem>.traj.xyz, <stem>.opt.xyz, <stem>.hessian.npz.

#SBATCH --partition=brorsenk-lab
#SBATCH --account=brorsenk-lab
#SBATCH -J hcnhf-cneo-dft
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH -t 1-00:00:00
#SBATCH --output=/home/gbtmdr/hpc-group-storage/gabbie/nqe-dimer-output/slurm-%x.%A-%a.out
#SBATCH --error=/home/gbtmdr/hpc-group-storage/gabbie/nqe-dimer-output/slurm-%x.%A-%a.err

# --- Locations (override in the environment) ---------------------------
# INPUT_DIR defaults to the inputs/ folder next to this script; it is
# exported to the array tasks, so no copying is needed.
OUTPUT_DIR="${OUTPUT_DIR:-$HOME/hpc-group-storage/gabbie/nqe-dimer-output}"
# Python interpreter whose environment holds the theorychemyang/pyscf fork
# (>= 2.14) and geomeTRIC.
PYTHON="${PYTHON:-/mnt/pixstor/data/gbtmdr/miniconda/envs/densth/bin/python3}"

GLOB_PREFIX="HCNHF.cneo-dft-*"
shopt -s nullglob

# --- Driver mode: no tier argument -> submit the three tiers -----------
if [ $# -eq 0 ]; then
    if [ -n "$SLURM_JOB_ID" ] && [ -n "$SLURM_SUBMIT_DIR" ] \
            && [ -f "$SLURM_SUBMIT_DIR/$(basename "$0")" ]; then
        self="$SLURM_SUBMIT_DIR/$(basename "$0")"
    else
        self="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
    fi
    export INPUT_DIR="${INPUT_DIR:-$(dirname "$self")/inputs}"
    export OUTPUT_DIR PYTHON
    mkdir -p "$OUTPUT_DIR"
    for tier in dz tz qz; do
        files=("$INPUT_DIR"/$GLOB_PREFIX."aug-cc-pv$tier".*.py)
        n=${#files[@]}
        if [ "$n" -eq 0 ]; then
            echo "tier $tier: no $GLOB_PREFIX inputs in $INPUT_DIR, skipping"
            continue
        fi
        case "$tier" in
            dz) opts=(--mem=16G  --cpus-per-task=8  -t 1-00:00:00) ;;
            tz) opts=(--mem=64G  --cpus-per-task=16 -t 2-00:00:00) ;;
            qz) opts=(--mem=128G --cpus-per-task=32 -t 2-00:00:00) ;;
        esac
        echo "tier $tier: $n files"
        # Pass the locations explicitly so task mode does not depend on the
        # cluster's default --export policy.
        sbatch --job-name="hcnhf-cneo-dft-$tier" "${opts[@]}" \
               --export="ALL,INPUT_DIR=$INPUT_DIR,OUTPUT_DIR=$OUTPUT_DIR,PYTHON=$PYTHON" \
               --output="$OUTPUT_DIR/slurm-%x.%A-%a.out" \
               --error="$OUTPUT_DIR/slurm-%x.%A-%a.err" \
               --array="0-$((n - 1))" "$self" "$tier"
    done
    exit 0
fi

# --- Task mode: pick this array task's input file ----------------------
tier="${1:?usage: sbatch ... submit_array.slurm <dz|tz|qz>}"
INPUT_DIR="${INPUT_DIR:?INPUT_DIR must be exported by the driver}"
case "$tier" in
    dz) FILES=("$INPUT_DIR"/$GLOB_PREFIX.aug-cc-pvdz.*.py); RUN_TIMEOUT=22h; export PYSCF_MAX_MEMORY=12000 ;;
    tz) FILES=("$INPUT_DIR"/$GLOB_PREFIX.aug-cc-pvtz.*.py); RUN_TIMEOUT=46h; export PYSCF_MAX_MEMORY=48000 ;;
    qz) FILES=("$INPUT_DIR"/$GLOB_PREFIX.aug-cc-pvqz.*.py); RUN_TIMEOUT=46h; export PYSCF_MAX_MEMORY=100000 ;;
    *)  echo "unknown tier '$tier' (expected dz|tz|qz)" >&2; exit 2 ;;
esac

if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "task mode ($tier) needs SLURM_ARRAY_TASK_ID -- submit with NO argument:" \
         "./submit_array.slurm" >&2
    exit 2
fi
if [ "$SLURM_ARRAY_TASK_ID" -ge "${#FILES[@]}" ]; then
    echo "array index $SLURM_ARRAY_TASK_ID out of range (${#FILES[@]} $tier files)" >&2
    exit 2
fi

input="${FILES[$SLURM_ARRAY_TASK_ID]}"
filename=$(basename "$input")
stem="${filename%.py}"
current_date=$(date +%d-%m-%Y)

mkdir -p "$OUTPUT_DIR"
exec 1>"$OUTPUT_DIR/$stem.$current_date.out" 2>"$OUTPUT_DIR/$stem.$current_date.err"

echo "Host:   $(hostname)"
echo "Job:    $SLURM_JOB_ID (array $SLURM_ARRAY_JOB_ID/$SLURM_ARRAY_TASK_ID, tier $tier)"
echo "Input:  $input"
echo "Cores:  $SLURM_CPUS_PER_TASK   Mem: $SLURM_MEM_PER_NODE MB   PYSCF_MAX_MEMORY: $PYSCF_MAX_MEMORY MB"
echo "Start:  $(date -Iseconds)"
echo "------------------------------------------------------------"

# --- Environment -------------------------------------------------------
module load openblas
module load gcc
module load openmpi/4.1.1_gcc_9.5.0

export PATH=$PATH:/home/gbtmdr/
export LD_LIBRARY_PATH=/cluster/software/SPACK/SPACK_v0.20_dev_a2/spack/opt/spack/linux-almalinux8-x86_64/gcc-9.5.0/openmpi-4.1.1-lzg4xxgoifgdslspsyr5w5iaisyz2r5r/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/gbtmdr/openmpi/:$LD_LIBRARY_PATH

module load miniconda3
source activate densth

# Scrub any login-shell PYTHONPATH so $PYTHON resolves packages from its
# own environment only.
unset PYTHONPATH

echo "Python: $PYTHON"
"$PYTHON" - <<'EOF' || { echo "environment sanity check failed" >&2; exit 1; }
import inspect, sys, numpy, pyscf, pyscf.neo, geometric
from pyscf.hessian import thermo
print("  ", sys.version.split()[0], "numpy", numpy.__version__,
      "pyscf", pyscf.__version__, "geometric", geometric.__version__)
print("   pyscf from", pyscf.__file__)
if "1e-9" not in inspect.getsource(thermo.rotation_const):
    print("   WARNING: this pyscf build lacks the linear-rotor fix in thermo.rotation_const"
          " (fork >= 2.14); a linear molecule may be given 3N-6 modes")
EOF

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

# --- Node-local scratch (one dir per array task) ------------------------
SCRATCH="/local/scratch/gbtmdr/$SLURM_JOB_ID"
mkdir -p "$SCRATCH"
export PYSCF_TMPDIR="$SCRATCH"
cp "$input" "$SCRATCH/"
cd "$SCRATCH" || exit 1

# --- Run (kill python ~2 h before the walltime so results copy back) ----
timeout "$RUN_TIMEOUT" "$PYTHON" -u "$filename"
status=$?

echo "------------------------------------------------------------"
echo "End:        $(date -Iseconds)"
echo "Exit code:  $status"
ls -l

cp "$SCRATCH"/*.xyz "$SCRATCH"/*.npz "$OUTPUT_DIR/" 2>/dev/null
rm -rf "$SCRATCH"

exit $status
