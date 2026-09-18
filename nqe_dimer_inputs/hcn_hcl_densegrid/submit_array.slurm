#!/bin/bash -l
#
# SLURM submission for the HCN...HCl CNEO-DFT (dense chlorine grid variant, cl400) re-optimization + Hessian
# inputs (HCNHCl.cneo-dft-*.py) on the requeue partition, brorsenk-lab
# account. One array task per input, one node per task.
#
# Dual-purpose script (same pattern as the mc-development submit scripts):
#   * Run it directly (./submit_array.slurm) on the login node and it
#     submits itself as three arrays, one per electronic-basis tier.
#   * Under SLURM (sbatch ... submit_array.slurm <dz|tz|qz>) it is the
#     array job script.
#
# Locations (override any of them in the environment before running):
#   INPUT_DIR   /home/gbtmdr/hpc-group-storage/gabbie/dimers/pyscf-master/dimer-alpha-inputs/cneodft/hcnhcl-densegrid
#   OUTPUT_DIR  /home/gbtmdr/hpc-group-storage/gabbie/dimers/pyscf-master/dimer-alpha-outputs/cneodft/hcnhcl-densegrid
#               (SLURM .out/.err, per-input .out/.err logs, .traj.xyz, .opt.xyz, .hessian.npz)
#   PYSCF_DIR   /home/gbtmdr/hpc-group-storage/gabbie/dimers/pyscf-master
#               (the directory that CONTAINS pyscf/; put on PYTHONPATH so the
#               source-tree build of the theorychemyang fork is imported)
#   PYTHON      the dimers conda environment's interpreter
#   CONDA_ENV   dimers
#
# Every tier finishes before 48 hours: walltime 1-23:00:00 (47 h) at most,
# and python runs under `timeout` two hours short of the walltime so a run
# that would hit the limit is stopped early and its partial results
# (trajectory, optimized geometry, log) still copy back to OUTPUT_DIR.
#
#   tier  files                            mem    cpus  walltime   timeout  PYSCF_MAX_MEMORY
#   dz    *.aug-cc-pvdz.*.py                16G     8   1-00:00:00   22h       8000
#   tz    *.aug-cc-pvtz.*.py                64G    16   1-23:00:00   45h      32000
#   qz    *.aug-cc-pvqz.*.py               128G    32   1-23:00:00   45h      64000
#
# The per-input .out log announces every stage ("stage: ...") and always
# ends with "script exit status N", so a run that stops early shows where.
# Fatal errors are written to both the .out and the .err file.
#
# The requeue partition can preempt jobs; --requeue puts a preempted task
# back in the queue and it restarts from the input (scratch is lost).
# Restarted attempts log to <stem>.<date>.r<N>.out/.err so the first
# attempt's log is kept.

#SBATCH --partition=requeue
#SBATCH --account=brorsenk-lab
#SBATCH --requeue
#SBATCH -J hcnhcl-densegrid-cneo-dft
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH -t 1-00:00:00
#SBATCH --output=/home/gbtmdr/hpc-group-storage/gabbie/dimers/pyscf-master/dimer-alpha-outputs/cneodft/hcnhcl-densegrid/slurm-%x.%A-%a.out
#SBATCH --error=/home/gbtmdr/hpc-group-storage/gabbie/dimers/pyscf-master/dimer-alpha-outputs/cneodft/hcnhcl-densegrid/slurm-%x.%A-%a.err

# --- Locations ---------------------------------------------------------
INPUT_DIR="${INPUT_DIR:-/home/gbtmdr/hpc-group-storage/gabbie/dimers/pyscf-master/dimer-alpha-inputs/cneodft/hcnhcl-densegrid}"
OUTPUT_DIR="${OUTPUT_DIR:-/home/gbtmdr/hpc-group-storage/gabbie/dimers/pyscf-master/dimer-alpha-outputs/cneodft/hcnhcl-densegrid}"
PYSCF_DIR="${PYSCF_DIR:-/home/gbtmdr/hpc-group-storage/gabbie/dimers/pyscf-master}"
PYTHON="${PYTHON:-/mnt/pixstor/data/gbtmdr/miniconda/envs/dimers/bin/python3}"
CONDA_ENV="${CONDA_ENV:-dimers}"

GLOB_PREFIX="HCNHCl.cneo-dft-*"
shopt -s nullglob

# --- Driver mode: no tier argument -> submit the three tiers -----------
if [ $# -eq 0 ]; then
    if [ -n "$SLURM_JOB_ID" ] && [ -n "$SLURM_SUBMIT_DIR" ] \
            && [ -f "$SLURM_SUBMIT_DIR/$(basename "$0")" ]; then
        self="$SLURM_SUBMIT_DIR/$(basename "$0")"
    else
        self="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
    fi
    if [ ! -d "$INPUT_DIR" ]; then
        echo "INPUT_DIR does not exist: $INPUT_DIR" >&2
        exit 1
    fi
    mkdir -p "$OUTPUT_DIR" || { echo "cannot create OUTPUT_DIR: $OUTPUT_DIR" >&2; exit 1; }
    for tier in dz tz qz; do
        files=("$INPUT_DIR"/$GLOB_PREFIX."aug-cc-pv$tier".*.py)
        n=${#files[@]}
        if [ "$n" -eq 0 ]; then
            echo "tier $tier: no $GLOB_PREFIX inputs in $INPUT_DIR, skipping"
            continue
        fi
        case "$tier" in
            dz) opts=(--mem=16G  --cpus-per-task=8  -t 1-00:00:00) ;;
            tz) opts=(--mem=64G  --cpus-per-task=16 -t 1-23:00:00) ;;
            qz) opts=(--mem=128G --cpus-per-task=32 -t 1-23:00:00) ;;
        esac
        echo "tier $tier: $n files"
        # Pass the locations explicitly so task mode does not depend on the
        # cluster's default --export policy.
        sbatch --job-name="hcnhcl-densegrid-cneo-dft-$tier" "${opts[@]}" \
               --export="ALL,INPUT_DIR=$INPUT_DIR,OUTPUT_DIR=$OUTPUT_DIR,PYSCF_DIR=$PYSCF_DIR,PYTHON=$PYTHON,CONDA_ENV=$CONDA_ENV" \
               --output="$OUTPUT_DIR/slurm-%x.%A-%a.out" \
               --error="$OUTPUT_DIR/slurm-%x.%A-%a.err" \
               --array="0-$((n - 1))" "$self" "$tier"
    done
    exit 0
fi

# --- Task mode: pick this array task's input file ----------------------
tier="${1:?usage: sbatch ... submit_array.slurm <dz|tz|qz>}"
case "$tier" in
    dz) FILES=("$INPUT_DIR"/$GLOB_PREFIX.aug-cc-pvdz.*.py); RUN_TIMEOUT=22h; export PYSCF_MAX_MEMORY=8000 ;;
    tz) FILES=("$INPUT_DIR"/$GLOB_PREFIX.aug-cc-pvtz.*.py); RUN_TIMEOUT=45h; export PYSCF_MAX_MEMORY=32000 ;;
    qz) FILES=("$INPUT_DIR"/$GLOB_PREFIX.aug-cc-pvqz.*.py); RUN_TIMEOUT=45h; export PYSCF_MAX_MEMORY=64000 ;;
    *)  echo "unknown tier '$tier' (expected dz|tz|qz)" >&2; exit 2 ;;
esac

if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "task mode ($tier) needs SLURM_ARRAY_TASK_ID -- submit with NO argument:" \
         "./submit_array.slurm" >&2
    exit 2
fi
if [ "$SLURM_ARRAY_TASK_ID" -ge "${#FILES[@]}" ]; then
    echo "array index $SLURM_ARRAY_TASK_ID out of range (${#FILES[@]} $tier files in $INPUT_DIR)" >&2
    exit 2
fi

input="${FILES[$SLURM_ARRAY_TASK_ID]}"
filename=$(basename "$input")
stem="${filename%.py}"
current_date=$(date +%d-%m-%Y)
attempt="${SLURM_RESTART_COUNT:-0}"
logstem="$stem.$current_date"
if [ "$attempt" -gt 0 ]; then
    logstem="$logstem.r$attempt"
fi

# From here on, stdout/stderr go to files named like the input file.
mkdir -p "$OUTPUT_DIR"
exec 1>"$OUTPUT_DIR/$logstem.out" 2>"$OUTPUT_DIR/$logstem.err"

# Diagnostics: stage markers in .out, fatal errors in both files, and a
# final line in .out whatever happens.
stage() { echo "$(date -Iseconds) stage: $*"; }
die()   { echo "$(date -Iseconds) FATAL: $*"; echo "FATAL: $*" >&2; exit 1; }
trap 'rc=$?; echo "$(date -Iseconds) script exit status $rc"' EXIT

echo "Host:     $(hostname)"
echo "Job:      $SLURM_JOB_ID (array $SLURM_ARRAY_JOB_ID/$SLURM_ARRAY_TASK_ID, tier $tier, attempt $attempt)"
echo "Input:    $input"
echo "Cores:    $SLURM_CPUS_PER_TASK   Mem: $SLURM_MEM_PER_NODE MB   PYSCF_MAX_MEMORY: $PYSCF_MAX_MEMORY MB"
echo "Timeout:  $RUN_TIMEOUT"
echo "Start:    $(date -Iseconds)"
echo "------------------------------------------------------------"

# --- Environment -------------------------------------------------------
stage "loading modules"
for m in openblas gcc openmpi/4.1.1_gcc_9.5.0 miniconda3; do
    if ! module load "$m" 2>>"$OUTPUT_DIR/$logstem.err"; then
        echo "warning: module load $m failed (continuing)"
    fi
done

export PATH=$PATH:/home/gbtmdr/
export LD_LIBRARY_PATH=/cluster/software/SPACK/SPACK_v0.20_dev_a2/spack/opt/spack/linux-almalinux8-x86_64/gcc-9.5.0/openmpi-4.1.1-lzg4xxgoifgdslspsyr5w5iaisyz2r5r/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/gbtmdr/openmpi/:$LD_LIBRARY_PATH

stage "activating conda environment $CONDA_ENV"
if ! source activate "$CONDA_ENV" 2>>"$OUTPUT_DIR/$logstem.err"; then
    echo "warning: 'source activate $CONDA_ENV' failed; relying on the absolute PYTHON path"
fi
[ -x "$PYTHON" ] || die "PYTHON is not an executable file: $PYTHON"

# Scrub any login-shell PYTHONPATH (a stale one mixed a 3.12 numpy into a
# 3.10 interpreter before), then point at the pyscf source tree.
unset PYTHONPATH
if [ -n "$PYSCF_DIR" ]; then
    [ -f "$PYSCF_DIR/pyscf/__init__.py" ] \
        || die "PYSCF_DIR does not contain pyscf/__init__.py; point it at the PARENT of pyscf/: $PYSCF_DIR"
    export PYTHONPATH="$PYSCF_DIR"
fi

stage "environment sanity check"
echo "Python:    $PYTHON"
echo "PYSCF_DIR: ${PYSCF_DIR:-<unset: using the pyscf installed in the environment>}"
"$PYTHON" - <<'EOF' || die "environment sanity check failed (see the .err file for the traceback)"
import inspect, sys, numpy, scipy, pyscf, pyscf.neo, geometric
from pyscf.hessian import thermo
print("  ", sys.executable)
print("  ", sys.version.split()[0], "numpy", numpy.__version__, "scipy", scipy.__version__,
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

# --- Scratch: node-local if available, otherwise a fallback -------------
stage "setting up scratch"
SCRATCH=""
for candidate in "/local/scratch/gbtmdr/$SLURM_JOB_ID" \
                 "${TMPDIR:-/tmp}/gbtmdr.$SLURM_JOB_ID" \
                 "$OUTPUT_DIR/scratch.$SLURM_JOB_ID"; do
    if mkdir -p "$candidate" 2>/dev/null && [ -w "$candidate" ]; then
        SCRATCH="$candidate"
        break
    fi
    echo "warning: cannot use scratch directory $candidate"
done
[ -n "$SCRATCH" ] || die "no writable scratch directory found"
echo "Scratch:   $SCRATCH"
export PYSCF_TMPDIR="$SCRATCH"
cp "$input" "$SCRATCH/" || die "cannot copy $input to $SCRATCH"
cd "$SCRATCH" || die "cannot cd to $SCRATCH"

# --- Run (kill python 2 h before the walltime so results copy back) -----
stage "running: timeout $RUN_TIMEOUT $PYTHON -u $filename"
timeout "$RUN_TIMEOUT" "$PYTHON" -u "$filename"
status=$?
if [ "$status" -eq 124 ]; then
    echo "python stopped by timeout ($RUN_TIMEOUT) before the walltime"
    echo "python stopped by timeout ($RUN_TIMEOUT) before the walltime" >&2
fi

echo "------------------------------------------------------------"
echo "End:        $(date -Iseconds)"
echo "Exit code:  $status"
stage "copying results to $OUTPUT_DIR"
ls -l
cp "$SCRATCH"/*.xyz "$SCRATCH"/*.npz "$OUTPUT_DIR/" 2>/dev/null
cd / && rm -rf "$SCRATCH"

exit $status
