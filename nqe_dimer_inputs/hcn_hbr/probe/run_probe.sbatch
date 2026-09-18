#!/bin/bash
# Run hcnhbr_avtz_hessian_probe2.py as a batch job with the HCN...HBr tz-tier
# resources (16 cores, 64 GB, 47 h) so the probe can finish; an interactive
# allocation killed the first probe in stage E on its time limit.
#
# Usage:
#   sbatch run_probe.sbatch                      # probe2 next to this script
#   sbatch run_probe.sbatch /path/to/probe2.py   # or give the path
#
# All output goes to $OUTPUT_DIR as probe.<jobid>.out (stdout, incl. the
# heartbeat and tracebacks) and probe.<jobid>.err (stderr).  15 minutes before
# the walltime SLURM sends USR1 to this script, which forwards it to python;
# the probe then dumps every thread's traceback and its memory state.
#
#SBATCH --job-name=hcnhbr-probe
#SBATCH --account=brorsenk-lab
#SBATCH --partition=requeue
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=1-23:00:00
#SBATCH --signal=B:USR1@900
#SBATCH --output=/home/gbtmdr/hpc-group-storage/gabbie/dimers/pyscf-master/dimer-alpha-outputs/cneodft/hcnhbr/probe.%j.out
#SBATCH --error=/home/gbtmdr/hpc-group-storage/gabbie/dimers/pyscf-master/dimer-alpha-outputs/cneodft/hcnhbr/probe.%j.err

set -u

OUTPUT_DIR="${OUTPUT_DIR:-/home/gbtmdr/hpc-group-storage/gabbie/dimers/pyscf-master/dimer-alpha-outputs/cneodft/hcnhbr}"
PYSCF_DIR="${PYSCF_DIR:-/home/gbtmdr/hpc-group-storage/gabbie/dimers/pyscf-master}"
PYTHON="${PYTHON:-/mnt/pixstor/data/gbtmdr/miniconda/envs/dimers/bin/python3}"
CONDA_ENV="${CONDA_ENV:-dimers}"
PROBE="${1:-${SLURM_SUBMIT_DIR:-.}/hcnhbr_avtz_hessian_probe2.py}"

export PYSCF_MAX_MEMORY="${PYSCF_MAX_MEMORY:-32000}"   # half of --mem, as in the production tiers
export PROBE_HEARTBEAT_MIN="${PROBE_HEARTBEAT_MIN:-5}"

die() { echo "ERROR: $*" >&2; echo "ERROR: $*"; exit 1; }

mkdir -p "$OUTPUT_DIR" || die "cannot create $OUTPUT_DIR"
[ -f "$PROBE" ] || die "probe script not found: $PROBE"

echo "============================================================"
echo "HCN...HBr aug-cc-pVTZ Hessian probe (batch)"
echo "Job:       ${SLURM_JOB_ID:-<none>} on ${SLURM_JOB_NODELIST:-<none>}"
echo "Restarts:  ${SLURM_RESTART_COUNT:-0}"
echo "CPUs:      ${SLURM_CPUS_PER_TASK:-?}   Memory: ${SLURM_MEM_PER_NODE:-?} MB   PYSCF_MAX_MEMORY: $PYSCF_MAX_MEMORY MB"
echo "Probe:     $PROBE"
echo "Start:     $(date -Iseconds)"
echo "============================================================"

# --- Environment (same as submit_array.slurm) --------------------------------
for m in openblas gcc openmpi/4.1.1_gcc_9.5.0 miniconda3; do
    module load "$m" 2>/dev/null || echo "warning: module load $m failed (continuing)"
done
export PATH=$PATH:/home/gbtmdr/
export LD_LIBRARY_PATH=/cluster/software/SPACK/SPACK_v0.20_dev_a2/spack/opt/spack/linux-almalinux8-x86_64/gcc-9.5.0/openmpi-4.1.1-lzg4xxgoifgdslspsyr5w5iaisyz2r5r/lib:${LD_LIBRARY_PATH:-}
export LD_LIBRARY_PATH=/home/gbtmdr/openmpi/:$LD_LIBRARY_PATH
source activate "$CONDA_ENV" 2>/dev/null || echo "warning: 'source activate $CONDA_ENV' failed; relying on the absolute PYTHON path"
[ -x "$PYTHON" ] || die "PYTHON is not an executable file: $PYTHON"

unset PYTHONPATH
[ -f "$PYSCF_DIR/pyscf/__init__.py" ] || die "PYSCF_DIR does not contain pyscf/__init__.py: $PYSCF_DIR"
export PYTHONPATH="$PYSCF_DIR"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

# --- Scratch -------------------------------------------------------------------
SCRATCH=""
for candidate in "/local/scratch/gbtmdr/${SLURM_JOB_ID:-probe}" \
                 "${TMPDIR:-/tmp}/gbtmdr.${SLURM_JOB_ID:-probe}" \
                 "$OUTPUT_DIR/scratch.${SLURM_JOB_ID:-probe}"; do
    if mkdir -p "$candidate" 2>/dev/null && [ -w "$candidate" ]; then
        SCRATCH="$candidate"
        break
    fi
done
[ -n "$SCRATCH" ] || die "no writable scratch directory found"
export PYSCF_TMPDIR="$SCRATCH"
echo "Scratch:   $SCRATCH"
cd "$SCRATCH" || die "cannot cd to $SCRATCH"

# --- Run python in the background so signals can be forwarded to it ------------
"$PYTHON" -u "$PROBE" &
PID=$!
trap 'echo "wrapper: $(date -Iseconds) USR1 from SLURM (15 min to walltime), forwarding to python"; kill -USR1 $PID 2>/dev/null' USR1
trap 'echo "wrapper: $(date -Iseconds) TERM from SLURM, forwarding to python"; kill -TERM $PID 2>/dev/null' TERM
rc=0
while kill -0 "$PID" 2>/dev/null; do
    wait "$PID"
    rc=$?
done

echo "------------------------------------------------------------"
echo "python exit status: $rc   end: $(date -Iseconds)"
cp "$SCRATCH"/*.npz "$OUTPUT_DIR/" 2>/dev/null
cd / && rm -rf "$SCRATCH"
sacct -j "${SLURM_JOB_ID:-0}" --format=JobID,State,Elapsed,MaxRSS,MaxVMSize,ReqMem,ExitCode -P 2>/dev/null
exit $rc
