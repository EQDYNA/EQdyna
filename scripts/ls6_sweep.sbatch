#!/bin/bash
# The full test sweep (`testsys/run.py all`) as one Lonestar6 job. From the
# repo root, after `./install-eqdyna.sh -e ls6`:
#   sbatch -A <allocation> scripts/ls6_sweep.sbatch
# Output: ls6_sweep_<jobid>.tgz (sweep log + the perf/profile rows it appended).
#SBATCH -J eqdyna-sweep
#SBATCH -o ls6_sweep_%j.log
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -t 02:00:00

cd "$SLURM_SUBMIT_DIR"
source ./install-eqdyna.sh -c ls6
export EQDYNA_TEST_MACHINE=ls6
git rev-parse HEAD; module list 2>&1; python3 -c "import jax, mpi4py; print('jax', jax.__version__, 'mpi4py', mpi4py.__version__)"

L0=$(wc -l < docs/perf_ledger.jsonl); P0=$(wc -l < docs/run_profiles.jsonl)
python3 testsys/run.py all
echo "run.py all exit $?"
tail -n +$((L0 + 1)) docs/perf_ledger.jsonl > ls6_sweep_$SLURM_JOB_ID.ledger.jsonl
tail -n +$((P0 + 1)) docs/run_profiles.jsonl > ls6_sweep_$SLURM_JOB_ID.profiles.jsonl
tar czf ls6_sweep_$SLURM_JOB_ID.tgz ls6_sweep_$SLURM_JOB_ID.*
