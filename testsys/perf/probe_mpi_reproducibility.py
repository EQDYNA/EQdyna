"""Full-precision reproducibility probe: per-rank sha256 of the float64
solver state after N steps. frt.txt is E18.7E4 (7 digits), so a byte-equal
frt does NOT prove the reduction was bit-stable -- this hashes the doubles."""
import hashlib, os, sys
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.environ.get('EQDYNAROOT') or os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0, os.path.join(_ROOT, 'src', 'python'))
from mpi4py import MPI
from eqdyna import eqdyna3d, driver, backend as B
case, nsteps = sys.argv[1], int(sys.argv[2])
comm = MPI.COMM_WORLD
S, mesh = eqdyna3d.build_solver_state(case)
out = driver.run_mpi(S, comm, nsteps=nsteps, verbose=False,
                     xp=B.array_module('jax'))
h = {k: hashlib.sha256(out[k].tobytes()).hexdigest()[:16]
     for k in ('velArr_local', 'dispArr_local', 'force_local',
               'fric', 'fnft')}
print('HASH rank %d %s' % (comm.Get_rank(), ' '.join('%s=%s' % kv for kv in sorted(h.items()))), flush=True)
