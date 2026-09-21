"""Full-precision reproducibility probe: per-rank sha256 of the float64
solver state after N steps. frt.txt is E18.7E4 (7 digits), so a byte-equal
frt does NOT prove the reduction was bit-stable -- this hashes the doubles."""
import hashlib, sys
sys.path.insert(0, '/home/utig5/dliu/wt-jaxshard/src/python')
from mpi4py import MPI
from eqdyna import eqdyna3d, driver, backend as B
case, nsteps = sys.argv[1], int(sys.argv[2])
comm = MPI.COMM_WORLD
S, mesh = eqdyna3d.build_solver_state(case)
out = driver.run_mpi(S, comm, nsteps=nsteps, verbose=False,
                     xp=B.array_module('jax'))
h = {k: hashlib.sha256(out[k].tobytes()).hexdigest()[:16]
     for k in ('velArr', 'dispArr', 'force', 'fric', 'fnft')}
print('HASH rank %d %s' % (comm.Get_rank(), ' '.join('%s=%s' % kv for kv in sorted(h.items()))), flush=True)
