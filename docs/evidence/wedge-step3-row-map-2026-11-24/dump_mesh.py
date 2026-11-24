import sys, os, tempfile, numpy as np
root, case, out = sys.argv[1:4]
sys.path.insert(0, root); sys.path.insert(0, os.path.join(root,'testsys','regression'))
import test_mesh_screening as T
d = tempfile.mkdtemp(dir=os.path.dirname(out))
T.make_case(case, os.path.join(d,'c'))
real = T.meshscreen.screen_mesh; got={}
def grab(*a):
    got['a']=a
    try: return real(*a)
    except Exception as e: print('screen raised', type(e).__name__); return None
T.meshscreen.screen_mesh = grab
T.eqdyna3d.build_solver_state(os.path.join(d,'c'), mesh_only=True)
a=got['a']; np.savez(out, meshCoor=a[0], conn=a[1], elem_type=a[2], nsmp=a[3], un=a[4])
print(out, [np.shape(x) for x in a[:5]])
