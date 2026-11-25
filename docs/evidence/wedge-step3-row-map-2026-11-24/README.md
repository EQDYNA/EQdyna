# Row 167 step 3: do the changed tpv36/37 reference rows equal the #281 rewiring? (2026-11-24)

Evidence only. Question: PR #284 regenerates the test.tpv36/37 references after PR #281 (one geometric side test in
`replaceSlaveWithMasterNode`). Is the change in the reference exactly the rewired elements, and nothing else?

Method (no solver cell run; python backend, mesh only, one core): `dump_mesh.py <root> <case> <out.npz>` builds the mesh at
the tree before #281 (`ecb8ccb`) and after (`a22b9dd`) and saves meshCoor, conn, elem_type, nsmp, un;
`map_rows.py pre.npz post.npz frt_old frt_new` lists the elements whose connectivity changed, the frt rows (matched on
rounded x,y,z) of every node of those elements (E), and the frt rows whose stress columns 11-13 changed by >1e-3 Pa (R).
frt_old = `origin/master:test.reference.results/<case>/frt.canonical.txt` (a22b9dd), frt_new = PR #284 head `3628339`.

Result (`map_rows_output.txt`): both cases, 92 elements rewired; meshCoor, element types and slave-master pairs bit-equal
(only conn changes); E = 105 fault rows of 3477, R = 100, R minus E = 0, E minus R = 5; coordinates and rupture time
identical; max stress change 231900 Pa. The tpv36 and tpv37 meshes are md5-identical (same geometry) while their reference
files differ, so the matching numbers are not a wrong-file read. A strict any-bit diff is 3196 / 3174 rows (1e-17 roundoff in
cols 4-10). The producer's "156 rows" is not reproduced on frt (100); not explained.

Scope: this maps the frt reference only. fault.dyna.r.nc and the station files were not mapped, and nothing here validates
the new references against an independent oracle (none exists for TPV36/37: no FaultMod dataset in the shared dataset index).
