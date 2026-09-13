# testsys/parity/ — Python-port parity tier (pathway_forward item 14)

## Default-build neutrality (the release-critical property)

`src/pydump.f90` generates the fixtures this tier compares against, but it
must not change the production `eqdyna` binary. Design: `driver.f90`/
`eqdyna3d.f90` call `pydump_step`/`pydump_state` **unconditionally** (no
`#ifdef`, no preprocessor). The seam is at LINK time: `src/makefile`'s
default `eqdyna` target links `pydump_noop.o` (two empty subroutines,
call+return); `make PYDUMP=1 eqdyna-pydump` links the real `pydump.o`.

This means `eqdyna3d.o`/`driver.o`'s **compiled bytes are identical**
between the default and `PYDUMP=1` builds (same source text, same compile
flags, verified by construction — no macro or flag differs between the
two compiles of these two files). The only difference between the two
final binaries is which `pydump_*.o` gets linked in.

**An `#ifdef PYDUMP` + `-cpp` design was tried first and rejected.**
Empirically, turning on `-cpp` — even with the macro never defined, so the
guarded block is always stripped — changed the compiled **object file
size** for `eqdyna3d.o` (not just embedded metadata like `.comment` or
`.note.gnu.build-id`, which were also checked and ruled out as the cause).
The root cause was not isolated in the time available (candidate: cpp's
`# <linenum> "<file>"` line markers interacting with Fortran free-form
`&` continuation syntax in a way that shifts something in codegen even at
`-O3`). Rather than ship a mechanism with an unexplained object-code
effect, the link-time no-op-stub seam above was used instead — it needs no
preprocessor at all.

**Residual, non-eliminable difference vs the TRUE pre-pydump pristine
source** (i.e., before `pydump.f90` existed at all, zero calls anywhere):
two extra unconditional `CALL` instructions in `eqdyna3d.o`/`driver.o`
(into an empty subroutine in the default build). This is not literal
byte-identity to pristine. It IS verified functionally neutral:

```
$ mpirun -np 1 -wdir <tpv8 case> <default-build-with-noop-stub>/eqdyna
$ mpirun -np 1 -wdir <same case> <PYDUMP=1 build>/eqdyna-pydump
$ cmp <case>/frt.txt0-from-default <case>/frt.txt0-from-pydump-build
# exit 0 -- byte-identical simulation output
```

This was run for tpv8 (114 steps, serial) during this tier's development:
`frt.txt0` from the default (no-op stub) binary and from the `PYDUMP=1`
binary are byte-identical. `testsys/parity/fixtures/eqdyna-default-for-
neutrality-check` (written by `make_fixtures.py`) is the default-build
binary saved specifically so this comparison can be re-run by hand or
folded into a future regression test.

## Tier contents

- `make_fixtures.py` — builds `eqdyna-pydump`, sets up a serial
  (`nx=ny=nz=1`) tpv8 case, runs it, leaves `frt.txt0` (golden oracle) +
  `pydump_*` (static-state + step 1-5 checkpoint dumps) in
  `testsys/parity/fixtures/test_tpv8_serial/`. Fixtures are regenerable,
  not committed as static blobs — this script IS their provenance.
- `run_parity.py` — loads the fixtures via `python/eqdyna/port.py`'s
  `load()`, runs the NumPy port and (if importable) the JAX port for the
  full fixture step count, and FAILS if any output column exceeds the
  thresholds documented in its own module docstring (2-5x the actual
  observed max abs diff from README-parity.md Updates 1-3/6). Columns the
  Fortran never writes for friclaw=1 are gated at exactly 0.
- Wired into `testsys/run.py` as `parity` (opt-in, not part of `all` —
  needs a Fortran build + fixtures a fresh checkout doesn't have).

## Running it

```
python3 testsys/parity/make_fixtures.py     # once, or whenever src/*.f90 changes
python3 testsys/run.py parity
```
