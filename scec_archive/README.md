# scec_archive

EQdyna's submissions to the [SCEC/USGS Spontaneous Rupture Code Verification
Project](https://strike.scec.org/cvws/), archived in the repository so a run's
published standing travels with the code that produced it.

## Layout

```
scec_archive/<benchmark>/<code-version>-<resolution>-<year>/
```
e.g. `scec_archive/tpv29/eqdyna-v3.1-100m-2015/`, holding the files exactly as
uploaded (raw SCEC format, no HTML wrapper): `cplot` (rupture-time contours),
`faultst*` / `body*` (on- and off-fault time series), and a `PROVENANCE.md`
naming the source, fetch date, and the submission label in the public record.

`INDEX.md` lists everything archived — one row per submission — together with
the enumeration it came from and the submissions the portal is known to hold
but will not serve. `scratch/scec_archive/` holds the tooling that produced it.

This directory is **gitignored** — it holds bulky historical run data, not
source. It is referenced as a historical reference by the repo (compset
READMEs, pathway_forward.md, and the TPV29 scoring tooling in
`scratch/tpv29/scoring/`), so keep the layout below stable even though git
does not track it. It is NOT reproducible from the repo: back it up
independently.

## Rules

- These are **archival records, read-only** (PROJECT_RULES.md rule 7). Never
  edit a file in place; a new run is a new directory.
- Every directory states the code version and resolution in its name, so a
  comparison never has to guess which EQdyna produced it.
- Public cross-code metrics and other codes' data are NOT mirrored here; they
  are fetched on demand (see `scratch/tpv29/scoring/`).
