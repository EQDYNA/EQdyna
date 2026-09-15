"""`python3 -m eqdyna <case_dir> [nsteps]` -- thin wrapper around
eqdyna3d.main() so the package is directly runnable, the same way
`./eqdyna` runs the Fortran binary."""
from eqdyna.eqdyna3d import main

if __name__ == '__main__':
    main()
