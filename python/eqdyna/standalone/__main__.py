"""`python3 -m eqdyna.standalone <case_dir> [nsteps]` entry point -- thin
wrapper around main.main() so the package is directly runnable."""
from eqdyna.standalone.main import main

if __name__ == '__main__':
    main()
