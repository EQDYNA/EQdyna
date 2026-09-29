#! /bin/bash
# Build EQdyna (src/fortran -> bin/eqdyna) and set its environment.
#
# Machines: ls6 (TACC Lonestar6), grace (TAMU Grace), ubuntu (22.04), macos.
#
#   ./install-eqdyna.sh -m <machine>   build
#   ./install-eqdyna.sh -e <machine>   install dependencies, then build (ubuntu, macos, ls6)
#   ./install-eqdyna.sh -c <machine>   set up the environment only, no build
#   source install-eqdyna.sh -c <m>    set modules, venv, EQDYNAROOT and PATH in this shell

while getopts "m:e:c:h" OPTION; do
    case $OPTION in
        m) MACH=$OPTARG ;;
        e) MACH=$OPTARG; ENV="True" ;;
        c) MACH=$OPTARG; CONFIG="True" ;;
        h) sed -n '2,10p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//' ;;
    esac
done

if [ -n "$MACH" ]; then
    export MACHINE=$MACH
    case $MACHINE in
        ls6)
            module load netcdf/4.6.2
            echo "netCDF: $TACC_NETCDF_INC $TACC_NETCDF_LIB"
            # Python (jax, mpi4py) lives in a venv on $WORK; -e creates it.
            VENV=${EQDYNA_VENV:-$WORK/eqdyna-venv}
            if [ -n "$ENV" ]; then
                module load python
                if ! python3 -c 'import sys; sys.exit(sys.version_info < (3, 10))'; then
                    echo "install-eqdyna.sh: $(python3 -V) is too old; jax needs 3.10+ (module spider python)" >&2
                    return 1 2>/dev/null || exit 1
                fi
                python3 -m venv "$VENV"
                "$VENV/bin/pip" install --upgrade pip
                "$VENV/bin/pip" install numpy scipy netCDF4 matplotlib xarray pytest jax
                # mpi4py must be built against the loaded MPI (Intel MPI), not a wheel's.
                MPICC=mpicc "$VENV/bin/pip" install --no-binary=mpi4py --no-cache-dir mpi4py
            fi
            if [ -f "$VENV/bin/activate" ]; then
                module load python
                source "$VENV/bin/activate"
            fi ;;
        grace)
            module load netCDF
            echo "netCDF: ${EBROOTNETCDF}/include ${EBROOTNETCDF}/lib64" ;;
        ubuntu)
            export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH
            if [ -n "$ENV" ]; then
                # gfortran is listed explicitly: mpich pulls in only its
                # runtime, and a bare ubuntu:22.04 has no compiler (2026-09-16).
                apt-get install git vim make mpich gfortran
                apt-get install libnetcdf-dev libnetcdff-dev
                apt-get install python3 python3-pip
                pip install numpy netCDF4 matplotlib xarray
                pip install --upgrade numpy
            fi ;;
        macos)
            export MACOS_NETCDF_INC=$(brew --prefix netcdf-fortran)/include
            export MACOS_NETCDF_LIB=$(brew --prefix netcdf)/lib
            export MACOS_NETCDFF_LIB=$(brew --prefix netcdf-fortran)/lib
            if [ -n "$ENV" ]; then
                brew install mpich python
                pip3 install --break-system-packages numpy netCDF4 matplotlib xarray
            fi ;;
        *)
            echo "install-eqdyna.sh: unknown machine '$MACHINE' (ls6, grace, ubuntu, macos)" >&2
            return 1 2>/dev/null || exit 1 ;;
    esac

    if [ -z "$CONFIG" ]; then
        (cd src/fortran && make) || { echo "install-eqdyna.sh: make failed" >&2; return 1 2>/dev/null || exit 1; }
        mkdir -p bin
        mv src/fortran/eqdyna bin
    fi

    # Activate the tracked git hooks (rules 21/21b). Relative, so every
    # worktree uses its own copy; idempotent.
    if git rev-parse --git-dir > /dev/null 2>&1; then
        git config core.hooksPath testsys/hooks
    fi

    # Only the entry points run directly, never the whole directory (rule 13).
    chmod 755 scripts/case.setup scripts/clean.py scripts/create.newcase \
        scripts/generateFaultInterface scripts/plotRuptureDynamics \
        scripts/plotSlipAndRPT
fi

export EQDYNAROOT=$(pwd)
export PATH=$(pwd)/bin:$(pwd)/scripts:$PATH
echo "EQDYNAROOT=$EQDYNAROOT"
