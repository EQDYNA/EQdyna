#! /bin/bash 

# The shell script is to set up environments for EQdyna and 
#	install it. It will call the makefile inside src/ and generate 
#	an executable eqdyna and move it to bin/.

# Currently, the machines supported are:
#	ls6:	Lonestar6 at TACC
#	ubuntu: Ubuntu 22.04

# Usage: install-eqdyna.sh [-h] [-m Machine_name] [-c Machine_name]

while getopts "m:c:h" OPTION; do
    case $OPTION in 
        m)
            MACH=$OPTARG
            ;;
        c)
            MACH=$OPTARG
            CONFIG="True"
            ;;
        h)
            echo "Usage: ./install-eqdyna.sh [-h] [-m Machine_name] [-c Machine_name] "
            echo "                                                                     "
            echo "Examples:                                                            "
            echo "                                                                     "
            echo "./install-eqdyna.sh -h                                              "
            echo " -----Display this help message                                      "
            echo "                                                                     "
            echo "./install-eqdyna.sh -m ls6                                          "
            echo " -----Install EQdyna on Lonestar6 at TACC                           "
            echo "                                                                     "
            echo "./install-eqdyna.sh -c ubuntu                                       "
            echo " -----Simply set up envs for EQdyna without installation            "
            echo " -----on ubuntu                                                      "
            echo "                                                                     "
            echo "source install-eqdyna.sh                                            "
            echo " -----Activate ENV VAR EQQUASIROOT and add exes to PATH              "
            echo "                                                                     "
            echo "Currently supported machines include:                                "
            echo " ls6/ubuntu                                                          "
            ;;
    esac
done 

if [ -n "$MACH" ]; then 
    export MACHINE=$MACH
    if [ $MACHINE == "ls6" ]; then 
        echo "Installing EQdyna on Lonestar6 at TACC ... ..."
        
        echo "Loading netcdf module ... ..."
        module load netcdf 
        ml
        
        echo "NETCDF INC and LIB PATH"
        echo $TACC_NETCDF_INC
        echo $TACC_NETCDF_LIB
        
    elif [ $MACHINE == "ubuntu" ]; then 
        echo "Installing EQdyna on Ubuntu 22.04 ... ..."
        export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH
    fi 
    
    if [ -n "$CONFIG" ]; then 
        echo "Simply configure EQdyna without installation ... ..."
    else
        cd src
        make
        cd ..
        mkdir bin
        mv src/eqdyna bin
    fi

    export EQDYNAROOT=$(pwd)
    export PATH=$(pwd)/bin:$PATH
    export PATH=$(pwd)/scripts:$PATH
    
    chmod -R 755 scripts
fi

export EQDYNAROOT=$(pwd)
export PATH=$(pwd)/bin:$PATH
export PATH=$(pwd)/scripts:$PATH

echo EQDYNAROOT
echo PATH 
