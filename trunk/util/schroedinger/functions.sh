#! /bin/sh
#
# This is a collection of functions used to carry out common tasks 
# in Schroedinger helper scripts.  It won't run as a stand-alone script.
#

set_mpi_and_compiler_flavor () {
    flavor="$1"; shift
    compiler="$(echo $flavor | (IFS=- read one two three; echo $one) )"
    if [ -z "$compiler" ]; then
        compiler=gcc443
    fi
    mpi="$(echo $flavor |  (IFS=- read one two three; echo $two) )"
    if [ -z "$mpi" ]; then
        mpi=ompi
    fi

    ## load modules
    source /panfs/panfs0.ften.es.hpcn.uzh.ch/share/software/Modules/default/init/sh

    supported_compilers='gcc412 gcc441 gcc443 gcc450 icc'
    supported_mpilibs='openmpi mvapich intel none'

     # load MPI - must match what the binary was compiled with!
    case "$mpi" in
        ompi|openmpi) 
            case "$compiler" in
                gcc450)    module load mpi/openmpi/gcc-4.5.0 ;;
                gcc*)      module load mpi/openmpi/gcc ;;
                icc|intel) module load mpi/openmpi/intel ;;
            esac
            ;;
        ompi*) # My own OpenMPI install in ${sw}/lib etc.
            ;;
        mvapich) # MVAPICH is not supported by modules, apparently
            case "$compiler" in
                gcc412)
                    export PATH=/panfs/panfs0.ften.es.hpcn.uzh.ch/mpi-libs/gcc/mvapich2-1.4rc2/bin:$PATH
                    export LD_LIBRARY_PATH=/panfs/panfs0.ften.es.hpcn.uzh.ch/mpi-libs/gcc/mvapich2-1.4rc2/lib:$LD_LIBRARY_PATH
                    export LD_RUN_PATH=$LD_LIBRARY_PATH
                    ;;
                gcc*)
                    export PATH=/panfs/panfs0.ften.es.hpcn.uzh.ch/mpi-libs/gcc-4.4.1/mvapich2-1.4rc2/bin:$PATH
                    export LD_LIBRARY_PATH=/panfs/panfs0.ften.es.hpcn.uzh.ch/mpi-libs/gcc-4.4.1/mvapich2-1.4rc2/lib:$LD_LIBRARY_PATH 
                    export LD_RUN_PATH=$LD_LIBRARY_PATH
                    ;;
            esac
            ;;
        mvapich*) # My own MVAPICH2 install in ${sw}/lib etc.
            ;;
        impi|intel) # Intel MPI is not supported by modules, apparently
            export PATH=/panfs/panfs0.ften.es.hpcn.uzh.ch/mpi-libs/intel_mpi/3.2.1.009/bin64:$PATH
            export LD_LIBRARY_PATH=/panfs/panfs0.ften.es.hpcn.uzh.ch/mpi-libs/intel_mpi/3.2.1.009/lib64:$LD_LIBRARY_PATH
            export LD_RUN_PATH=$LD_LIBRARY_PATH
            ;;
        none) # no MPI
            while type -a mpicc >/dev/null 2>&1; do
                p=$(dirname $(which mpicc) )
                echo 1>&2 "MPI compiler 'mpicc' found in PATH at '$p'; removing since no MPI was requested ..."
                export PATH=$(echo $PATH | sed -e "s|:${p}||;s|${p}:||;")
            done
            ;;
        *) 
            die 1 "Unknown MPI library '${mpi}' - please choose one of: $supported_mpilibs"
        ;;
    esac

    # load the compiler and libs
    module load binutils/2.20.1
    case "$compiler" in
        gcc412)
            module load gcc/4.1.2
            export CC=gcc
            export CXX=g++
            export OMPI_CC=`which gcc`
            export OMPI_CXX=`which g++`
            cflags='-march=nocona'
            cxxflags='-march=nocona'
            std_cflags='-O3'
            std_cxxflags='-O3'
            toolset=gcc
            ;;
        gcc441)
            module load gcc/4.4.1
            export CC=gcc-4.4.1
            export CXX=g++-4.4.1
            export OMPI_CC=`which gcc-4.4.1`
            export OMPI_CXX=`which g++-4.4.1`
            cflags='-march=nocona'
            cxxflags='-march=nocona'
            std_cflags='-O3'
            std_cxxflags='-O3'
            toolset=gcc
            ;;
        gcc443)
            module load gcc/4.4.3
            export CC=gcc
            export CXX=g++
            export OMPI_CC=`which gcc`
            export OMPI_CXX=`which g++`
            cflags='-march=native'
            cxxflags='-march=native'
            std_cflags='-O3'
            std_cxxflags='-O3'
            toolset=gcc
            ;;
        gcc450)
            module load gcc/4.5.0
            export CC=gcc
            export CXX=g++
            export OMPI_CC=`which gcc`
            export OMPI_CXX=`which g++`
            cflags='-march=native'
            cxxflags='-march=native'
            std_cflags='-O3'
            std_cxxflags='-O3'
            toolset=gcc
            ;;
        icc|intel)
            module load intel/comp/11.1.064
            export CC=icc
            export CXX=icpc
            export OMPI_CC=`which icc`
            export OMPI_CXX=`which icpc`
            export I_MPI_CC=`which icc`
            export I_MPI_CXX=`which icpc`
            cflags='-xHOST -O3 -ipo -no-prec-div'
            cxxflags='-xHOST -O3 -ipo -no-prec-div'
            std_cflags='-O3'
            std_cxxflags='-O3'
            toolset=intel
            ;;
        *) 
            die 1 "Unknown compiler flavor '${compiler}' - please choose one of: $supported_compilers"
            ;;
    esac

    # enable local sw
    sw="$HOME/sw/${compiler}-${mpi}"
    PATH=${sw}/bin:$PATH; export PATH
    LD_LIBRARY_PATH=${sw}/lib:$LD_LIBRARY_PATH; export LD_LIBRARY_PATH
}


show_build_information () {
    echo === running info ===
    echo flavor: $flavor
    echo compiler: $compiler
    echo mpi: $mpi
    echo node: $(hostname)
    echo === compiler information ===

    which ${CXX}
    set -e # exit on error here, in case the compiler does not have enough licences
    ${CXX} --version
    set +e

    if [ "x${mpi}" != 'xnone' ]; then
        which mpicxx
        mpicxx --version
        echo
    fi

    which as
    as --version
    echo
}


prepare_openmpi_environment () {
    require_environment_var PE_HOSTFILE
    require_environment_var TMPDIR
    
    cat $PE_HOSTFILE | \
        (while read hostname nslots queue rest; 
            do echo "$hostname slots=$nslots"; done) > $TMPDIR/hostfile

    set -e
    eval `ssh-agent -s`
    ssh-add $HOME/.ssh/id_dsa
    set +e
   
    gmca="--gmca plm_rsh_agent $HOME/bin/qrsh.sh"
    hostfile="--hostfile $TMPDIR/hostfile"
}


## generic functions
die () {
  rc="$1"
  shift
  (echo -n "$PROG: ERROR: ";
      if [ $# -gt 0 ]; then echo "$@"; else cat; fi) 1>&2
  exit $rc
}

have_command () {
  type "$1" >/dev/null 2>/dev/null
}

require_command () {
  if ! have_command "$1"; then
    die 1 "Could not find required command '$1' in system PATH. Aborting."
  fi
}

require_environment_var () {
    if [ -z "$(eval echo '$'$1)" ]; then
        die 1 "Require environment variable '$1' is not defined or has empty value. Aborting."
    fi
}

is_absolute_path () {
    expr match "$1" '/' >/dev/null 2>/dev/null
}

_ () {
    echo
    echo ==== "$@" ...;
}

