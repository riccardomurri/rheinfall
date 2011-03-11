#! /bin/bash

flavor="$1"; shift
if expr match "$flavor" 'run/' >/dev/null; then
    flavor=$(echo $flavor | cut -c5-)
fi

bindir="run/${flavor}"
if test ! -d "$bindir"; then
    echo 1>&2 "Cannot read directory: '$bindir'"
    exit 1
fi

libdir=$(pwd)/util/schroedinger
source $libdir/functions.sh \
    || { echo 1>&2 "Cannot load 'functions.sh' - aborting."; exit 1; }
alias module=true # work around lack of "modules" on idesl4-g
set_mpi_and_compiler_flavor "$flavor"

flavor="${rev}-${compiler}-${mpi}"

for prog in \
    rank-int \
    rank-mod \
    crank-int \
    ; 
    #rank-int-omp \
    #rank-mod-omp \
do
    if test ! -x "${bindir}/${prog}"; then
        echo 1>&2 "Skipping executable '$prog': not found in '$bindir'"
        continue
    fi
    qsub \
        -cwd \
        -S /bin/bash \
        -N "${prog}.${flavor}.single" \
        -l s_cpu=$(( 24 * 3600 )) \
        -pe openmpi2 8 \
        -j y \
        run.sh "${bindir}/${prog}" \
        inputs/M05-D{4,5,6,7,8}.txt \
        inputs/M06-D{5,6,7,8,11,10}.txt \
        ;
done
