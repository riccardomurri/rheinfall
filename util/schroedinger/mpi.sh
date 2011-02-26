#! /bin/bash

flavor="$1"; shift

rev="$(echo $flavor | cut -d/ -f1)"
compiler_and_mpilib="$(echo $flavor | cut -d/ -f2)"
if [ -z "$compiler_and_mpilib" ]; then
    compiler_and_mpilib="$rev"
    rev="$( cat .bzr/branch/last-revision | (read revno rest; echo r$revno) )"
fi

compiler=$(echo $compiler_and_mpilib | cut -d- -f1)
if [ -z "$compiler" ]; then
    compiler=gcc443
fi
mpi="$(echo $compiler_and_mpilib | cut -d- -f2)"
if [ -z "$mpi" ]; then
    mpi=ompi
fi

bindir="run/${flavor}"
if test ! -d "$bindir"; then
    echo 1>&2 "Cannot read directory: '$bindir'"
    exit 1
fi

flavor="${rev}-${compiler}-${mpi}"

nps='16 32 64' # 128 256 512'
ws='4096 1024 256 64 16' # 4 1'

for np in $nps; do
for w in $ws; do
for prog in \
    rank-int-mpi \
    rank-mod-mpi \
    rank-int-mpi-omp \
    rank-mod-mpi-omp \
    ; 
do
    if test ! -x "${bindir}/${prog}"; then
        echo 1>&2 "Skipping executable '$prog': not found in '$bindir'"
        continue
    fi
    qsub \
        -cwd \
        -S /bin/bash \
        -N "mpi-n${np}-w${w}.${prog}.${flavor}" \
        -l s_rt=$(( 24 * 3600 )) \
        -pe openmpi2 $np \
        -j y \
        run.sh "${bindir}/${prog}" \
        -w $w \
        inputs/M05-D{4,5,6,7,8}.txt \
        inputs/M06-D{5,6,7,8,11,10}.txt \
        ;
done
done
done
