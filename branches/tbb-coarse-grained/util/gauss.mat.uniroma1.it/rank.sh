#! /bin/bash
#
PROG="$(basename $0)"

usage () {
cat <<EOF
Usage: $PROG [-x /path/to/rank-program] MATRIX.sms

Compute the rank of the MATRIX (an SMS format file)
using the Rheinfall 'rank-*' program.

Options:

  --help, -h  Print this help text.

  -x PATH     Path to the 'rank-*' program to use.
              (default: src.c++/rank-int)

EOF
}


## helper functions
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

is_absolute_path () {
    expr match "$1" '/' >/dev/null 2>/dev/null
}


## parse command-line 

exe='src.c++/rank-int'

if [ "x$(getopt -T)" == 'x--' ]; then
    # old-style getopt, use compatibility syntax
    set -- $(getopt 'hx:' "$@")
else
    # GNU getopt
    set -- $(getopt --shell sh -l 'help' -o 'hx:' -- "$@")
fi
while [ $# -gt 0 ]; do
    case "$1" in
        -x) shift; eval exe=$1 ;;
        --help|-h) usage; exit 0 ;;
        --) shift; break ;;
    esac
    shift
done

eval matrix=$1

## consistency check

if ! test -e "$exe"; then 
    die 126 "Program file '$exe' does not exist (use the '-x' option to set a different path)."
fi
if ! test -r "$exe"; then 
    die 126 "Cannot read program file '$exe' (use the '-x' option to set a different path)."
fi
if ! test -x "$exe"; then 
    die 126 "Cannot execute program '$exe' (use the '-x' option to set a different path)."
fi

if ! test -e "$matrix"; then 
    die 125 "Matrix file '$matrix' does not exist."
fi
if ! test -r "$matrix"; then 
    die 125 "Cannot read matrix file '$matrix'."
fi


## main

# exit on first error
set -e

stem=$(echo "$matrix" | tr '/' '~' | sed -e 's/.sms.gz//')

# uncompress matrix to a temporary file...
wkf="${stem}.sms"
gzip -dc < "$matrix" > "$wkf"

# differentiate output file by matrix name; need to keep full path in
# here, because matrices in different collections may have the same
# base name
output=${stem}.${SLURM_JOBID:-out}.txt

# disable core dumps
ulimit -c 0

# use max 30GB of memory
ulimit -m $(expr 30 '*' 1024 '*' 1024)
ulimit -v $(expr 30 '*' 1024 '*' 1024)

# use Google's TCMalloc library
export LD_PRELOAD=sw/lib/libtcmalloc.so

# run!
srun -N 1 --propagate --time=$(expr 48 '*' 60) \
    --cpu_bind=verbose,cores --mem_bind=verbose,local \
    -o "$output" $exe "$wkf"
rc=$?

# remove temp file
rm -f "$wkf"

# exit with srun return code
exit $rc
