#! /bin/bash
#
PROG="$(basename $0)"

usage () {
cat <<EOF
Usage: $PROG [-r /path/to/rank.sh] [-x /path/to/rank-program] DIRECTORY

Submit a SLURM job (using rank.sh) for each '.sms' matrix file found
in DIRECTORY (descend it recursively).

Options:

  --help, -h  Print this help text.

  -r PATH     Path to the 'rank.sh' script to use.
              (Default: ./rank.sh)

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

rank_sh='rank.sh'
exe='src.c++/rank-int'

if [ "x$(getopt -T)" == 'x--' ]; then
    # old-style getopt, use compatibility syntax
    set -- $(getopt 'hr:x:' "$@")
else
    # GNU getopt
    set -- $(getopt --shell sh -l 'help' -o 'hr:x:' -- "$@")
fi
while [ $# -gt 0 ]; do
    case "$1" in
        -r) shift; eval rank_sh=$1 ;;
        -x) shift; eval exe=$1 ;;
        --help|-h) usage; exit 0 ;;
        --) shift; break ;;
    esac
    shift
done

eval dir=$1

## consistency check

require_command find
require_command sbatch

if ! test -e "$rank_sh"; then 
    die 126 "Script file '$rank_sh' does not exist (use the '-r' option to set a different path)."
fi
if ! test -r "$rank_sh"; then 
    die 126 "Cannot read script file '$rank_sh' (use the '-r' option to set a different path)."
fi
if ! test -x "$rank_sh"; then 
    die 126 "Cannot execute script '$rank_sh' (use the '-r' option to set a different path)."
fi

if ! test -e "$exe"; then 
    die 126 "Program file '$exe' does not exist (use the '-x' option to set a different path)."
fi
if ! test -r "$exe"; then 
    die 126 "Cannot read program file '$exe' (use the '-x' option to set a different path)."
fi
if ! test -x "$exe"; then 
    die 126 "Cannot execute program '$exe' (use the '-x' option to set a different path)."
fi

if ! test -d "$dir"; then 
    die 125 "Path '$dir' is not a directory."
fi
if ! test -r "$dir"; then 
    die 125 "Cannot read directory '$dir'."
fi


## main

# disable core dumps
ulimit -c 0

# submit them all
find "$dir" -name '*.sms.gz' -or -name '*.sms' \
  | (while read file; do
      sbatch -N 1 --mail-type=FAIL \
          --job-name=$(basename "$file" .sms) \
          "$rank_sh" -x "$exe" "$file"
      done)
