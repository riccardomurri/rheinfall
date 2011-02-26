#! /usr/bin/env python
#
"""
Import run data from rank*/lbrank* output.
"""
__docformat__ = 'reStructuredText'

import logging
from optparse import OptionParser
import os
import os.path
import re
import sys

# patch sys.path to include likely locations for the `namedtuple` and
# `texttable` modules
p = os.path.dirname(sys.argv[0])
sys.path.append(p)
sys.path.append(os.path.join(p, 'util'))

try:
    from collections import namedtuple
except ImportError:
    from namedtuple import namedtuple
from texttable import Texttable


cmdline = OptionParser()
cmdline.add_option('-c', '--compare', dest='compare', action='append', metavar='FIELD', default=[],
                   help="Print comparison of the specified fields."
                   " Repeat multiple times to compare several fields.")
cmdline.add_option('-D', '--defaults', dest='defaults', action='append', metavar="SPEC",
                   help="Provide defaults for inserted fields, in the form"
                   " of a comma-separated list of KEY:VALUE pairs.")
cmdline.add_option('-v', '--verbose', dest='verbose', action='count', default=0,
                   help="Report more details on program execution.")
(options, args) = cmdline.parse_args()

if len(options.compare) == 0:
    options.compare = ['cputime', 'wctime']


class Logger(object):
    def __init__(self, verbose, progname=None):
        if progname is None:
            self._logger = logging.getLogger()
        else:
            self._logger = logging.getLogger(progname)
        self._logger.setLevel(logging.WARNING - verbose)
    def log(self, lvl, *args, **kw):
        self._logger.log(logging.WARNING - lvl, *args, **kw)
logging.basicConfig(level=logging.DEBUG,
                    format='%(name)s: %(message)s')
logger = Logger(options.verbose, os.path.basename(sys.argv[0]))


# build record defaults
defaults = { 
    'compiler':None,
    'mpi':None,
    'mpilib':None,
    'omp':None,
    'revno':None,
    'w':1,
    }
if options.defaults is not None:
    for kv in str.join(',', options.defaults).split(','):
        k, v = kv.split(':')
        defaults[k] = v




def nullable(ctor):
    def ctor_or_None(x):
        if x is None or x == '':
            return None
        else:
            return ctor(x)
    return ctor_or_None
    
schema = {
    'exitcode':   int,            # UNIX exit code
    'matrix':     str,            # matrix file processed
    'rows':       int,            # no. of matrix rows
    'cols':       int,            # no. of matrix columns
    'nonzero':    int,            # no. of nonzero entries
    'rank':       int,            # computed matrix rank
    'cputime':    float,          # CPU time consumed on MPI rank 0
    'wctime':     float,          # wall-clock time on MPI rank 0
    'omp':        nullable(int),  # no. of OpenMP threads (0 = no OpenMP)
    'mpi':        nullable(int),  # no. of MPI ranks (0 = no MPI)
    'w':          int,            # width of a band of columns
    'compiler':   nullable(str),  # compiler used for compiling the C++ source
    'mpilib':     nullable(str),  # MPI library used 
    'revno':      nullable(int),  # BZR revision no. of the C++ source
    'provenance': str,            # file name where this info was extracted from
}

# these have been renamed over time
translate = {
    'time': 'cputime',
}

rank_line_re = re.compile(r'^[^ ]+/(?:[idqz]?rank(-int|-mod|-double)?([_-]mpi|[_-]omp){0,2}|lbrank_(?:bb|se_linear|se_none)) file:')
jobid_re = re.compile(r'.[eo]([0-9]+)')
whitespace_re = re.compile(r'\s+', re.X)
garbage_re = re.compile(r'(MPI process \(rank: \d+\) terminated unexpectedly|\[[0-9a-z]+:\d+\] \[ *\d+\]|.*SIG[A-Z]+[0-9]*|/[/a-z0-9_]+: line [0-9]+: [0-9]+ Aborted|\[[a-z0-9]+:[0-9]+\] \*\*\* Process received signal).*')

def process_input_file(input_file_name, keys):
    result = { }
    input_file = open(input_file_name, 'r')
    logger.log(2, "Processing file %s ...", input_file_name)
    for line in input_file.readlines():
        line = line.strip()
        if not rank_line_re.match(line):
            continue
        line = garbage_re.sub('', line)
        words = line.split(' ')[1:]
        outcome = dict((k,v) for k,v in defaults.items() if k in keys)
        outcome.update(dict(w.split(':',1) for w in words))
        if not outcome.has_key('rank'):
            # computation failed, go to next line
            continue
        # extract matrix name from file name
        outcome['matrix'] = os.path.splitext(os.path.basename(outcome['file']))[0]
        # format according to schema
        outcome_is_valid = True
        for k in outcome:
            if k in translate:
                k_ = translate[k]
            else:
                k_ = k
            if k_ in schema:
                try:
                    outcome[k_] = schema[k_](outcome[k])
                    if k_ != k:
                        del outcome[k]
                except:
                    # invalid input
                    outcome_is_valid = False
                    break
        if outcome_is_valid:
            result[outcome['matrix']] = outcome
    return result


def process_set(input_files, keys):
    result = None
    for input_file in input_files:
        outcome = process_input_file(input_file, keys)
        if result is None:
            result = outcome
        else:
            for m in outcome:
                for k in keys:
                    result[m][k] = min(result[m][k], outcome[m][k])
    return result


## main 

# ':' separates sets of files
if ':' in args:
    filesets = []
    while ':' in args:
        brk = args.index(':')
        filesets.append(args[:brk])
        del args[:brk]
    filesets.append(args)
else:
    # each file is a set
    filesets = [ [filename] for filename in args ]

# process each fileset
logger.log(2, "Comparing: %s (length %d)" 
           % (str.join(",", options.compare), len(options.compare)))
results = [ process_set(fileset, options.compare) for fileset in filesets ]

# build list of processed matrices
matrices = set()
for result in results:
    for m in result:
        matrices.add(m)

# create output table
t = Texttable(0)
t.set_deco(Texttable.BORDER | Texttable.HEADER)
hdr = ["Matrix"]
for feature in options.compare:
    if 2 == len(filesets):
        hdr += [feature, "", "(var%)"]
    else:
        hdr_ = [feature] + ["..."] * (len(filesets)-1)
        for n in range(len(hdr_)):
            hdr_[n] += "(%d)" % (n+1)
        hdr += hdr_
t.header(hdr)
t.set_cols_align(['l'] + (['l'] * max(len(filesets), 3) * len(options.compare)))
t.set_cols_align(['m'] + (['m'] * max(len(filesets), 3) * len(options.compare)))


# hack to get M06-D10 come after M06-D9
digits_re = re.compile(r'(\d+)')
def cmp_numerically(a,b):
    parts_a = digits_re.split(a)
    for pos in range(len(parts_a)):
        try:
            parts_a[pos] = int(parts_a[pos])
        except ValueError:
            pass
    parts_b = digits_re.split(b)
    for pos in range(len(parts_b)):
        try:
            parts_b[pos] = int(parts_b[pos])
        except ValueError:
            pass
    return cmp(parts_a, parts_b)


def variation(v1, v2, fmt="%+.2f"):
    if v1 == 0 or v2 == 0:
        return ""
    elif v1 == v2:
        return ""
    else:
        return fmt % (100 * ((v2 / v1) - 1.0))

def has_valid_value(d, k):
    return d.has_key(k) and d[k] is not None

def value(d, k, fmt="%.2f"):
    if has_valid_value(d, k):
        return fmt % d[k]
    else:
        return "***"

def pprint_values(m, k, rs):
    all_valid_values = True
    values = [ ]
    for r in rs:
        values.append(value(r[m], k))
        if not has_valid_value(r[m], k):
            all_valid_values = False
    if len(values) == 2:
        if all_valid_values:
            values.append(variation(rs[0][m][k], rs[1][m][k]))
        else:
            values.append("")
    return values

for m in sorted(matrices, cmp=cmp_numerically):
    for result in results:
        if not has_valid_value(result, m):
            result[m] = dict() # simulate empty result set
    values = [ m ]
    for k in options.compare:
        values += pprint_values(m, k, results)
    t.add_row(values)

print(t.draw())
if options.verbose > 0:
    for n, fileset in enumerate(filesets):
        print ("    (%d): %s" % (n+1, str.join(' ', fileset)))
