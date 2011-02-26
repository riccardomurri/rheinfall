#! /usr/bin/env python
#
"""
Import run data from irank* output.
"""
__docformat__ = 'reStructuredText'

import logging
from optparse import OptionParser
import os
import os.path
import re
import sys

try:
    from collections import namedtuple
except ImportError:
    from namedtuple import namedtuple


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




def nullable(fn):
    def fn_or_None(x):
        if x is None or x == '':
            return None
        else:
            return fn(x)
    return fn_or_None
    
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

rank_line_re = re.compile(r'^[^ ]+/[idqz]?rank(-int|-mod|-double)?([_-]mpi|[_-]omp){0,2} file:')
jobid_re = re.compile(r'.[eo]([0-9]+)')
whitespace_re = re.compile(r'\s+', re.X)
garbage_re = re.compile(r'(MPI process \(rank: \d+\) terminated unexpectedly|\[[0-9a-z]+:\d+\] \[ *\d+\]).*')

def process_input_file(input_file_name, keys):
    result = { }
    input_file = open(input_file_name, 'r')
    logger.log(1, "Processing file %s ...", input_file_name)
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
            if k in schema:
                try:
                    outcome[k] = schema[k](outcome[k])
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

if len(args) == 2:
    set1 = [ args[0] ]
    set2 = [ args[1] ]
else:
    # ':' separates the two sets of files
    brk = args.index(':')
    set1 = args[:brk]
    set2 = args[(brk+1):]

logger.log(1, "Comparing: %s (length %d)" 
           % (str.join(",", options.compare), len(options.compare)))
result1 = process_set(set1, options.compare)
result2 = process_set(set2, options.compare)

# build list of processed matrices
matrices = set()
for m in result1:
    matrices.add(m)
for m in result2:
    matrices.add(m)


def center(text, width):
    if len(text) > width:
        return text
    else:
        padding = ((width - len(text)) / 2) * " "
        return  padding + text + padding

hdr = ("=" * 16) + " " + (" " + ("=" * 8) + " " + ("=" * 8) + " " + ("=" * 7) + " ") * len(options.compare)
print(hdr)
print(("%-16s " % center("Matrix", 16))
      + str.join('', [(" %25s " % center(k, 25)) for k in options.compare]))
print((16 * '-') + " " + (" " + (25 * '-') + " ") * len(options.compare))
print((16 * " ") + " " 
      + (" %8s %8s %7s " % (center("A", 8), center("B", 8), center("var%", 7))) * len(options.compare))
print(hdr)

def has_valid_value(d, k):
    return d.has_key(k) and d[k] is not None

def variation(v1, v2):
    if v1 == 0 or v2 == 0:
        return ""
    elif v1 == v2:
        return ""
    else:
        return "%+.2f" % (100 * ((v2 / v1) - 1.0))

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


for m in sorted(matrices, cmp=cmp_numerically):
    if has_valid_value(result1, m) and has_valid_value(result2, m):
        r1 = result1[m]
        r2 = result2[m]
        values = [ m ]
        for k in options.compare:
            if has_valid_value(r1, k) and has_valid_value(r2, k):
                values += ["%.2f" % r1[k], 
                           "%.2f" % r2[k], 
                           variation(r1[k], r2[k])]
            elif has_valid_value(r1, k):
                values += ["%.2f" % r1[k], "***", ""]
            elif has_valid_value(r2, k):
                values += ["***", "%.2f" % r2[k], ""]
    elif has_valid_value(result1, m):
        values = [ m ]
        r1 = result1[m]
        for k in options.compare:
            values += ["%.2f" % r1[k], "***", ""]
    elif has_valid_value(result2, m):
        values = [ m ]
        r2 = result2[m]
        for k in options.compare:
            values += ["***", "%.2f" % r2[k], ""]
    print(("%-16s " + (" %8s %8s %7s ") * len(options.compare)) % tuple(values))

print(hdr)
