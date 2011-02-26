#! /usr/bin/env python
#
"""
Compare a list of matrices with run data from rank*/lbrank* output 
and output the list of matrices for which running time is not known.
"""
__docformat__ = 'reStructuredText'

import logging
import math
from optparse import OptionParser
import os
import os.path
import re
import sys


cmdline = OptionParser()
cmdline.add_option('-d', '--duration', dest='duration', action='store', default="140000",
                   help="Set duration to fake for missing matrices (default: %default)")
cmdline.add_option('-v', '--verbose', dest='verbose', action='count', default=0,
                   help="Report more details on program execution.")
(options, args) = cmdline.parse_args()

try:
    options.duration = float(options.duration)
    if options.duration <= 0:
        raise ValueError(str(options.duration))
except:
    sys.stderr.write("Argument to option -d/--duration must be a positive number.")

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


## parse command line

if ':' in args:
    i = args.index(':')
    list_files = args[:i]
    result_files = args[(i+1):]
else:
    list_files = [ args[0] ]
    result_files = args[1:]


## make list of matrices

matrices = [ ]
for list_file in list_files:
    ls = open(list_file, 'r')
    for line in ls.readlines():
        line = line.strip()
        if not os.path.isabs(line):
            line = os.path.join(os.getcwd(), line)
        matrices.append(line)


## process result files

rank_line_re = re.compile(r'^[^ ]+/(?P<program>[idqz]?rank(-int|-mod|-double)?([_-]mpi|[_-]omp){0,2}|lbrank_(?:bb|se_linear|se_none)) file:')
whitespace_re = re.compile(r'\s+', re.X)
garbage_re = re.compile(r'(MPI process \(rank: \d+\) terminated unexpectedly|\[[0-9a-z]+:\d+\] \[ *\d+\]).*')

def process_input_file(input_file_name):
    """
    Return a sequence of triples `(program, file, wctime)`
    extracted from the given output file.
    """
    input_file = open(input_file_name, 'r')
    logger.log(1, "Processing file %s ...", input_file_name)
    for line in input_file.readlines():
        line = line.strip()
        match = rank_line_re.match(line)
        if not match:
            continue
        program = match.group('program')
        line = garbage_re.sub('', line)
        words = line.split(' ')[1:]
        outcome = dict(w.split(':',1) for w in words)
        if not outcome.has_key('rank'):
            # computation failed, go to next line
            continue
        if outcome.has_key('wctime'):
            wctime = float(outcome['wctime'])
        elif outcome.has_key('cputime'):
            wctime = float(outcome['cputime'])
        elif outcome.has_key('time'):
            wctime = float(outcome['time'])
        else:
            logger.warning("Cannot read CPU/wall-clock time from input line '%s'"
                           " - ignoring." % line)
            continue
        yield (program, outcome['file'], wctime)
    input_file.close()

data = [ ]
for filename in result_files:
    data.extend(list(process_input_file(filename)))
        

## bin processed data according to filename/matrix

D = { }
for program, filename, wctime in data:
    if not D.has_key(program):
        D[program] = { }
    m = D[program]
    if m.has_key(filename):
        # keep largest running time
        if wctime > m[filename]:
            m[filename] = wctime
    else:
        m[filename] = wctime

## output results

for program in D.keys():
    for filename in matrices:
        if not m.has_key(filename):
            print ("run/gcc443-r66/%s file:%s cputime:%d" 
                   % (program, filename, int(options.duration)))
