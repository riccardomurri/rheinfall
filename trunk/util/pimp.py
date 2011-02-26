#! /usr/bin/env python
#
"""
Import run data from rank.cpp output.
"""
__docformat__ = 'reStructuredText'

import csv
import logging
from optparse import OptionParser
import os
import os.path
import re
import subprocess
import sys

try:
    from collections import namedtuple
except ImportError:
    from namedtuple import namedtuple


cmdline = OptionParser()
cmdline.add_option('-f', '--file', action='store', dest='db_file', default='stats.csv',
                   metavar='FILE',
                   help="Add records into the specified CSV file.")
cmdline.add_option('-D', '--defaults', dest='defaults', action='append', metavar="SPEC",
                   help="Provide defaults for inserted fields, in the form"
                   " of a comma-separated list of KEY:VALUE pairs.")
cmdline.add_option('-v', '--verbose', dest='verbose', action='count', default=0,
                   help="Report more details on program execution.")
(options, args) = cmdline.parse_args()


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
    }
if options.defaults is not None:
    for kv in str.join(',', options.defaults).split(','):
        k, v = kv.split(':')
        defaults[k] = v


class StatsDb(object):
    """
    In-memory manipulation of a CSV file database.
    """
    def __init__(self, csvfile, schema, defaults={}):
        self._data = set()
        self._schema = schema
        self._csvfile = csvfile
        self._converters = { }
        self._fieldnames = [ ]
        for fieldname, converter in schema:
            self._converters[fieldname] = converter
            self._fieldnames.append(fieldname)
        self._recordtype = namedtuple('record_%s' % self.__class__.__name__, 
                                      self._fieldnames)
        self._defaults = defaults
    def add(self, row):
        for k,v in self._defaults.items():
            if not row.has_key(k):
                row[k] = v
        self._data.add(self._recordtype(**row))
    def load(self, filename=None):
        if filename is None:
            filename = self._csvfile
        csvfile = open(filename, 'r')
        for row in csv.DictReader(csvfile):
            for f, conv in self._converters.items():
                row[f] = conv(row[f])
            self.add(row)
        csvfile.close()
    def save(self, filename=None):
        if filename is None:
            filename = self._csvfile
        if os.path.exists(filename):
            csvfile = open(filename, 'a')
            write_header = False
        else:
            csvfile = open(filename, 'w')
            write_header = True
        writer = csv.writer(csvfile)
        if write_header:
            writer.writerow(self._fieldnames)
        writer.writerows(self._data)
        csvfile.close()
    def close(self):
        pass

def nullable(fn):
    def fn_or_None(x):
        if x is None or x == '':
            return None
        else:
            return fn(x)
    return fn_or_None

db = StatsDb(options.db_file, (
        ('exitcode', int), # UNIX exit code
        ('matrix', str),   # matrix file processed
        ('rows', int),     # no. of matrix rows
        ('cols', int),     # no. of matrix columns
        ('nonzero', int),  # no. of nonzero entries
        ('rank', nullable(int)),     # computed matrix rank
        ('cputime', nullable(float)),# CPU time consumed on MPI rank 0
        ('wctime', nullable(float)), # wall-clock time on MPI rank 0
        ('omp', nullable(int)),      # no. of OpenMP threads (0 = no OpenMP)
        ('mpi', nullable(int)),      # no. of MPI ranks (0 = no MPI)
        ('compiler', nullable(str)), # compiler used for compiling the C++ source
        ('mpilib', nullable(str)),   # MPI library used 
        ('revno', nullable(int)),    # BZR revision no. of the C++ source
        ('provenance', str),         # file name where this info was extracted from
        ), defaults)
if os.path.exists(options.db_file):
    db.load(options.db_file)


## main 

rank_line_re = re.compile(r'^[^ ]+/rank(_mpi|_omp){0,2} file:')
jobid_re = re.compile(r'.[eo]([0-9]+)')
whitespace_re = re.compile(r'\s+', re.X)
garbage_re = re.compile(r'(MPI process \(rank: \d+\) terminated unexpectedly|\[[0-9a-z]+:\d+\] \[ *\d+\]).*')

for input_file_name in args:
    input_file = open(input_file_name, 'r')
    logger.log(1, "Processing file %s ...", input_file_name)
    for line in input_file.readlines():
        line = line.strip()
        if not rank_line_re.match(line):
            continue
        line = garbage_re.sub('', line)
        words = line.split(' ')[1:]
        outcome = { 'provenance':input_file_name }
        outcome.update(dict([(k,v) for k,v in defaults.items()]))
        outcome.update(dict(w.split(':',1) for w in words))
        # extract matrix name from file name
        outcome['matrix'] = os.path.splitext(os.path.basename(outcome['file']))[0]
        del outcome['file']
        if not outcome.has_key('rank'):
            # computation failed, try to get exitcode from SGE
            match = jobid_re.search(input_file_name)
            jobid = match.group(1)
            logger.log(2, "Running 'qacct -j %s' ...", jobid)
            qacct = subprocess.Popen(["qacct", "-j", jobid], 
                                     stdout=subprocess.PIPE).communicate()[0]
            for kv in qacct.split('\n')[1:]:
                kv = kv.strip()
                k, v = whitespace_re.split(kv, 1)
                if k == 'exit_status':
                    outcome['exitcode'] = int(v)
                    break
            # supply NULL values for unknown fields
            outcome['rank'] = None
            outcome['cputime'] = None
            outcome['wctime'] = None
        else:
            outcome['exitcode'] = 0
        db.add(outcome)
db.save()
db.close()

        
