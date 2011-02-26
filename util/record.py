#! /usr/bin/env python
#
"""
Import run data from irank* output.
"""
__docformat__ = 'reStructuredText'

import csv
import logging
from optparse import OptionParser
import os
import os.path
import re
import sqlite3
import subprocess
import sys

try:
    from collections import namedtuple
except ImportError:
    from namedtuple import namedtuple


cmdline = OptionParser()
cmdline.add_option('-a', '--csv', action='store', dest='csv_file', default=None,
                   metavar='FILE',
                   help="Add records into the specified CSV file.")
cmdline.add_option('-b', '--sql', action='store', dest='sql_file', default=None,
                   metavar='FILE',
                   help="Add records into the specified SQLite3 file.")
cmdline.add_option('-D', '--defaults', dest='defaults', action='append', metavar="SPEC",
                   help="Provide defaults for inserted fields, in the form"
                   " of a comma-separated list of KEY:VALUE pairs.")
cmdline.add_option('-v', '--verbose', dest='verbose', action='count', default=0,
                   help="Report more details on program execution.")
(options, args) = cmdline.parse_args()

if (options.csv_file is not None and options.sql_file is not None) \
        or (options.csv_file is None and options.sql_file is None):
    sys.stderr.write("Please specify one and only one of options '-a' and '-b'.")
    sys.exit(1)

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


class BaseStatsDb(object):
    """
    Abstract base class for tabular format data.

    Derived classes should implement the `load`, `save` and `close` methods.
    """
    def __init__(self, schema, defaults={}):
        self._data = set()
        self._schema = schema
        self._aux = { }
        self._fieldnames = [ ]
        for fieldname, aux in schema:
            self._aux[fieldname] = aux
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
        raise NotImplementedError("Abstract method `BaseStatsDb.load` called.")
    def save(self, filename=None):
        raise NotImplementedError("Abstract method `BaseStatsDb.save` called.")
    def close(self):
        raise NotImplementedError("Abstract method `BaseStatsDb.close` called.")

class CsvStatsDb(BaseStatsDb):
    """
    In-memory manipulation of a CSV file database.
    """
    def __init__(self, csvfile, schema, defaults={}):
        BaseStatsDb.__init__(self, schema, defaults)
        self._csvfile = csvfile
    def load(self, filename=None):
        if filename is None:
            filename = self._csvfile
        csvfile = open(filename, 'r')
        for row in csv.DictReader(csvfile):
            for f, conv in self._aux.items():
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

class SqliteStatsDb(BaseStatsDb):
    """
    Manipulate tabular data backed by a SQL storage
    """
    def __init__(self, database, schema, defaults={}, tablename='rundata'):
        BaseStatsDb.__init__(self, schema, defaults)
        self._conn = sqlite3.connect(database)
        self._conn.row_factory = sqlite3.Row
        self._sql = self._conn.cursor()
        self._tablename = tablename
    def load(self, tablename=None):
        if tablename is None:
            tablename = self._tablename
        self._sql.execute("CREATE TABLE IF NOT EXISTS %s (%s);"
                          % (tablename, 
                             str.join(", ", 
                                      [ ("%s %s" % (fieldname, fieldtype))
                                        for fieldname,fieldtype
                                        in self._aux.items() ])))
        self._sql.execute("SELECT %s FROM %s;"
                          % (str.join(',', self._fieldnames), tablename))
        for row in self._sql.fetchall():
            self.add(row)
    def save(self, database=None, tablename=None):
        if tablename is None:
            tablename = self._tablename
        for row in self._data:
            self._sql.execute("INSERT INTO %s (%s) VALUES (%s);"
                              % (tablename, 
                                 str.join(",", self._fieldnames),
                                 str.join(',', ['?'] * len(self._fieldnames))),
                              row)
    def close(self):
        self._conn.commit()
        self._conn.close()


if options.csv_file:
    def nullable(fn):
        def fn_or_None(x):
            if x is None or x == '':
                return None
            else:
                return fn(x)
        return fn_or_None
    db = CsvStatsDb(options.csv_file, (
            ('exitcode',   int),            # UNIX exit code
            ('matrix',     str),            # matrix file processed
            ('rows',       int),            # no. of matrix rows
            ('cols',       int),            # no. of matrix columns
            ('nonzero',    int),            # no. of nonzero entries
            ('rank',       nullable(int)),  # computed matrix rank
            ('cputime',    nullable(float)) # CPU time consumed on MPI rank 0
            ('wctime',     nullable(float)) # wall-clock time on MPI rank 0
            ('omp',        nullable(int)),  # no. of OpenMP threads (0 = no OpenMP)
            ('mpi',        nullable(int)),  # no. of MPI ranks (0 = no MPI)
            ('w',          int),            # width of a band of columns
            ('compiler',   nullable(str)),  # compiler used for compiling the C++ source
            ('mpilib',     nullable(str)),  # MPI library used 
            ('revno',      nullable(int)),  # BZR revision no. of the C++ source
            ('provenance', str),            # file name where this info was extracted from
            ), defaults)
    if os.path.exists(options.csv_file): 
        db.load(options.csv_file)
elif options.sql_file:
    db = SqliteStatsDb(options.sql_file, (
            ('exitcode',   'INTEGER'), # UNIX exit code
            ('matrix',        'TEXT'), # matrix file processed
            ('rows',       'INTEGER'), # no. of matrix rows
            ('cols',       'INTEGER'), # no. of matrix columns
            ('nonzero',    'INTEGER'), # no. of nonzero entries
            ('rank',       'INTEGER'), # computed matrix rank
            ('cputime',      'FLOAT'), # CPU time consumed on MPI rank 0
            ('wctime',       'FLOAT'), # wall-clock time on MPI rank 0
            ('omp',        'INTEGER'), # no. of OpenMP threads (0 = no OpenMP)
            ('mpi',        'INTEGER'), # no. of MPI ranks (0 = no MPI)
            ('w',          'INTEGER'), # width of a band of columns
            ('compiler',   'INTEGER'), # compiler used for compiling the C++ source
            ('mpilib',     'INTEGER'), # MPI library used 
            ('revno',      'INTEGER'), # BZR revision no. of the C++ source
            ('provenance',    'TEXT'), # file name where this info was extracted from
            ), defaults)
    if os.path.exists(options.sql_file): 
        db.load()


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

        
