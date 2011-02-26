#! /usr/bin/env python
#
"""
Analyze run data from rank*/lbrank* output and output a set of files
such that each bacth of matrices will be processed in a given amount
of time.
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
cmdline.add_option('-q', '--duration', dest='durations', action='append', metavar='FIELD', default=[],
                   help="Allow a batch to have the specified duration."
                   " Repeat multiple times for mutliple different durations.")
cmdline.add_option('-o', '--output-dir', dest='output_dir', action='store', metavar="SPEC", default='/tmp',
                   help="Create batch files in the specified directory.")
cmdline.add_option('-s', '--safety-factor', dest='factor', action='store', metavar="NUM", default='1.2',
                   help="Multiply the CPU time by this factor to get a safe estimation of wall-clock time.")
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


try:
    factor = float(options.factor)
    if factor < 1:
        raise TypeError()
except TypeError:
    logger.error("Argument to option -s/--safety-factor must be a number greater than 1.")
    sys.exit(1)

durations_ = [ ]
for duration in options.durations:
    if ',' in duration:
        durations_.extend(duration.split(','))
    else:
        durations_.append(duration)
def str2seconds(timespec):
    if ':' in timespec:
        parts = timespec.split(':')
        if len(parts) == 2:
            h,m = parts
            s = 0
        elif len(parts) == 3:
            h,m,s = parts
        return h*3600 + m*60 + s
    elif timespec.endswith('m'):
        return 60*int(timespec[:-1])
    elif timespec.endswith('h'):
        return 3600*int(timespec[:-1])
    elif timespec.endswith('d'):
        return 24*3600*int(timespec[:-1])
durations = [ str2seconds(d) for d in durations_ ]
logger.log(1, "Will group files in batches of sizes %s" 
           % str.join(", ", [str(dur) for dur in durations]))

## process input files
rank_line_re = re.compile(r'^[^ ]+/(?P<program>[idqz]?rank(-int|-mod|-double)?([_-]mpi|[_-]omp){0,2}|lbrank_(?:bb|se_linear|se_none)) file:[^ ]+ ([a-z]+:[0-9]+ *)+ *((cpu|wc)?time:[0-9]+\.[0-9]+ *)+ *$')
whitespace_re = re.compile(r'\s+', re.X)
garbage_re = re.compile(r'(MPI process \(rank: \d+\) terminated unexpectedly|\[[0-9a-z]+:\d+\] \[ *\d+\]|.*SIG[A-Z]+[0-9]*|/[/a-z0-9_]+: line [0-9]+: [0-9]+ Aborted|\[[a-z0-9]+:[0-9]+\] \*\*\* Process received signal).*')

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
        if not (outcome.has_key('time') 
                or outcome.has_key('cputime') 
                or outcome.has_key('wctime')):
            # computation failed, go to next line
            continue
        if outcome.has_key('wctime'):
            wctime = 1 + int(float(outcome['wctime']))
        elif outcome.has_key('cputime'):
            wctime = 1 + int(factor*float(outcome['cputime']))
        elif outcome.has_key('time'):
            wctime = 1 + int(factor*float(outcome['time']))
        else:
            logger.warning("Cannot read CPU/wall-clock time from input line '%s'"
                           " - ignoring." % line)
            continue
        yield (program, outcome['file'], wctime)
    input_file.close()

data = [ ]
for filename in args:
    data.extend(list(process_input_file(filename)))


## group data according to program, then keep only max wctime for
## each filename
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


## bin processed files according to duration

class Bag(set):
    def __new__(cls, duration):
        self = set.__new__(cls)
        self.max_duration = duration
        self.current_duration = 0
        return self
    def __init__(self, duration):
        pass
    def add(self, item):
        if item[1] + self.current_duration < self.max_duration:
            self.current_duration += item[1]
            set.add(self, item)
        else:
            raise ValueError("Adding item '%s' would exceed maximum duration %ds" 
                             % (item, self.max_duration))

# Browse items in decreasing duration order; bags are sorted in
# increasing capacity order; this should allow us to just add an item
# to the first bag that accepts it.  Definitely not efficient, but it
# should work.
batches = { }
for program, data in D.items():
    if not batches.has_key(program):
        batches[program] = list()
    bags = batches[program]
    def cmp_by_wctime(a, b):
        return cmp(a[1], b[1])
    for item in sorted(data.items(), cmp=cmp_by_wctime, reverse=True):
        logger.log(3, "Processing file '%s', duration '%d'" % (item[0], item[1]))
        added = False
        for bag in bags:
            try: 
                bag.add(item)
                added = True
                break
            except ValueError:
                continue
        if not added:
            try:
                q = min([ d for d in durations if d > item[1] ])
            except ValueError:
                logger.log(0, "Computation of matrix '%s' by program '%s' does not fit in any specified duration"
                           " (needs: %ds).  Ignoring." % (os.path.basename(item[0]), program, item[1]))
                continue
            logger.log(1, "Adding one more bag of size %ds" % q)
            bag = Bag(q)
            bag.add(item)
            inserted = False
            for n in range(len(bags)):
                if q < bags[n].max_duration:
                    bags.insert(n, bag)
                    inserted = True
                    break
            if not inserted:
                bags.append(bag)


## create output files

if not os.path.exists(options.output_dir):
    os.mkdir(options.output_dir)

n = 1 + int(math.log10(len(bags)))
fmt = os.path.join(options.output_dir, 
                   str.join('', ['%s.b%0' , str(n), 'd', ':s_rt=%d']))

for program, batch in batches.items():
    logger.log(1, "Outputting %d batches for program '%s' ..." % (len(batch), program))
    for n, bag in enumerate(batch):
        if len(bag) == 0:
            continue
        # the bag may be half-full, so recompute the needed duration here
        q = min([ dur for dur in durations if dur > bag.current_duration ])
        name = fmt % (program, n, int(1.05 * q)-1)
        logger.log(1, "Creating output file '%s' with %d items (duration: %ds)" 
                   % (name, len(bag), bag.current_duration))
        output = open(name, 'w+')
        for item in bag:
            output.write("%s\n" % item[0])
        output.close()
