#!/usr/bin/python

"""
pycolator -- Percolator output manipulation

USAGE:
    python pycolator.py [option] [input file] [output file]
EXAMPLE:
    python pycolator.py splittd input.xml output.xml
"""

import argparse
import os
import drivers


def parser_file_exists(currentparser, fn):
    if not os.path.exists(fn):
        currentparser.error('Input file %s not found!' % fn)
    else:
        return fn

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-c', dest='command', type=str, help='How to manipulate the input:\n'
'splittd - splits target and decoy data, multiple inputs to multiple outputs\n'
'merge - merges xml files. nothing else.\n'
'mergebest - merges xml files and only includes best scoring unique peptides\n'
'filterknown - Filters out peptides that are found in a certain FASTA search\n'
'space which is passed using the -b flag. Then merges xml files and only\n'
'includes best scoring unique peptides',
required=True
)

parser.add_argument('-i', dest='infile', nargs='+', 
        type=lambda x:parser_file_exists(parser, x), 
        required=True, help='Input file(s)')

parser.add_argument('-d', dest='outdir', required=True, help='Directory to output in')
parser.add_argument('-s', dest='score', help='Score to filter unique peptides '
'on (only for command mergebest and filterknownmerge)', default='svm')
parser.add_argument('-b', dest='database', help='Database file(s). Make sure'
' they are included when filterknown command is used, since they will be'
' used to exclude peptides from.', nargs='+', 
                     type=lambda x:parser_file_exists(parser, x))

# FIXME make db files required after we figure out if supplying raw db files is
# ok performance wise. If too slow, we may switch to sqlite db.

#parser.add_argument('-o', dest='outfiles', nargs='+', help='Output file(s)')

args = parser.parse_args()

commandmap = {
    'splittd'    : drivers.SplitDriver,
    'merge'      : drivers.MergeDriver,
    'mergebest'  : drivers.MergeUniquePeptides,
    'filterknown': drivers.MergeUniqueAndFilterKnownPeptides,
    }


command = commandmap[args.command](**vars(args))
command.run()

