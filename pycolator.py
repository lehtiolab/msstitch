#!/usr/bin/python
# FIXME docstring bad command format
"""
pycolator -- Percolator output manipulation

USAGE:
    python pycolator.py [option] [input file] [output file]
EXAMPLE:
    python pycolator.py -c splittd input.xml output.xml
"""

import argparse
import os
from app.drivers import pycolator as drivers


def parser_file_exists(currentparser, fn):
    if not os.path.exists(fn):
        currentparser.error('Input file %s not found!' % fn)
    else:
        return fn

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-c', dest='command', type=str, help='How to manipulate the input:\n'
'splittd        - Splits target and decoy data, multiple inputs to multiple outputs\n'
'merge          - Merges xml files. nothing else.\n'
'mergebest      - Merges xml files and only includes best scoring unique peptides\n'
'filterknown    - Filters out peptides that are found in a certain FASTA search\n'
'                 space which is passed using the -b flag. Then merges xml files and only\n'
'                 includes best scoring unique peptides'
'qvality        - Runs qvality on two inputfiles: one containing target and \n'
'                 containing decoy data. It is assumed that the first file \n'
'                 specified under -i is the target and the second is the\n'
'                 decoy.'
'reassign - Reassigns statistics from a qvality output file onto a single'
'           percolator input file. Needs -q flag.',
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
parser.add_argument('-f', dest='feattype', help='Feature type to use for '
                    'qvality. Can either be psm or peptide.')
parser.add_argument('-o', dest='options', nargs='+',
                    help='Extra options that may be passed to qvality.')


parser.add_argument('-q', dest='qvalityout', help='Qvality output file. '
                    'Required when using the reassign command.',
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
    'qvality'    : drivers.QvalityDriver,
    'reassign'   : drivers.ReassignmentDriver,
    }


command = commandmap[args.command](**vars(args))
command.run()

