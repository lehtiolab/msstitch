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

def parser_file_exists(parser, fn):
    if not os.path.exists(fn):
        parser.error('Input file %s not found!' % fn)
    else:
        return fn

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-c', dest='command', type=str, help='How to manipulate the input:\n'
'splittd - splits target and decoy data, multiple inputs to multiple outputs\n'
'mergebest - merges xml files and only includes best scoring unique peptides',
required=True
)

parser.add_argument('-i', dest='infile', nargs='+', 
        type=lambda x:parser_file_exists(parser, x), 
        required=True, help='Input file(s)')

parser.add_argument('-d', dest='outdir', required=True, help='Directory to output in')
parser.add_argument('-s', dest='score', help='Score to filter unique peptides '
'on (only for command mergebest)')
#parser.add_argument('-o', dest='outfiles', nargs='+', help='Output file(s)')

args = parser.parse_args()

commandmap = {
    'splittd'   : drivers.split_target_decoy,
    'mergebest' : drivers.merge_unique_best_scoring_peptides,
    }


commandmap[args.command](args.infile, args.outdir)

