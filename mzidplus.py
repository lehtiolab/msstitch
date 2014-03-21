#!/usr/bin/python
# FIXME docstring wrong
"""
mzidplus -- Modifying MSGF+ mzIdentML output

USAGE:
   mzidplus.py [option] [input files]
EXAMPLE:
   mzidplus.py -c percotsv -i psms.mzid -d /data -o psmquant.tsv
"""

import argparse
import os
from app.drivers import mzidplus as drivers


def parser_file_exists(currentparser, fn):
    if not os.path.exists(fn):
        currentparser.error('Input file %s not found!' % fn)
    else:
        return fn

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-c', dest='command', type=str,
                    help='How to manipulate the input:\n'
                    'percotsv       - Convert MSGF+ output (in mzIdentML '
                    'format) to a tab separated file with '
                    'PSMs, their statistics and their quantification data. '
                    'This includes Percolator-generated data from pout2mzid.',
                    required=True
                    )
parser.add_argument('-d', dest='outdir', required=True,
                    help='Directory to output in')
parser.add_argument('--ids', dest='ids', help='mzIdentML file',
                    type=lambda x: parser_file_exists(parser, x))

# not supported yet
#parser.add_argument('--allpsms', dest='allpsms', action='store_true',
#                    help='All PSMs from a single scan should be included, '
#                    'not only the best scoring one.')

args = parser.parse_args()

commandmap = {
    'percotsv': drivers.MzidPercoTSVDriver,
}

command = commandmap[args.command](**vars(args))
command.run()
