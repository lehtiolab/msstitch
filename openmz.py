#!/usr/bin/python
# FIXME docstring wrong
"""
openmz -- Connecting mzIdentML with openMS

USAGE:
   openmz.py [option] [input files]
EXAMPLE:
   openmz.py -c quanttsv --mzidtsv psms.tsv --spectra file1.mzML file2.mzML -d /home/data/output
"""

import argparse
import os
from app.drivers.openmz import openmz as drivers


def parser_file_exists(currentparser, fn):
    if not os.path.exists(fn):
        currentparser.error('Input file %s not found!' % fn)
    else:
        return fn

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-c', dest='command', type=str,
                    help='How to manipulate the input:\n'
                    'quanttsv       - Add quantitative data from openMS '
                    'consensusXML to a tab separated file with '
                    'PSMs, their statistics and their quantification data. '
                    'Needs except these also the corresponding mzML spectra '
                    'files to correlate retention time to scan nrs.',
                    required=True
                    )
parser.add_argument('-i', dest='infile',
                    type=lambda x: parser_file_exists(parser, x),
                    required=True, help='Mzid TSV file')
parser.add_argument('-d', dest='outdir', required=True,
                    help='Directory to output in',
                    type=lambda x: parser_file_exists(parser, x))
parser.add_argument('--quants', dest='quants', help='Quants from OpenMS in '
                    'consXML format', nargs='+',
                    type=lambda x: parser_file_exists(parser, x))
parser.add_argument('--spectra', dest='spectra', help='mzML files', nargs='+',
                    type=lambda x: parser_file_exists(parser, x))


# TODO to be supported
#parser.add_argument('--allpsms', dest='allpsms', action='store_true',
#                    help='All PSMs from a single scan should be included, '
#                    'not only the best scoring one.')

args = parser.parse_args()

commandmap = {
    'quanttsv': drivers.TSVQuantDriver,
}

command = commandmap[args.command](**vars(args))
command.run()
