#!/usr/bin/python3
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
                    'percotsv       - Add percolator data to a  TSV with \n'
                    'MSGF+ output. Specify TSV file with -i, mzid file with \n'
                    '--mzid.\n'
                    'mergetsv       - Merges multiple TSV tables of MSGF+ \n'
                    'output. Make sure headers are same in all files.\n'
                    'quanttsv       - Add quantitative data from openMS\n'
                    'consensusXML to a tab separated file with\n'
                    'PSMs, their statistics and their quantification data.\n'
                    'Needs except these also the corresponding mzML spectra\n'
                    'files to correlate retention time to scan nrs.\n',
                    required=True
                    )
parser.add_argument('-i', dest='infile', help='TSV table of mzIdentML',
                    type=lambda x: parser_file_exists(parser, x),
                    required=True)
parser.add_argument('--multifiles', dest='multifile_input', nargs='+',
                    type=lambda x: parser_file_exists(parser, x),
                    help='Multiple input files for use in e.g. merging data.'
                    ' PSMs will be picked from these and e.g. merged in '
                    'the file specified with -i.')
parser.add_argument('-d', dest='outdir', required=True,
                    help='Directory to output in',
                    type=lambda x: parser_file_exists(parser, x))
parser.add_argument('--mzid', dest='mzid', help='mzIdentML file',
                    type=lambda x: parser_file_exists(parser, x))
parser.add_argument('--quants', dest='quants', help='Quants from OpenMS in '
                    'consXML format', nargs='+',
                    type=lambda x: parser_file_exists(parser, x))
parser.add_argument('--spectra', dest='spectra', help='mzML files', nargs='+',
                    type=lambda x: parser_file_exists(parser, x))

# not supported yet
#parser.add_argument('--allpsms', dest='allpsms', action='store_true',
#                    help='All PSMs from a single scan should be included, '
#                    'not only the best scoring one.')

args = parser.parse_args()

commandmap = {
    'percotsv': drivers.MzidPercoTSVDriver,
    'mergetsv': drivers.MzidTSVConcatenateDriver,
    'quanttsv': drivers.TSVQuantDriver,
}

command = commandmap[args.command](**vars(args))
command.run()
