#!/usr/bin/python
"""
ompstitch -- Combining OpenMS, mzIdent, percolator data

USAGE:
    ompstitch.py [option] [input files] [output file]
EXAMPLE:
    ompstitch.py -o quanttsv --spectra spectra.mzML --ids psms.mzid --quants quant.consXML psmquant.tsv
"""

import argparse
import os
from drivers import ompstitch as drivers


def parser_file_exists(currentparser, fn):
    if not os.path.exists(fn):
        currentparser.error('Input file %s not found!' % fn)
    else:
        return fn

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-c', dest='command', type=str,
                    help='How to manipulate the input:\n'
                    'quanttsv       -   Combine spectra (in mzML format), '
                    'identifications (in mzIdentML format) and quantification'
                    ' (consensusXML) to output a tab separated file with '
                    'PSMs, their statistics and their quantification data.\n'
                    'CURRENTLY ONLY WORKING WITH msgfplus-percolator mzidML',
                    required=True
                    )
parser.add_argument('-d', dest='outdir', required=True,
                    help='Directory to output in')
parser.add_argument('--spectra', dest='spectra',
                    type=lambda x: parser_file_exists(parser, x),
                    required=True, help='mzML spectra file')
parser.add_argument('--ids', dest='ids', help='mzIdentML file',
                    type=lambda x: parser_file_exists(parser, x))
parser.add_argument('--quants', dest='quants',
                    help='Quantifications in consensusXML format.',
                    type=lambda x: parser_file_exists(parser, x))

args = parser.parse_args()

commandmap = {
    'quanttsv': drivers.QuantTSVDriver,
}

command = commandmap[args.command](**vars(args))
command.run()
