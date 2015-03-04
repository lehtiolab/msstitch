#!/usr/bin/python3

"""
mslookup - Creating SQLite lookups for internal and external use
"""

import argparse
import os
import app.drivers.lookup as drivers


def parser_file_exists(currentparser, fn):
    if not os.path.exists(fn):
        currentparser.error('Input file %s not found!' % fn)
    else:
        return fn


def parser_value_in_list(currentparser, value, valuelist):
    if not value in valuelist:
        currentparser.error('Value {0} should be one of {1}'.format(
            value,
            ', '.join(valuelist)))
    else:
        return value


parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-c', dest='command', type=str,
                    help='How to manipulate the input:\n'
                    'spectra - Create SQLite lookup of spectra in mzML \n'
                    'format. Requires passing an mzML file to -i, but\n'
                    'neither --spectra nor --dbfile.\n'
                    'isoquant - Create lookup of isobaric quant data in\n'
                    'OpenMS consensusXML format. Use requires --spectra,\n'
                    '--dbfile with an sqlite lookup of spectra, and passing\n'
                    'a consensusXML file to -i.\n'
                    'ms1quant - Creates lookup of precursor quant data in \n'
                    'OpenMS featureXML format. Use requires --spectra,\n'
                    '--dbfile wiht an sqlite lookup of spectra, and passing '
                    'a featureXML file to -i',
                    required=True
                    )
parser.add_argument('-i', dest='infile', nargs='+',
                    help='The input files to create the lookup from. Can be\n'
                    'of type consensusXML, featureXML, percolator out XML, \n'
                    'tab separated PSM table. If order is important then it\n'
                    'is taken from the input order at the command line.',
                    type=lambda x: parser_file_exists(parser, x),
                    required=True)
parser.add_argument('-d', dest='outdir', required=True,
                    help='Directory to output in',
                    type=lambda x: parser_file_exists(parser, x))
parser.add_argument('--dbfile', dest='lookup', help='File containing an\n'
                    'SQLite database to build on. Some commands need this\n'
                    'in order to build on.',
                    type=lambda x: parser_file_exists(parser, x))
parser.add_argument('--spectra', dest='spectra', help='Spectra files in mzML\n'
                    'format. Multiple files can be specified, if order is\n'
                    'important, e.g. when matching them with quant data, the\n'
                    'order will be their input order at the command line.',
                    type=lambda x: parser_file_exists(parser, x), nargs='+')

args = parser.parse_args()

commandmap = {
    'spectra': drivers.SpectraLookupDriver,
    'isoquant': drivers.IsobaricQuantLookupDriver,
    'ms1quant': drivers.PrecursorQuantLookupDriver,
}

command = commandmap[args.command](**vars(args))
command.run()
