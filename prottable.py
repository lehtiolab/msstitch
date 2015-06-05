#!/usr/bin/python3

"""
prottable -- Creating and modifying protein tables
"""

import argparse
import os
import app.drivers.prottable as drivers


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
                    'addprotdata - Add protein data (description, coverage,\n'
                    '# PSMs, etc.) to a table with protein accessions\n'
                    'Use with --dbfile\n\n'
                    'buildquant - Create protein quant data from a lookup\n'
                    'database. E.g. when multiple protein quant tables have\n'
                    'been read into the lookup and will be combined. Use \n'
                    'with --dbfile, --proteindata but NOT with -i.\n\n'
                    'addms1quant - Add MS1 quantification data from a\n'
                    'peptide table containing precursor quant areas. Needs\n'
                    '--peptable.',
                    required=True
                    )
parser.add_argument('-i', dest='infile',
                    help='Protein accession table if applicable\n',
                    type=lambda x: parser_file_exists(parser, x),
                    )
parser.add_argument('-d', dest='outdir', required=True,
                    help='Directory to output in',
                    type=lambda x: parser_file_exists(parser, x))
parser.add_argument('--dbfile', dest='lookup', help='Protein group '
                    'lookup database in SQLite format. Can be created using '
                    'mslookup.py command.',
                    type=lambda x: parser_file_exists(parser, x))
parser.add_argument('--peptable', dest='pepfile', help='Peptide table file '
                    'containing data for protein table, for example precursor '
                    'area amounts.',
                    type=lambda x: parser_file_exists(parser, x))
parser.add_argument('--proteindata', dest='proteindata', help='Include '
                    'protein group data such as coverage in output. Flag.',
                    action='store_const', default=False, const=True)

args = parser.parse_args()

commandmap = {
    'addprotdata': drivers.AddProteinInfoDriver,
    'buildquant': drivers.BuildProteinTableDriver,
    'addms1quant': drivers.AddPrecursorAreaDriver,
    #'createprottable': drivers.CreateProteinTableDriver,
}

command = commandmap[args.command](**vars(args))
command.run()
