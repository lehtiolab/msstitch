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
                    'Use with --fasta and --protgroupdb\n',
                    required=True
                    )
parser.add_argument('-i', dest='infile',
                    help='Protein accession table, or if using \n'
                    'createprottable an MzidTSV file with protein groups',
                    type=lambda x: parser_file_exists(parser, x),
                    required=True)
parser.add_argument('-d', dest='outdir', required=True,
                    help='Directory to output in',
                    type=lambda x: parser_file_exists(parser, x))
parser.add_argument('--workdir', dest='workdir',
                    help='Working directory to output temporary files in. '
                    'Temporary files will be removed after finishing.'
                    )
parser.add_argument('--fasta', dest='fasta', help='FASTA sequence database, '
                    'to optionally use with proteingrouping to enable sorting '
                    'on coverage, and in case of UNIPROT FASTA, evidence '
                    'levels.', type=lambda x: parser_file_exists(parser, x))
parser.add_argument('--protgroupdb', dest='protgroupdb', help='Protein group '
                    'lookup database in SQLite format. Can be created using '
                    'mzidplus.py command.',
                    type=lambda x: parser_file_exists(parser, x))

args = parser.parse_args()

commandmap = {
    'addprotdata': drivers.AddProteinInfoDriver,
    #'createprottable': drivers.CreateProteinTableDriver,
}

command = commandmap[args.command](**vars(args))
command.run()
