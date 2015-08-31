#!/usr/bin/python3

"""
prottable -- Creating and modifying protein tables
"""

import argparse
import os
from app.drivers.prottable import (probability, info, merge, precursorarea,
                                   create_empty, create_labelfree, qvality,
                                   fdr)


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
                    'with --dbfile, --proteindata, --precursor, --isobaric\n'
                    'but NOT with -i.\n\n'
                    'addms1quant - Add MS1 quantification data from a\n'
                    'PSM table containing precursor quant areas. Needs\n'
                    'a psmtable specified with -i.\n\n'
                    'createlabelfree - Create protein table from PSM table\n'
                    'containing precursor area quant information. Needs\n'
                    'input from --psmtable.\n\n'
                    'emptyprottable - Create protein table from PSM table\n'
                    'containing no quant data, resulting in one column with\n'
                    'master proteins only.\n\n'
                    'addprob - Add protein probabilities from a\n'
                    'peptide table posterior error probabilities. Needs\n'
                    '--peptable, and probabilities are calculated \n'
                    'as in Nesvizhskii et al. (2003) Anal.Chem., eq 3.\n\n'
                    'protqvality - Run qvality on protein tables containing \n'
                    'target (-i) proteins and decoy (--decoy) proteins.\n\n'
                    'addfdr - Add protein FDR to protein table by comparing\n'
                    'protein probability with qvality lookup table. Needs \n'
                    'to have qvality output file specified with --qvality',
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
parser.add_argument('--psmtable', dest='psmfile', help='PSM table file '
                    'containing data for protein table, for example precursor '
                    'area amounts.',
                    type=lambda x: parser_file_exists(parser, x))
parser.add_argument('--peptable', dest='pepfile', help='Peptide table file '
                    'containing data for protein table, for example '
                    'peptide probabilities.',
                    type=lambda x: parser_file_exists(parser, x))
parser.add_argument('--decoy', dest='decoyfn', help='Protein table containing '
                    'decoy proteins for running qvality',
                    type=lambda x: parser_file_exists(parser, x))
parser.add_argument('--feattype', dest='feattype', help='Score type to use for '
                    'qvality, reassign features or pout2tsv. Can either be '
                    'psm or peptide.')
parser.add_argument('-o', dest='options', nargs='+',
                    help='Extra options that may be passed to qvality.'
                    'Option form: -o ***flag value ***flag ***flag value')
parser.add_argument('--qvality', dest='qvalityfile', help='Qvality output '
                    'table where q value and PEP can be looked up in.',
                    type=lambda x: parser_file_exists(parser, x))
parser.add_argument('--setname', dest='setname', help='Name of biological '
                    'set which to use when adding protein info to table. '
                    'Must be used with addprotdata.',
                    default=False, type=str)
parser.add_argument('--proteindata', dest='proteindata', help='Include '
                    'protein group data such as coverage in output. Flag.',
                    action='store_const', default=False, const=True)
parser.add_argument('--precursor', dest='precursor', help='Build protein '
                    'table which contains precursor quant data. Flag.',
                    action='store_const', default=False, const=True)
parser.add_argument('--isobaric', dest='isobaric', help='Build protein '
                    'table which contains isobaric quant data. Flag.',
                    action='store_const', default=False, const=True)
parser.add_argument('--probability', dest='probability', help='Build protein '
                    'table which contains probability data. Flag.',
                    action='store_const', default=False, const=True)

args = parser.parse_args()

commandmap = {
    'addprotdata': info.AddProteinInfoDriver,
    'buildquant': merge.BuildProteinTableDriver,
    'addms1quant': precursorarea.AddPrecursorAreaDriver,
    'addprob': probability.AddProteinProbability,
    'createlabelfree': create_labelfree.CreateLabelfreeProteinDriver,
    'emptyprottable': create_empty.CreateEmptyDriver,
    'protqvality': qvality.ProttableQvalityDriver,
    'protfdr': fdr.ProttableFDRDriver,
}

command = commandmap[args.command](**vars(args))
command.run()
