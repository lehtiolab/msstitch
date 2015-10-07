#!/usr/bin/python3

"""
prottable -- Creating and modifying protein tables
"""

import argparse
import os
from app.drivers.prottable import (probability, info, merge, isoquant,
                                   precursorarea, create_empty, bestpeptide,
                                   qvality, fdr)


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

                    'addisoquant - Add isobaric quantification data from a\n'
                    'proteintable containing this. Needs a second table\n'
                    'specified with --quantfile and --isobquantcolpattern\n\n'

                    'addms1quant - Add MS1 quantification data from a\n'
                    'PSM table containing precursor quant areas. Needs\n'
                    'a psmtable specified with --psmtable.\n\n'
                    'emptyprottable - Create protein table from PSM table\n'
                    'containing no quant data, resulting in one column with\n'
                    'master proteins only.\n\n'

                    'addprob - Add protein probabilities from a\n'
                    'peptide table posterior error probabilities. Needs\n'
                    '--peptable, and probabilities are calculated \n'
                    'as in Nesvizhskii et al. (2003) Anal.Chem., eq 3.\n\n'

                    'bestpeptide - Given the protein table and corresponding\n'
                    'peptide table, fetch the best scoring peptide for each\n'
                    'protein and annotates that score in the protein table.\n'
                    'Use with --scorecolpattern, --peptable, --logscore.\n\n'

                    'protqvality - Run qvality on protein (or tsv) tables\n'
                    'containing target (-i) proteins and decoy (--decoy)\n'
                    'proteins. Use with --feattype to use either protein\n'
                    'error probability (Nesvizhskii 2003) or Q score from\n'
                    'Savitski 2014 MCP.\n\n'

                    'pickqvality - Run qvality on protein tables\n'
                    'containing target (-i) proteins and decoy (--decoy)\n'
                    'proteins. Targets and decoys will be filtered according\n'
                    'Savitski et al. 2014 MCP, where the best scoring of a\n'
                    'pair of target/decoy proteins is retained and the other\n'
                    'protein discarded. Matching (reversed, tryptic reversed,'
                    ' scrambled) target and decoy FASTA files are needed to\n'
                    'determine the pairs, use --targetfasta, --decoyfasta \n\n'

                    'addfdr - Add protein FDR to protein table by comparing\n'
                    'score (peptide q-value, protein probability, etc)\n'
                    'with qvality lookup table. Needs \n'
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
parser.add_argument('--dbfile', dest='lookup', help='File containing a'
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
parser.add_argument('--quantfile', dest='quantfile', help='Protein table file '
                    'containing isobaric quant data to add to protein table.',
                    type=lambda x: parser_file_exists(parser, x))
parser.add_argument('--isobquantcolpattern', dest='isobquantcolpattern',
                    help='Unique text pattern to identify isobaric quant \n'
                    'columns in protein table.',
                    type=str, required=False)
parser.add_argument('--scorecolpattern', dest='scorecolpattern',
                    help='Regular expression pattern to find column where\n'
                    'score to filter on is written.',
                    type=str, required=False)
parser.add_argument('--logscore', dest='logscore',
                    help='Flag. When using proteinbest, e.g. q-values will\n'
                    'be converted to -log10 values.',
                    action='store_const', const=True, default=False)
parser.add_argument('--decoy', dest='decoyfn', help='Protein table containing '
                    'decoy proteins for running qvality',
                    type=lambda x: parser_file_exists(parser, x))
parser.add_argument('--targetfasta', dest='targetfasta',
                    help='FASTA file with target proteins to determine best\n'
                    'scoring proteins of target/decoy pairs for pickqvality.')
parser.add_argument('--decoyfasta', dest='decoyfasta',
                    help='FASTA file with decoy proteins to determine best\n'
                    'scoring proteins of target/decoy pairs for pickqvality.')
parser.add_argument('--feattype', dest='feattype', help='Score type to use for'
                    ' qvality. Can either be probability or qvalue.')
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
parser.add_argument('--fdr', dest='fdr', help='Build protein '
                    'table which contains protein FDR data. Flag.',
                    action='store_const', default=False, const=True)
parser.add_argument('--pep', dest='pep', help='Build protein '
                    'table which contains protein PEP data. Flag.',
                    action='store_const', default=False, const=True)

args = parser.parse_args()

commandmap = {
    'addprotdata': info.AddProteinInfoDriver,
    'buildquant': merge.BuildProteinTableDriver,
    'addisoquant': isoquant.AddIsobaricQuantDriver,
    'addms1quant': precursorarea.AddPrecursorAreaDriver,
    'addprob': probability.AddProteinProbability,
    'emptyprottable': create_empty.CreateEmptyDriver,
    'bestpeptide': bestpeptide.BestPeptidePerProtein,
    'protqvality': qvality.ProttableQvalityDriver,
    'pickqvality': qvality.PickedQvalityDriver,
    'protfdr': fdr.ProttableFDRDriver,
}

command = commandmap[args.command](**vars(args))
command.run()
