#!/usr/bin/python3
# FIXME docstring bad command format
"""
pycolator -- Percolator output manipulation

USAGE:
    python pycolator.py [option] [input file] [output file]
EXAMPLE:
    python pycolator.py -c splittd input.xml output.xml
"""

import argparse
import os
from app.drivers.pycolator import splitmerge
from app.drivers.pycolator import filters
from app.drivers.pycolator import lookup
from app.drivers.pycolator import stats
from app.drivers.pycolator import converters


def parser_file_exists(currentparser, fn):
    if not os.path.exists(fn):
        currentparser.error('Input file %s not found!' % fn)
    else:
        return fn

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-c', dest='command', type=str, help='How to manipulate the input:\n'
'splittd        - Splits target and decoy data, multiple inputs to multiple outputs\n'
'merge          - Merges percolator xml files. nothing else. Use -i for\n'
'                 base file, and specify files with --multifiles\n'
'trypticlookup  - Creates a lookup DB from a FASTA file for use with e.g. filterknown\n'
'                 --cutproline and --ntermwildcards can be used.\n'
'filteruni      - Only includes best scoring unique peptides in a (merged) file\n'
'filterknown    - Filters out peptides that are found in a certain lookup DB\n'
'                 passed using the -b flag. Use with or without\n'
'                 --cutproline and --ntermwildcards\n'
'filterlen      - Filters out peptides that exceed --maxlen and --minlen\n'
'qvality        - Runs qvality on an inputfile: target and decoy data.\n'
'                 When using separate files for target and decoy, \n'
'                 use --decoy to specify the decoy input file, and -f\n'
'                 to specify feature type (psm or peptide)\n'
'reassign       - Reassigns statistics from a qvality output file onto a single'
'                 percolator input file. Needs -q flag.\n'
'pout2tsv       - Converts a percolator output file to tab separated format\n'
'                 Use -f to specify feature type (psm or peptide), otherwise\n'
'                 both psm and peptide tsv files are written.',
required=True
)

parser.add_argument('-i', dest='infile',
        type=lambda x: parser_file_exists(parser, x),
        required=True, help='Input file')
parser.add_argument('--multifiles', dest='multifile_input', nargs='+',
                    type=lambda x: parser_file_exists(parser, x),
                    help='Multiple input files for use in e.g. merging data.'
                    ' Features will be picked from these and e.g. merged in '
                    'the file specified with -i.')
parser.add_argument('-d', dest='outdir', required=True, help='Directory to output in')
parser.add_argument('-s', dest='score', help='Score to filter unique peptides '
'on (only for command mergebest and filterknownmerge)', default='svm')
parser.add_argument('--maxlen', dest='maxlength', help='Maximum length of '
                    'peptide to be included in filtered data.', default=None)
parser.add_argument('--minlen', dest='minlength', help='Minimum length of '
                    'peptide to be included in filtered data.', default=None)

parser.add_argument('-b', dest='database', help='Database file. Make sure'
' they are included when filterknown command is used, since they will be'
' used to exclude peptides from.',
                     type=lambda x:parser_file_exists(parser, x))
parser.add_argument('--cutproline', dest='proline', help='With flag, trypsin is '
                    'considered to cut before a proline residue. The filter '
                    'known will filter against both cut and non-cut peptides.',
                    action='store_const', const=True, default=False)
parser.add_argument('--ntermwildcards', dest='falloff', help='With flag, the filter '
                    'known will filter against both intact peptides and those '
                    'that match to the C-terminal part of a tryptic peptide '
                    'from the database. Database should be built with this'
                    'flag in order for the lookup to work, since sequences'
                    'will be stored and looked up reversed',
                    action='store_const', const=True, default=False)

parser.add_argument('--decoy', dest='decoyfn',
                    type=lambda x: parser_file_exists(parser, x),
                    help='Decoy input file (percolator out XML) for qvality')
parser.add_argument('--target', dest='targetfn',
        type=lambda x: parser_file_exists(parser, x),
        help='Target input file for qvality')
parser.add_argument('-f', dest='feattype', help='Feature type to use for '
                    'qvality or pout2tsv. Can either be psm or peptide.')
parser.add_argument('-o', dest='options', nargs='+',
                    help='Extra options that may be passed to qvality.'
                    'Option form: -o ***flag value ***flag ***flag value')
parser.add_argument('-q', dest='qvalityout', help='Qvality output file. '
                    'Required when using the reassign command.',
                    type=lambda x:parser_file_exists(parser, x))

# FIXME make db files required after we figure out if supplying raw db files is
# ok performance wise. If too slow, we may switch to sqlite db.

#parser.add_argument('-o', dest='outfiles', nargs='+', help='Output file(s)')

args = parser.parse_args()

commandmap = {
    'splittd'    : splitmerge.SplitDriver,
    'merge'      : splitmerge.MergeDriver,
    'trypticlookup': lookup.CreateLookup,
    'filteruni'  : filters.FilterUniquePeptides,
    'filterlen'  : filters.FilterPeptideLength,
    'filterknown': filters.FilterKnownPeptides,
    'qvality'    : stats.QvalityDriver,
    'reassign'   : stats.ReassignmentDriver,
    'pout2tsv'   : converters.Pout2TSVDriver,
    }


command = commandmap[args.command](**vars(args))
command.run()

