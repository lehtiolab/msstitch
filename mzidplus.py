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
import app.drivers.mzidtsv.percolator as percodrivers
import app.drivers.mzidtsv.proteingrouping as pgdrivers
import app.drivers.mzidtsv.quant as quantdrivers
import app.drivers.mzidtsv.merge as mergedrivers


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
                    'percotsv       - Add percolator data to a  TSV with \n'
                    'MSGF+ output. Specify TSV file with -i, mzid file with \n'
                    '--mzid.\n'
                    'mergetsv       - Merges multiple TSV tables of MSGF+ \n'
                    'output. Make sure headers are same in all files.\n'
                    'quanttsv       - Add quantitative data from openMS\n'
                    'consensusXML to a tab separated file with\n'
                    'PSMs, their statistics and their quantification data.\n'
                    'Needs except these also the corresponding mzML spectra\n'
                    'files to correlate retention time to scan nrs.\n'
                    'proteingrouplookup  - Groups proteins from mzid2tsv\n'
                    'output. With flags --confidence-lvl, --confidence-col,\n'
                    '--confidence-better, --fasta\n'
                    'proteingroup   - Takes lookup SQLite result from \n'
                    'proteingrouplookup, uses it to output mzidtsv file with\n'
                    'protein groups. Same flags as proteingrouplookup, and\n'
                    '--protgroupdb \n',
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
parser.add_argument('--workdir', dest='workdir',
                    help='Working directory to output temporary files in. '
                    'Temporary files will be removed after finishing.'
                    )
parser.add_argument('--mzid', dest='mzid', help='mzIdentML file',
                    type=lambda x: parser_file_exists(parser, x))
parser.add_argument('--quants', dest='quants', help='Quants from OpenMS in '
                    'consXML format', nargs='+',
                    type=lambda x: parser_file_exists(parser, x))
parser.add_argument('--spectra', dest='spectra', help='mzML files', nargs='+',
                    type=lambda x: parser_file_exists(parser, x))
parser.add_argument('--fasta', dest='fasta', help='FASTA sequence database, '
                    'to optionally use with proteingrouping to enable sorting '
                    'on coverage, and in case of UNIPROT FASTA, evidence '
                    'levels.', type=lambda x: parser_file_exists(parser, x))
parser.add_argument('--confidence-col', dest='confcol', help='Confidence '
                    'column number or name in the tsv file. First column has'
                    ' number 1.')
parser.add_argument('--confidence-lvl', dest='conflvl', help='Confidence '
                    'cutoff level as a floating point number', type=float)
parser.add_argument('--confidence-better', dest='conftype', help='Confidence '
                    'type to define if higher or lower score is better. One '
                    'of [higher, lower]',
                    type=lambda x: parser_value_in_list(parser, x, ['higher',
                                                                    'lower']))
parser.add_argument('--protgroupdb', dest='protgroupdb', help='Protein group '
                    'lookup database in SQLite format. Can be created using '
                    'mzidplus.py command.',
                    type=lambda x: parser_file_exists(parser, x))
parser.add_argument('--unroll', dest='unroll', help='Flag. The tsv input file '
                    'from Mzid2TSV contains either one PSM per line with all '
                    'the proteins of that shared peptide on the same line (not'
                    ' unrolled, default), or one PSM/protein match per line '
                    'where each protein from that shared peptide gets its own '
                    'line (unrolled).',
                    action='store_const', const=True, default=False)
# not supported yet
#parser.add_argument('--allpsms', dest='allpsms', action='store_true',
#                    help='All PSMs from a single scan should be included, '
#                    'not only the best scoring one.')

args = parser.parse_args()

commandmap = {
    'percotsv': percodrivers.MzidPercoTSVDriver,
    'mergetsv': mergedrivers.MzidTSVConcatenateDriver,
    'quanttsv': quantdrivers.TSVQuantDriver,
    'proteingroup': pgdrivers.ProteinGroupDriver,
    'proteingrouplookup': pgdrivers.ProteinGroupLookupDriver,
}

command = commandmap[args.command](**vars(args))
command.run()
