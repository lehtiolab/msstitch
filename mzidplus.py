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
import app.drivers.mzidtsv.spectra as spectradrivers
import app.drivers.mzidtsv.percolator as percodrivers
import app.drivers.mzidtsv.proteingrouping as pgdrivers
import app.drivers.mzidtsv.quant as quantdrivers
import app.drivers.mzidtsv.splitmerge as splitmergedrivers


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
                    'spectratsv - Add spectra data such as retention time,\n'
                    'biological set name to tsv\n'
                    'percotsv       - Add percolator data to a  TSV with \n'
                    'MSGF+ output. Specify TSV file with -i, mzid file with \n'
                    '--mzid.\n\n'
                    'quanttsv       - Add quantitative data from openMS\n'
                    'consensusXML to a tab separated file with\n'
                    'PSMs. Needs to be passed a lookup db with --dbfile,\n'
                    'which has to contain quant information, and\n'
                    'optionally --isobaric, --precursor, --rttol, --mztol,\n'
                    '--spectracol changes the column where the spectra\n'
                    'file names are in from the standard #SpecFile column.\n\n'
                    'proteingroup   - Takes lookup SQLite result, uses it\n'
                    'to output mzidtsv file with protein groups\n'
                    'With flags --confidence-lvl, --confidence-col,\n'
                    '--confidence-better, --fasta, --dbfile, --spectracol\n\n'
                    'mergetsv       - Merges multiple TSV tables of MSGF+ \n'
                    'output. Make sure headers are same in all files.\n\n'
                    'splittsv       - Splits an MSGF TSV PSM table into\n'
                    'multiple new tables based. Use with flags --bioset or\n'
                    '--splitcol, and optionally --set-rename and\n'
                    '--rename-cols.',
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
parser.add_argument('--precursor', dest='precursor', help='Flag. Specifies\n'
                    'quant data from MS1 precursors should be added from\n'
                    'lookup DB to tsv file.',
                    action='store_const', const=True, default=False)
parser.add_argument('--isobaric', dest='isobaric', help='Flag. Specifies\n'
                    'quant data from isobaric tags should be added from\n'
                    'lookup DB to tsv file.',
                    action='store_const', const=True, default=False)
parser.add_argument('--rttol', dest='rttol', help='Specifies tolerance\n'
                    'in seconds for retention time when mapping MS1 feature\n'
                    'quant info to identifications in the PSM table.',
                    type=float)
parser.add_argument('--mztol', dest='mztol', help='Specifies tolerance\n'
                    'in mass-to-charge when mapping MS1 feature quant info\n'
                    'to identifications in the PSM table.', type=float)
parser.add_argument('--mztoltype', dest='mztoltype', help='Type of tolerance\n'
                    'in mass-to-charge when mapping MS1 feature quant info\n'
                    'to identifications in the PSM table. One of ppm, Da.',
                    type=lambda x: parser_value_in_list(parser, x, ['ppm',
                                                                    'Da']))
parser.add_argument('--spectracol', dest='speccol', help='Column number\n'
                    'in which spectra file names are, in case some framework\n'
                    'has changed the file names. First column number is 1.',
                    type=int, required=False)
parser.add_argument('--dbfile', dest='lookup', help='Lookup database in '
                    'SQLite format, to be created using mslookup.py.',
                    type=lambda x: parser_file_exists(parser, x))
parser.add_argument('--unroll', dest='unroll', help='Flag. The tsv input file '
                    'from Mzid2TSV contains either one PSM per line with all '
                    'the proteins of that shared peptide on the same line (not'
                    ' unrolled, default), or one PSM/protein match per line '
                    'where each protein from that shared peptide gets its own '
                    'line (unrolled).',
                    action='store_const', const=True, default=False)
parser.add_argument('--bioset', dest='bioset', help='Flag. When using '
                    'splittsv, this enables automatic splitting on '
                    'biological set names, for which a a column specifying '
                    'these must exist.',
                    action='store_const', const=True, default=False)
parser.add_argument('--splitcol', dest='splitcol', help='Column number on '
                    'which to split a TSV PSM table', type=int, required=False)
parser.add_argument('--rename-col', dest='renamecols', help='Column numbers '
                    'to rename with name of e.g. set used for splitting. '
                    'Rename pattern: setname_oldcolumnname.', nargs='+')


# not supported yet
#parser.add_argument('--allpsms', dest='allpsms', action='store_true',
#                    help='All PSMs from a single scan should be included, '
#                    'not only the best scoring one.')

args = parser.parse_args()

commandmap = {
    'spectratsv': spectradrivers.TSVSpectraDriver,
    'percotsv': percodrivers.MzidPercoTSVDriver,
    'mergetsv': splitmergedrivers.MzidTSVConcatenateDriver,
    'splittsv': splitmergedrivers.MzidTSVSplitDriver,
    'quanttsv': quantdrivers.TSVQuantDriver,
    'proteingroup': pgdrivers.ProteinGroupDriver,
}

command = commandmap[args.command](**vars(args))
command.run()
