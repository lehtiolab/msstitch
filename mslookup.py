#!/usr/bin/python3

"""
mslookup - Creating SQLite lookups for internal and external use
"""

import argparse
import os
from app.drivers.mslookup import (spectra, quant, proteingroups, biosets,
                                  proteinquant)


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
                    'biosets - Create SQLite lookup of mzML input files and\n'
                    'biological set names. Input files are passed to -i,\n'
                    'respective set names are passed to --setnames\n\n'

                    'spectra - Create lookup of spectra in mzML \n'
                    'format. Requires passing mzML files to -i, but\n'
                    'neither --spectra nor --dbfile. Biological set names\n'
                    'for each file should be specified using\n'
                    '--setnames\n\n'

                    'proteingrouplookup  - Groups proteins from mzid2tsv\n'
                    'output (single file passed to -i). With flags \n'
                    '--confidence-lvl, --confidence-col,\n'
                    '--confidence-better, --fasta, --spectracolumn\n\n'

                    'isoquant - Create lookup of isobaric quant data in\n'
                    'OpenMS consensusXML format. Use requires --spectra,\n'
                    '--dbfile with an sqlite lookup of spectra, and passing\n'
                    'a consensusXML file to -i.\n\n'

                    'ms1quant - Creates lookup of precursor quant data in \n'
                    'OpenMS featureXML format. Use requires --spectra,\n'
                    '--dbfile with an sqlite lookup of spectra, --quanttype\n'
                    'to determine quant output and passing \n'
                    'a featureXML or kronik file to -i\n\n'

                    'protquant - Creates lookup of protein quantification\n'
                    'data in tab separated format. Header should include\n'
                    'quantification channel names, and if possible the\n'
                    'number of peptides quantified for each protein in the\n'
                    'respective channels. Lookup should already include '
                    'proteins.',
                    required=True
                    )
parser.add_argument('-i', dest='infile', nargs='+',
                    help='The input files to create the lookup from. Can be\n'
                    'type mzML, consensusXML, featureXML, percolatoroutXML,\n'
                    'tab separated PSM table or protein table. If order is\n'
                    'important then it is taken from the input order at the\n'
                    'command line.',
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
parser.add_argument('--setnames', dest='setnames', help='Names of biological\n'
                    'sets. Can be specified with quotation marks if spaces\n'
                    'are used', nargs='+')
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
parser.add_argument('--spectracol', dest='speccol', help='Column number\n'
                    'in which spectra file names are, in case some framework\n'
                    'has changed the file names. First column number is 1.',
                    type=int, required=False)
parser.add_argument('--unroll', dest='unroll', help='Flag. The tsv input file '
                    'from Mzid2TSV contains either one PSM per line with all '
                    'the proteins of that shared peptide on the same line (not'
                    ' unrolled, default), or one PSM/protein match per line '
                    'where each protein from that shared peptide gets its own '
                    'line (unrolled).',
                    action='store_const', const=True, default=False)
parser.add_argument('--protcol', dest='protcol', help='Column number\n'
                    'of protein table in which protein accessions are \n'
                    'stored. First column number is 1.',
                    type=int, required=False)
parser.add_argument('--quantcolpattern', dest='quantcolpattern',
                    help='Unique text pattern to identify quant column in \n'
                    'protein table.',
                    type=str, required=False)
parser.add_argument('--quanttype', dest='quanttype',
                    help='Filetype of precursor quants to store. Choose from\n'
                    'kronik or openms.',
                    type=str, required=False)

args = parser.parse_args()

commandmap = {
    'biosets': biosets.BioSetLookupDriver,
    'spectra': spectra.SpectraLookupDriver,
    'isoquant': quant.IsobaricQuantLookupDriver,
    'ms1quant': quant.PrecursorQuantLookupDriver,
    'proteingrouplookup': proteingroups.ProteinGroupLookupDriver,
    'protquant': proteinquant.ProteinQuantLookupDriver,
}

command = commandmap[args.command](**vars(args))
command.run()
