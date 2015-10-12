#!/usr/bin/python3

"""
mslookup - Creating SQLite lookups for internal and external use
"""

import argparse
import os
from app.drivers.mslookup import (spectra, quant, proteingroups, biosets,
                                  proteinquant, pepquant, psms, seqspace)


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

                    'psms - Loads PSM table into lookup. Important for\n'
                    'several steps later on, such as protein grouping and\n'
                    'PSM quantitation. PSM TSV table passed to -i, \n'
                    'With flags --unroll, --spectracol, and specify a FASTA\n'
                    'file with --fasta OR use --protcol to get proteins.\n\n'

                    'proteingrouplookup  - Groups proteins from mzid2tsv\n'
                    'output (single file passed to -i).\n\n'

                    'isoquant - Create lookup of isobaric quant data in\n'
                    'OpenMS consensusXML format. Use requires --spectra,\n'
                    '--dbfile with an sqlite lookup of spectra, and passing\n'
                    'a consensusXML file to -i.\n\n'

                    'ms1quant - Creates lookup of precursor quant data in \n'
                    'OpenMS featureXML format. Use requires --spectra,\n'
                    '--dbfile with an sqlite lookup of spectra, --quanttype\n'
                    'to determine quant output, --mztol, --rttol, --mztoltype\n'
                    'for tolerance specification, and passing \n'
                    'a featureXML or kronik file to -i\n\n'

                    'protquant - Creates lookup of protein quantification\n'
                    'data in tab separated format. Header should include\n'
                    'quantification channel names, and if possible the\n'
                    'number of peptides quantified for each protein in the\n'
                    'respective channels. Lookup should already include '
                    'proteins.\n\n'

                    'seqspace - Creates lookup DB from FASTA file for use\n'
                    'with e.g. pycolator.py -c filterknown. You may  use \n'
                    'flags --cutproline, --notrypsin and --ntermwildcards.',
                    required=True
                    )
parser.add_argument('-i', dest='infile', nargs='+',
                    help='The input files to create the lookup from. Can be\n'
                    'type mzML, consensusXML, featureXML, percolatoroutXML,\n'
                    'tab separated PSM table, protein table or FASTA file.\n'
                    'If order is important then it is taken from the input\n'
                    'order at the command line. FASTA input for seqspace \n'
                    'is a single file.',
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
                    'of protein/PSM table in which protein accessions are \n'
                    'stored. First column number is 1.',
                    type=int, required=False)
parser.add_argument('--fasta', dest='fasta', help='FASTA sequence database, '
                    'to use when extracting gene names to the PSM table from '
                    'proteins.', type=lambda x: parser_file_exists(parser, x))
parser.add_argument('--pepcol', dest='pepcol', help='Column number\n'
                    'of peptide table in which peptide sequences are \n'
                    'stored. First column number is 1.',
                    type=int, required=False)
parser.add_argument('--ms1quantcolpattern', dest='precursorquantcolpattern',
                    help='Unique text pattern to identify precursor quant \n'
                    'column in protein table.',
                    type=str, required=False)
parser.add_argument('--isobquantcolpattern', dest='isobquantcolpattern',
                    help='Unique text pattern to identify isobaric quant \n'
                    'column in protein table.',
                    type=str, required=False)
parser.add_argument('--psmnrcolpattern', dest='psmnrcolpattern',
                    help='Unique text pattern to identify number-of-psms \n'
                    'column in protein table.',
                    type=str, required=False)
parser.add_argument('--probcolpattern', dest='probcolpattern',
                    help='Unique text pattern to identify protein probability\n'
                    'column in protein table.',
                    type=str, required=False)
parser.add_argument('--fdrcolpattern', dest='fdrcolpattern',
                    help='Unique text pattern to identify protein FDR\n'
                    'column in protein table.',
                    type=str, required=False)
parser.add_argument('--pepcolpattern', dest='pepcolpattern',
                    help='Unique text pattern to identify protein PEP\n'
                    'column in protein table.',
                    type=str, required=False)
parser.add_argument('--quanttype', dest='quanttype',
                    help='Filetype of precursor quants to store. Choose from\n'
                    'kronik or openms.',
                    type=str, required=False)
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
parser.add_argument('--cutproline', dest='proline', help='With flag, trypsin is '
                    'considered to cut before a proline residue. The filter '
                    'known will filter against both cut and non-cut peptides.',
                    action='store_const', const=True, default=False)
parser.add_argument('--notrypsin', dest='notrypsin', help='With flag, no \n'
                    'trypsinization is performed. User is expected to deliver\n'
                    'pretrypsinized FASTA file.',
                    action='store_const', const=False, default=True)
parser.add_argument('--ntermwildcards', dest='falloff', help='With flag, the filter '
                    'known will filter against both intact peptides and those '
                    'that match to the C-terminal part of a tryptic peptide '
                    'from the database. Database should be built with this'
                    'flag in order for the lookup to work, since sequences'
                    'will be stored and looked up reversed',
                    action='store_const', const=True, default=False)


args = parser.parse_args()

commandmap = {
    'biosets': biosets.BioSetLookupDriver,
    'spectra': spectra.SpectraLookupDriver,
    'psms': psms.PSMLookupDriver,
    'isoquant': quant.IsobaricQuantLookupDriver,
    'ms1quant': quant.PrecursorQuantLookupDriver,
    'proteingrouplookup': proteingroups.ProteinGroupLookupDriver,
    'pepquant': pepquant.PeptideQuantLookupDriver,
    'protquant': proteinquant.ProteinQuantLookupDriver,
    'seqspace': seqspace.SeqspaceLookupDriver,
}

command = commandmap[args.command](**vars(args))
command.run()
