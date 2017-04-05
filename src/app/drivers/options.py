import argparse

shared_options = {
    'fn': {'driverattr': 'fn', 'dest': 'infile', 'type': 'file', 'clarg': '-i',
           'help': 'Input file of {} format'},
    'outfile': {'driverattr': 'outfile', 'dest': 'outfile', 'type': str,
                'clarg': '-o', 'help': 'Output file', 'required': False},
    'outdir': {'driverattr': 'outdir', 'dest': 'outdir', 'clarg': '-d',
               'help': 'Directory to output in', 'type': 'file',
               'required': False},
    'multifiles': {'driverattr': 'fn', 'dest': 'infile', 'clarg': '-i',
                   'help': 'Multiple input files of {} format',
                   'type': 'file', 'nargs': '+'},
    'lookupfn': {'driverattr': 'lookupfn', 'clarg': '--dbfile',
                 'type': 'file', 'help': 'Database lookup file'},
    'setnames': {'driverattr': 'setnames', 'dest': 'setnames',
                 'type': str, 'nargs': '+', 'clarg': '--setnames',
                 'help': 'Names of biological sets. Can be '
                 'specified with quotation marks if spaces are '
                 'used'},
    'spectracol': {'driverattr': 'spectracol', 'dest': 'spectracol',
                   'type': int, 'clarg': '--spectracol', 'help':
                   'Column number in which spectra file names are, '
                   'in case some framework has changed the file '
                   'names. First column number is 1.', 'required': False,
                   'default': 1},
    'decoyfn': {'driverattr': 'decoyfn', 'dest': 'decoyfn',
                'help': 'Decoy input file (percolator out XML) for qvality',
                'type': 'file', 'clarg': '--decoyfn'},
    'proline': {'driverattr': 'proline', 'dest': 'proline', 'required': False,
                'clarg': '--cutproline', 'action': 'store_const',
                'const': True, 'default': False, 'help': 'Flag to make '
                'trypsin before a proline residue. Then filtering will be '
                'done against both cut and non-cut peptides.',
                },
    'fasta': {'driverattr': 'fasta', 'dest': 'fasta',
              'type': 'file', 'help': 'FASTA sequence database',
              'required': False, 'default': False, 'clarg': '--fasta'},
    'trypsinize': {'driverattr': 'trypsinize', 'dest': 'trypsinize',
                   'clarg': '--notrypsin', 'required': False,
                   'action': 'store_const', 'const': False, 'default': True,
                   'help': 'Do not trypsinize. User is expected to deliver a'
                   'pretrypsinized FASTA file'
                   },
    'featuretype': {'driverattr': 'featuretype', 'dest': 'featuretype',
                    'help': 'Feature type to use for qvality. Can either be '
                    'psm or peptide.', 'clarg': '--feattype',
                    'type': 'pick', 'picks': ['psm', 'peptide']},
    'unroll': {'driverattr': 'unroll', 'clarg': '--unroll', 'const': True,
               'action': 'store_const', 'default': False, 'help': 'PSM table '
               'from Mzid2TSV contains either one PSM per line with all '
               'the proteins of that shared peptide on the same line (not'
               ' unrolled, default), or one PSM/protein match per line '
               'where each protein from that shared peptide gets its own '
               'line (unrolled).', 'required': False},
    'genecentric': {'driverattr': 'genecentric', 'dest': 'genecentric',
                    'clarg': '--genecentric', 'type': 'pick',
                    'picks': ['assoc', 'genes'], 'required': False,
                    'default': False, 'help': 'Do not include protein group '
                    'data in output. Should be one of [genes, assoc]. '
                    'With assoc, associated gene IDs are used from e.g. '
                    'Biomart rather than the ones found in the FASTA db used '
                    'for PSM search. These need to have been stored when '
                    'creating a PSM lookup.'},
    'isobaric': {'driverattr': 'isobaric', 'clarg': '--isobaric',
                 'action': 'store_const', 'const': True, 'default': False,
                 'help': 'Specifies to add isobaric quant data from lookup DB '
                 'to output table', 'required': False,
                 },
    'precursor': {'driverattr': 'precursor', 'clarg': '--precursor',
                  'action': 'store_const', 'const': True, 'default': False,
                  'help': 'Specifies to add precursor quant data from lookup '
                  'DB to output table', 'required': False,
                  },
    'quantcolpattern': {'driverattr': 'quantcolpattern',
                        'clarg': '--isobquantcolpattern', 'type': str,
                        'default': None, 'required': False,
                        'help': 'Unique text pattern to identify '
                        'isobaric quant columns in input table.'},
    'precursorquantcolpattern': {'driverattr': 'precursorquantcolpattern',
                                 'type': str, 'required': False,
                                 'dest': 'precursorquantcolpattern',
                                 'clarg': '--ms1quantcolpattern',
                                 'default': None,
                                 'help': 'Unique text pattern to identify '
                                 'precursor quant column in input table.'},
    'quantacccolpattern': {'driverattr': 'quantacccolpattern',
                           'clarg': '--qaccpattern', 'type': str,
                           'help': 'Unique text pattern to identify '
                           'accession column in table containing quant info.'},
    'qvalityout': {'driverattr': 'qvalityout', 'dest': 'qvalityout',
                   'help': 'Qvality output file to fetch q-values and PEP '
                   'from', 'type': 'file', 'clarg': ['-q', '--qvality']},
    'proteincol': {'driverattr': 'proteincol', 'clarg': '--protcol',
                   'type': int, 'required': False, 'help': 'Column number in '
                   'table in which protein or gene accessions are. '
                   'stored. First column number is 1. Use in case of not '
                   'using standard {} column'},
    'pcolpattern': {'driverattr': 'pcolpattern', 'clarg': '--protcolpattern',
                    'type': str, 'required': False, 'help': 'Text pattern to '
                    'identify column in table in which protein or gene '
                    'accessions are. Use in case of not using standard '
                    '{} column', 'default': False},
    'fdrcolpattern': {'driverattr': 'fdrcolpattern', 'dest': 'fdrcolpattern',
                      'clarg': '--fdrcolpattern', 'type': str,
                      'required': False, 'default': None,
                      'help': 'Unique text pattern to identify '
                      'FDR column in input table.'},
    'fastadelim': {'driverattr': 'fastadelim', 'clarg': '--fastadelim',
                   'dest': 'fastadelim', 'required': False, 'type': 'pick',
                   'picks': ['tab', 'pipe', 'semicolon'],
                   'help': 'Delimiter in FASTA header, used to parse gene '
                   'names in case of non-ENSEMBL/Uniprot'},
    'genefield': {'driverattr': 'genefield', 'clarg': '--genefield',
                  'dest': 'genefield', 'required': False, 'type': int,
                  'help': 'Field nr (first=1) in FASTA that contains gene '
                  'name when using --fastadelim to parse the gene names'},
    'minlength': {'driverattr': 'minlength', 'dest': 'minlength', 'default': 0,
                  'help': 'Minimum length of peptide to be included',
                  'type': int, 'clarg': '--minlen', 'required': False},
}

mslookup_options = {
    'falloff': {'driverattr': 'falloff', 'dest': 'falloff',
                'clarg': '--insourcefrag', 'default': False,
                'action': 'store_const', 'const': True, 'help': 'Apply '
                'filter against both intact peptides and those '
                'that match to the C-terminal part of a tryptic peptide '
                'from the database, resulting from in-source fragmentation, '
                'where some amino acids will be missing from the N-terminus. '
                'Specify the max number of amino acids that may be missing. '
                'Database should be built with this '
                'flag in order for the lookup to work, since sequences '
                'will be stored and looked up reversed', 'required': False
                },
    'mapfn': {'driverattr': 'mapfn', 'dest': 'mapfn',
              'type': 'file', 'clarg': '--map',
              'required': False, 'help': 'File that contains '
              'a map obtained from ENSEMBL BioMart which '
              'should contain mappings from protein accession '
              'to Gene ENSG and Symbol.'},
    'decoy': {'driverattr': 'decoy', 'dest': 'decoy', 'clarg': '--decoy',
              'action': 'store_const', 'const': True,
              'default': False, 'help': 'Specifies lookup is '
              'for decoy PSMs, use with --map in case there '
              'are no decoy symbols in the FASTA used to '
              'search.', 'required': False},
    'spectrafns': {'driverattr': 'spectrafns', 'dest': 'spectra',
                   'type': str, 'help': 'Spectra files in mzML '
                   'format. Multiple files can be specified, if '
                   'order is important, e.g. when matching them '
                   'with quant data, the order will be their input '
                   'order at the command line.', 'clarg': '--spectra',
                   'nargs': '+'},
    'quantfiletype': {'driverattr': 'quantfiletype', 'dest': 'quanttype',
                      'clarg': '--quanttype', 'type': 'pick', 'help':
                      'Filetype of '
                      'precursor quants to store. One of kronik or openms.',
                      'picks': ['kronik', 'openms']},
    'rttol': {'driverattr': 'rt_tol', 'dest': 'rttol', 'clarg': '--rttol',
              'type': float, 'help': 'Specifies tolerance in seconds for '
              'retention time when mapping MS1 feature quant info to '
              'identifications in the PSM table.'},
    'mztol': {'driverattr': 'mz_tol', 'dest': 'mztol', 'clarg': '--mztol',
              'type': float, 'help': 'Specifies tolerance in mass-to-charge '
              'when mapping MS1 feature quant info to identifications in '
              'the PSM table.'},
    'mztoltype': {'driverattr': 'mz_toltype', 'dest': 'mztoltype',
                  'type': 'pick', 'picks': ['ppm', 'Da'],
                  'clarg': '--mztoltype',
                  'help': 'Type of tolerance in mass-to-charge when mapping '
                  'MS1 feature quant info to identifications in the PSM table.'
                  ' One of ppm, Da.'},
    'peptidecol': {'driverattr': 'peptidecol', 'dest': 'peptidecol',
                   'type': int, 'clarg': '--peptidecol', 'help':
                   'Column nr of peptide table where peptide sequences are '
                   'stored. First column is nr. 1'},
    'psmnrcolpattern': {'driverattr': 'psmnrcolpattern',
                        'dest': 'psmnrcolpattern',
                        'clarg': '--psmnrcolpattern',
                        'type': str, 'default': None, 'required': False,
                        'help': 'Unique text pattern to identify '
                        'number-of-psms column in input table.'},
    'probcolpattern': {'driverattr': 'probcolpattern',
                       'dest': 'probcolpattern', 'clarg': '--probcolpattern',
                       'type': str, 'required': False,
                       'default': None,
                       'help': 'Unique text pattern to identify '
                       'protein probability column in input table.'},
    'pepcolpattern': {'driverattr': 'pepcolpattern', 'dest': 'pepcolpattern',
                      'clarg': '--pepcolpattern', 'type': str,
                      'required': False, 'default': None,
                      'help': 'Unique text pattern to identify '
                      'protein PEP column in input table.'},
}
mslookup_options['proteincol'] = {k: v for k, v
                                  in shared_options['proteincol'].items()}
mslookup_options['proteincol'].update(
    {'default': 1, 'help':
     mslookup_options['proteincol']['help'].format('first')})
mslookup_options['fasta'] = {k: v for k, v in shared_options['fasta'].items()}
mslookup_options['fasta']['help'] = ('FASTA sequence database to use when '
                                     'extracting gene names to the PSM '
                                     'table from proteins')

pycolator_options = {
    'maxlength': {'driverattr': 'maxlength', 'dest': 'maxlength',
                  'default': None, 'required': False,
                  'help': 'Maximum length of peptide to be included in '
                  'filtered data.', 'type': int, 'clarg': '--maxlen'},
    'score': {'driverattr': 'score', 'dest': 'score', 'default': 'svm',
              'help': 'Score to filter unique peptides on, e.g. svm',
              'type': 'pick', 'picks': ['svm', 'q', 'pep', 'p'],
              'clarg': '--score', 'required': False},
    'qoptions': {'driverattr': 'qoptions', 'dest': 'qoptions',
                 'required': False, 'clarg': '--qoptions', 'default': None,
                 'help': 'Extra options that may be passed to qvality. '
                 'Option form: --qoptions ***flag value ***flag ***flag value',
                 'nargs': '+', 'type': str},
    'protheaders': {'driverattr': 'protheaders', 'clarg': '--protheaders',
                    'nargs': '+', 'type': str,
                    'help': 'Specify protein FASTA headers to split on. '
                    'Multiple headers of the same split-type can be grouped '
                    'with semicolons. E.g. --protheaders ENSP;sp '
                    'PSEUDOGEN;ncRNA would split into ENSP/swissprot peptides '
                    'and pseudogenes/non-coding RNA peptides.'},
    'deamidate': {'driverattr': 'deamidate', 'clarg': '--deamidate',
                  'action': 'store_const', 'default': False, 'const': True,
                  'help': 'Filter against both normal peptides and deamidated '
                  'peptides where a D->N transition has occurred.',
                  'required': False},
    'forcetryp': {'driverattr': 'forcetryp', 'clarg': '--enforce-tryptic',
                  'action': 'store_const', 'default': False, 'const': True,
                  'help': 'When filtering peptides against whole proteins, '
                  'filter out peptides that match a whole protein but only '
                  'if they are fully tryptic, i.e. the protein needs K,R 1 '
                  'position upstream of the peptide, and the peptide '
                  'C-terminal should also be K,R. Useful when discerning '
                  'from pseudogenes', 'required': False},
    'falloff': {'driverattr': 'falloff', 'dest': 'falloff',
                'clarg': '--insourcefrag',
                'type': int, 'default': 0, 'help': 'Apply '
                'filter against both intact peptides and those '
                'that match to the C-terminal part of a tryptic peptide '
                'from the database, resulting from in-source fragmentation, '
                'where some amino acids will be missing from the N-terminus. '
                'Specify the max number of amino acids that may be missing. '
                'Database should be built with this '
                'flag in order for the lookup to work, since sequences '
                'will be stored and looked up reversed', 'required': False
                },
}

mzidtsv_options = {
    'confcol': {'driverattr': 'confcol', 'clarg': '--confidence-col',
                'help': 'Confidence column number or name in the tsv file. '
                'First column has number 1.', 'type': int, 'required': False},
    'confpattern': {'driverattr': 'confpattern', 'clarg': '--confcolpattern',
                    'type': str, 'required': False, 'help': 'Text pattern to '
                    'identify column in table on which confidence filtering '
                    'should be done. Use in case of not using standard '
                    '{} column', 'default': False},
    'conflvl': {'driverattr': 'conflvl', 'clarg': '--confidence-lvl',
                'help': 'Confidence cutoff level as a floating point number',
                'type': float},
    'conftype': {'driverattr': 'conftype', 'clarg': '--confidence-better',
                 'help': 'Confidence type to define if higher or lower score '
                 'is better. One of [higher, lower]', 'type': 'pick',
                 'picks': ['higher', 'lower']},
    'medianpsms': {'driverattr': 'medianpsms', 'clarg': '--medianpsms',
                   'help': 'In case of using a separate PSM table with more '
                   'data to generate more robust medians (i.e. a superset of '
                   'input PSMs), specify that file here. The normalization '
                   'factors will be calculated on this file, and the PSMs '
                   'in the input will be adjusted using those factors rather '
                   'than factors derived from solely their own quantification '
                   'data', 'type': 'file', 'required': False},
    'normalize': {'driverattr': 'normalize', 'clarg': '--normalize',
                  'type': str, 'default': False, 'required': False,
                  'help': 'Normalization method for isobaric '
                  'quant data on protein or peptide level. Currently only '
                  'median centering is used. Use "--normalize median"'},
    'normalizeratios': {'driverattr': 'normalizeratios', 'type': 'file',
                        'clarg': '--norm-ratios', 'required': False,
                        'default': False, 'help': 'In case of using a '
                        'separate table to generate channel medians for '
                        'normalizing, specify that file here. The '
                        'normalization factors will be calculated from this '
                        'file, and the features in the input will be adjusted '
                        'using those factors rather than factors derived from '
                        'their own quantification data'},
    'targettable': {'driverattr': 'targettable', 'clarg': '--targettable',
                    'help': 'Table to output PSM or other feature quant data '
                    'to. Used when calculating PSM isobaric intenstity ratios '
                    'for proteins, peptides, genes. Leaving empty will output '
                    'to a new table, or when no --protcol is specified, '
                    'pastes ratios to the PSM table they are fetched from.',
                    'type': 'file', 'required': False},
    'mzidfn': {'driverattr': 'mzidfn', 'clarg': '--mzid', 'help': 'mzIdentML',
               'type': 'file'},
    'bioset': {'driverattr': 'bioset', 'clarg': '--bioset', 'const': True,
               'action': 'store_const', 'default': False,
               'help': 'this enables automatic splitting on '
               'biological set names, for which a a column specifying '
               'these must exist.', 'required': False},
    'splitcol': {'driverattr': 'splitcol', 'clarg': '--splitcol', 'type': int,
                 'help': 'Column number to split a PSM table on. First column '
                 'is number 1', 'required': False, 'default': None,
                 },
    'denompatterns': {'driverattr': 'denompatterns', 'required': False,
                      'clarg': '--denompatterns', 'type': str, 'nargs': '+',
                      'help': 'Regex patterns to detect denominator channels '
                      'when creating a PSM table with normalized ratios. If '
                      'both patterns and column numbers are given then column '
                      'numbers are used. Usage e.g. --denompattern _126 _131. '
                      'Also possible: --denompattern _12[6-7] to detect '
                      'multiple columns.'
                      },
    'denomcols': {'driverattr': 'denomcols', 'clarg': '--denomcols',
                  'type': int, 'nargs': '+', 'required': False,
                  'help': 'Column numbers of denominator channels when '
                  'creating a PSM table with normalized ratios',
                  },
    'minint': {'driverattr': 'minint', 'clarg': '--minint', 'type': float,
               'help': 'Intensity threshold of PSMs when calculating '
               'isobaric ratios. Values below threshold will be set to NA.',
               'required': False, 'default': -1,
               },
}
mzidtsv_options['quantcolpattern'] = {
    k: v for k, v in shared_options['quantcolpattern'].items()}
mzidtsv_options['quantcolpattern']['required'] = True

pepprottable_options = {
    # a mock infile to make sure we don't show or need an infile, e.g. in
    # case of building something from lookup
    'mock_infn': {'driverattr': 'fn', 'clarg': '-i', 'required': False,
                  'help': argparse.SUPPRESS, 'type': 'file'},
    'fdr': {'driverattr': 'fdr', 'clarg': '--fdr', 'action': 'store_const',
            'default': False, 'const': True, 'required': False,
            'help': 'Output FDR data to table'},
    'pep': {'driverattr': 'pep', 'clarg': '--pep', 'action': 'store_const',
            'default': False, 'const': True, 'required': False,
            'help': 'Output posterior error probabilities (PEP) to table.'},
    'quantfile': {'driverattr': 'quantfile', 'clarg': '--quantfile',
                  'type': 'file', 'help': 'File containing isobaric quant '
                  'data to add to table.'},
    'scorecolpattern': {'driverattr': 'scorecolpattern', 'type': str,
                        'clarg': '--scorecolpattern', 'help': 'Regular '
                        'expression pattern to find column where score '
                        'to use (e.g. percolator svm-score) is written.'},

}


prottable_options = {k: v for k, v in pepprottable_options.items()}
prottable_options['proteincol'] = {k: v for k, v
                                   in shared_options['proteincol'].items()}
prottable_options['proteincol'].update(
    {'default': False, 'help':
     prottable_options['proteincol']['help'].format('Master protein')})
prottable_options.update({
    'setname': {'driverattr': 'setname', 'clarg': '--setname', 'type': str,
                'help': 'Name of biological set to use when adding protein '
                'info to the table'},
    'probability': {'driverattr': 'probability', 'clarg': '--probability',
                    'action': 'store_const', 'default': False, 'const': True,
                    'help': 'Output protein probability to table.',
                    'required': False},
    'psmfile': {'driverattr': 'psmfile', 'clarg': '--psmtable', 'type': 'file',
                'help': 'PSM table file containing precursor quant data to '
                'add to table.'},
    'pepfile': {'driverattr': 'pepfile', 'clarg': '--peptable', 'type': 'file',
                'help': 'Peptide table file'},
    'minlogscore': {'driverattr': 'minlogscore', 'clarg': '--logscore',
                    'action': 'store_const', 'default': False, 'const': True,
                    'required': False, 'help': 'Score, e.g. q-values will '
                    'be converted to -log10 values.'},
    't_fasta': {'driverattr': 't_fasta', 'clarg': '--targetfasta',
                'type': 'file', 'help': 'FASTA file with target proteins '
                'to determine best scoring proteins of target/decoy pairs '
                'for pickqvality. In case using --picktype fasta',
                'required': False},
    'd_fasta': {'driverattr': 'd_fasta', 'clarg': '--decoyfasta',
                'type': 'file', 'help': 'FASTA file with decoy proteins '
                'to determine best scoring proteins of target/decoy pairs '
                'for pickqvality. In case using --picktype fasta',
                'required': False},
    'picktype': {'driverattr': 'picktype', 'clarg': '--picktype',
                 'type': 'pick', 'picks': ['fasta', 'result'],
                 'help': 'Feature type to use for qvality. Can be one of '
                 '[fasta, result].'},
    'mergecutoff': {'driverattr': 'mergecutoff', 'clarg': '--mergecutoff',
                    'type': float, 'default': False, 'help': 'FDR cutoff when '
                    'building merged protein table, to use when a cutoff has '
                    'been used before storing the table to lookup. FDR values '
                    'need to be stored in the lookup', 'required': False},
})


peptable_options = {k: v for k, v in pepprottable_options.items()}
peptable_options.update({
    'spectracol': {'driverattr': 'spectracol', 'dest': 'spectracol',
                   'type': int, 'clarg': '--spectracol', 'help':
                   'Specify this column number (first col. is 1) '
                   'containing PSM table spectrafiles (e.g. mzML) '
                   'if you want to track PSMs when creating peptide '
                   'tables', 'required': False},
    'noncentric': {'driverattr': 'noncentric', 'dest': 'noncentric',
                   'clarg': '--noncentric', 'action': 'store_const',
                   'const': True, 'default': False, 'required': False,
                   'help': 'Do not include protein group '
                   'or gene data in output, just use accessions. '},
    'genecentric': {'driverattr': 'genecentric', 'dest': 'genecentric',
                    'clarg': '--genecentric', 'action': 'store_const',
                    'const': True, 'default': False, 'required': False,
                    'help': 'Do not include protein group '
                    'data in output, but use gene names instead. '
                    'These need to have been stored when '
                    'creating a PSM lookup.'},
})
