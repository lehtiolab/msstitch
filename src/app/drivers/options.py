import argparse

shared_options = {
    'fn': {'driverattr': 'fn', 'type': 'file', 'clarg': '-i',
           'help': 'Input file of {} format'},
    'outfile': {'driverattr': 'outfile', 'type': str,
                'clarg': '-o', 'help': 'Output file', 'required': False},
    'outdir': {'driverattr': 'outdir', 'clarg': '-d',
               'help': 'Directory to output in', 'type': 'file',
               'required': False},
    'multifiles': {'driverattr': 'fn', 'clarg': '-i',
                   'help': 'Multiple input files of {} format',
                   'type': 'file', 'nargs': '+'},
    'lookupfn': {'driverattr': 'lookupfn', 'clarg': '--dbfile',
                 'type': 'file', 'help': 'Database lookup file'},
    'setnames': {'driverattr': 'setnames',
                 'type': str, 'nargs': '+', 'clarg': '--setnames',
                 'help': 'Names of biological sets. Can be '
                 'specified with quotation marks if spaces are '
                 'used'},
    'spectracol': {'driverattr': 'spectracol',
                   'type': int, 'clarg': '--spectracol', 'help':
                   'Column number in which spectra file names are, '
                   'in case some framework has changed the file '
                   'names. First column number is 1.', 'required': False,
                   'default': 1},
    'proline': {'driverattr': 'proline', 'required': False,
                'clarg': '--cutproline', 'action': 'store_const',
                'const': True, 'default': False, 'help': 'Flag to make '
                'trypsin before a proline residue. Then filtering will be '
                'done against both cut and non-cut peptides.',
                },
    'fasta': {'driverattr': 'fasta',
              'type': 'file', 'help': 'FASTA sequence database',
              'required': False, 'default': False, 'clarg': '--fasta'},
    'unroll': {'driverattr': 'unroll', 'clarg': '--unroll', 'const': True,
               'action': 'store_const', 'default': False, 'help': 'PSM table '
               'from Mzid2TSV contains either one PSM per line with all '
               'the proteins of that shared peptide on the same line (not'
               ' unrolled, default), or one PSM/protein match per line '
               'where each protein from that shared peptide gets its own '
               'line (unrolled).', 'required': False},
    'isobaric': {'driverattr': 'isobaric', 'clarg': '--isobaric',
                 'action': 'store_const', 'const': True, 'default': False,
                 'help': 'Specifies to add isobaric quant data from lookup DB '
                 'to output table', 'required': False,
                 },
    'precursor': {'driverattr': 'precursor', 'clarg': '--ms1quant',
                  'action': 'store_const', 'const': True, 'default': False,
                  'help': 'Specifies to add precursor quant data from lookup '
                  'DB to output table', 'required': False,
                  },
    'quantcolpattern': {'driverattr': 'quantcolpattern',
                        'clarg': '--isobquantcolpattern', 'type': str,
                        'default': False, 'required': False,
                        'help': 'Unique text pattern to identify '
                        'isobaric quant columns in input table.'},
    'precursorquantcolpattern': {'driverattr': 'precursorquantcolpattern',
                                 'type': str, 'required': False,
                                 'clarg': '--ms1quantcolpattern',
                                 'default': None,
                                 'help': 'Unique text pattern to identify '
                                 'precursor quant column in input table.'},
    'quantacccolpattern': {'driverattr': 'quantacccolpattern',
                           'clarg': '--qaccpattern', 'type': str,
                           'help': 'Unique text pattern to identify '
                           'accession column in table containing quant info.'},
    'featcol': {'driverattr': 'featcol', 'clarg': '--featcol',
                   'type': int, 'required': False, 'help': 'Column number in '
                   'table in which desired accessions are. '
                   'stored. First column number is 1. Use in case of not '
                   'using default column {}'},
    'featcolpattern': {'driverattr': 'featcolpattern', 'clarg': '--featcolpattern',
                    'type': str, 'required': False, 'help': 'Text pattern to '
                    'identify column in table in which protein or gene '
                    'accessions are. Use in case of not using standard '
                    '{} column', 'default': False},
    'fdrcolpattern': {'driverattr': 'fdrcolpattern',
                      'clarg': '--fdrcolpattern', 'type': str,
                      'default': None, 'help': 'Unique text pattern to identify '
                      'FDR column in input table.'},
    'flrcolpattern': {'driverattr': 'flrcolpattern', 'required': False,
                      'clarg': '--flrcolpattern', 'type': str,
                      'default': None, 'help': 'Unique text pattern to identify '
                      'FLR (peptide PTM false localization rate) column in input table.'},
    'fastadelim': {'driverattr': 'fastadelim', 'clarg': '--fastadelim',
                   'required': False, 'type': str,
                   'choices': ['tab', 'pipe', 'semicolon'],
                   'help': 'Delimiter in FASTA header, used to parse gene '
                   'names in case of non-ENSEMBL/Uniprot'},
    'genefield': {'driverattr': 'genefield', 'clarg': '--genefield',
                  'required': False, 'type': int,
                  'help': 'Field nr (first=1) in FASTA that contains gene '
                  'name when using --fastadelim to parse the gene names'},
    'minlength': {'driverattr': 'minlength', 'default': 0,
                  'help': 'Minimum length of peptide to be included',
                  'type': int, 'clarg': '--minlen', 'required': False},
    'addbioset': {'driverattr': 'addbioset',
            'clarg': '--addbioset', 'required': False, 'action': 'store_const',
            'default': False, 'const': True,
            'help': 'Add biological setname from DB lookup to PSM table',
               },
    'addmiscleav': {'driverattr': 'addmiscleav',
            'clarg': '--addmiscleav', 'required': False, 'action': 'store_const',
            'default': False, 'const': True, 'help': 'Add missed cleavages to PSM table',
               },
    'fullprotein': {'driverattr': 'fullprotein', 'clarg': '--fullprotein',
        'default': False, 'action': 'store_const', 'const': True, 'help':
        'Store full protein sequences (at a minimum-match length) in the '
        'SQLite file rather than tryptic sequences', 'required': False},
    'minint': {'driverattr': 'minint', 'clarg': '--minint', 'type': float,
               'help': 'Intensity threshold of PSMs when calculating '
               'isobaric ratios. Values below threshold will be set to NA. '
               'Defaults to no threshold.', 'required': False, 'default': 0,
               },
    'denompatterns': {'driverattr': 'denompatterns', 'required': False,
                      'clarg': '--denompatterns', 'type': str, 'nargs': '+',
                      'help': 'Regex patterns to detect denominator channels '
                      'in a PSM table, when creating a table with summarized '
                      'feature isobaric ratios. '
                      'If both --denompatterns and --denomcol are given then column '
                      'numbers are used. Usage e.g. --denompattern _126 _131. '
                      'Also possible: --denompattern _12[6-7] to detect '
                      'multiple columns.'
                      },
    'denomcols': {'driverattr': 'denomcols', 'clarg': '--denomcols',
                  'type': int, 'nargs': '+', 'required': False,
                  'help': 'Column numbers of denominator channels when '
                  'creating a summarized feature table with isobaric ratios '
                  'from PSMs'
                  },
    'mediansweep': {'driverattr': 'mediansweep', 'clarg': '--mediansweep',
            'action': 'store_const', 'const': True, 'default': False,
            'help': 'Instead of choosing a denominator channel, use the median intensity of each '
            'PSM as its denominator.', 'required': False,
            },
    'medianintensity': {'driverattr': 'medianintensity', 'clarg': '--medianintensity',
            'action': 'store_const', 'const': True, 'default': False,
            'help': 'Instead of choosing a denominator channel or median-sweeping, '
            'report the the median intensity of each summarized feat per channel. This '
            'results in reported intensities rather than ratios.', 'required': False,
            },
    'median_or_avg': {'driverattr': 'median_or_avg', 'clarg': '--summarize-average',
            'action': 'store_const', 'const': 'average', 'default': 'median',
            'required': False, 'help': 'Use average isobaric quantification values for summarizing '
            'quant from PSMs, instead of default PSM median values'},
    'keep_psms_na': {'driverattr': 'keepnapsms', 'clarg': '--keep-psms-na-quant',
            'action': 'store_const', 'const': True, 'default': False, 'required': False,
            'help': 'When summarizing isobaric quantification data, also use '
            'the PSMs that have an NA in any channel, even if these may contain '
            'overly noisy quant data in the other channels. Normally these PSMs '
            'would be skipped in quantification'},
    'logisoquant': {'driverattr': 'logisoquant', 'clarg': '--logisoquant',
        'required': False, 'action': 'store_const', 'const': True, 'help':
        'Output log2 values for isoquant ratios. This log2-transforms input PSM data '
        'prior to summarizing and optional normalization. Ratios will '
        'be calculated subtracted rather than divided, obviously.'},
    'mediannormalize': {'driverattr': 'mediannormalize',
        'clarg': '--median-normalize', 'default': False, 
        'required': False, 'action': 'store_const', 'const': True,
        'help': 'Median-centering normalization for isobaric quant data on protein or '
        'peptide level. This median-centers the data for each channel by '
        'dividing with the median channel value (or subtracting in case of '
        'log data), for each channel in output features, e.g. proteins.'},
    'ms1mediannormalize': {'driverattr': 'ms1mediannormalize',
        'clarg': '--median-normalize-ms1', 'default': False, 
        'required': False, 'action': 'store_const', 'const': True,
        'help': 'Median-centering normalization for MS1 quant data on protein or '
        'peptide level. This median-centers the data by dividing output MS1 quant '
        'values with the median output MS1 quant value '},
}

sequence_options = {
    'scramble': {
        'driverattr': 'scramble', 'clarg': '--scramble',
        'help': 'Decoy scrambling method, use: '
        '"tryp_rev": tryptic reverse, or "prot_rev": full (protein) reverse.',
        'required': False, 'default': 'tryp_rev'},
    'ignoretarget': {
        'driverattr': 'ignoretarget', 'clarg': '--ignore-target-hits',
        'help': 'Do not remove tryptic peptides from sequence where they match target DB',
        'required': False, 'action': 'store_const', 'const': True, 'default': False},
    'trypsinize': {'driverattr': 'trypsinize',
               'clarg': '--notrypsin', 'required': False,
               'action': 'store_const', 'const': False, 'default': True,
               'help': 'Do not trypsinize. User is expected to deliver a'
               'pretrypsinized FASTA file'
               },
    'max_shuffle': {'driverattr': 'max_shuffle',
               'clarg': '--maxshuffle', 'required': False, 'type': int, 'default': 10,
               'help': 'Amount of times to attempt to shuffle a decoy reversed peptide '
               'to make it not match target peptides, before discarding it.'
               ' Used when using tryptic peptide reversal (not protein reversal)'},
    'miss_cleavage': {'driverattr': 'miss_cleavage',
               'clarg': '--miscleav', 'required': False, 'type': int, 'default': 0,
               'help': 'Amount of missed cleavages to allow when trypsinizing, '
               'default is 0',
               },
        }


lookup_options = {
    'falloff': {'driverattr': 'falloff',
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
    'spectrafns': {'driverattr': 'spectrafns',
                   'type': str, 'help': 'Spectra files in mzML '
                   'format. Multiple files can be specified, if '
                   'order is important, e.g. when matching them '
                   'with quant data, the order will be their input '
                   'order at the command line.', 'clarg': '--spectra',
                   'nargs': '+'},
    'isobaric': {'driverattr': 'isobaricfns', 'clarg': '--isobaric', 
        'type': str, 'nargs': '+', 'required': False,
        'help': 'Isobaric quant output files from '
        'OpenMS in consensusXML '
        'format. Multiple files can be specified, '
        'and matching order with spectra files is important.',
        },
    'kronik': {'driverattr': 'kronikfns', 'clarg': '--kronik', 
        'type': str, 'nargs': '+', 'required': False,
        'help': 'MS1 persisting peptide quant output files from Kronik in text format.'
        'Multiple files can be specified, '
        'and matching order with spectra files is important.',
        },
    'dinosaur': {'driverattr': 'dinosaurfns', 'clarg': '--dinosaur', 
        'type': str, 'nargs': '+', 'required': False,
        'help': 'MS1 persisting peptide output files from Dinosaur in text format.'
        'Multiple files can be specified, '
        'and matching order with spectra files is important.',
        },
    'sum_or_apex': {'driverattr': 'sum_or_apex', 'clarg': '--apex',
        'action': 'store_const', 'const': 'apex', 'default': 'sum',
        'required': False, 'help': 'Use MS1 peak envelope apex instead of '
        'peak sum when storing quant data.'},
    'rttol': {'driverattr': 'rt_tol', 'clarg': '--rttol',
        'conditional_required': ['kronik'], 'type': float, 
        'help': 'Specifies tolerance in seconds for '
        'retention time when mapping MS1 feature quant info to '
        'identifications in the PSM table.'},
    'mztol': {'driverattr': 'mz_tol', 'clarg': '--mztol',
              'type': float, 'help': 'Specifies tolerance in mass-to-charge '
              'when mapping MS1 feature quant info to identifications in '
              'the PSM table.', 'conditional_required': ['kronik']},
    'mztoltype': {'driverattr': 'mz_toltype', 
                  'type': str, 'choices': ['ppm', 'Da'],
                  'clarg': '--mztoltype', 'conditional_required': ['kronik'],
                  'help': 'Type of tolerance in mass-to-charge when mapping '
                  'MS1 feature quant info to identifications in the PSM table.'
                  ' One of ppm, Da.'},
    'peptidecol': {'driverattr': 'peptidecol',
        'type': int, 'clarg': '--peptidecol', 'default': 1, 'required': False,
        'help': 'Column nr of peptide table where peptide sequences are '
        'stored. First and default column is nr. 1'},
}
lookup_options['featcol'] = {k: v for k, v 
        in shared_options['featcol'].items()}
lookup_options['featcol'].update(
        {'default': 1, 'help': lookup_options['featcol']['help'].format('1')})

lookup_options['fasta'] = {k: v for k, v in shared_options['fasta'].items()}
lookup_options['fasta']['help'] = ('FASTA sequence database to use when '
                                     'extracting gene names to the PSM '
                                     'table from proteins')

percolator_options = {
    'protheaders': {'driverattr': 'protheaders', 'clarg': '--protheaders',
                    'nargs': '+', 'type': str,
                    'help': 'Specify protein FASTA headers to split on. '
                    'Multiple headers of the same split-type can be grouped '
                    'with semicolons. E.g. --protheaders \'ENSP;sp '
                    'PSEUDOGEN;ncRNA\' would split into ENSP/swissprot peptides '
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
    'falloff': {'driverattr': 'falloff',
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

psmtable_options = {
    'genes': {'driverattr': 'genes', 'clarg': '--genes', 'action': 'store_const', 
        'const': True, 'default': False, 'help': 'Specifies to add genes to PSM '
        'table', 'required': False,
        },
    'proteingroup': {'driverattr': 'proteingroup', 'clarg': '--proteingroup',
        'action': 'store_const', 'const': True, 'default': False, 
        'help': 'Specifies to add protein groups to PSM table', 'required': False,
        },
    'oldpsmfile': {'driverattr': 'oldpsmfile', 'clarg': '--oldpsms', 'type': 'file',
                'help': 'PSM table file containing previously analysed PSMs to '
                'append new PSM table to.', 'required': False},
    'filtpep': {'driverattr': 'filtpep', 'clarg': '--filtpep',
                'help': 'Peptide q-value cutoff level as a floating point number',
                'type': float, 'required': False},
    'filtpsm': {'driverattr': 'filtpsm', 'clarg': '--filtpsm',
                'help': 'PSM q-value cutoff level as a floating point number',
                'type': float, 'required': False},
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
                 'is better. One of [higher, lower]', 'type': str,
                 'choices': ['higher', 'lower']},
    'minpurity': {'driverattr': 'min_purity', 'clarg': '--min-precursor-purity',
        'type': float, 'required': False, 'default': 0,
        'help': 'Minimum purity of precursor required to output isobaric '
        'quantification data, MS2 scans with purity below this value '
        'will be assigned NA in isobaric channels'},
    'medianpsms': {'driverattr': 'medianpsms', 'clarg': '--medianpsms',
                   'help': 'In case of using a separate PSM table with more '
                   'data to generate more robust medians (i.e. a superset of '
                   'input PSMs), specify that file here. The normalization '
                   'factors will be calculated on this file, and the PSMs '
                   'in the input will be adjusted using those factors rather '
                   'than factors derived from solely their own quantification '
                   'data', 'type': 'file', 'required': False},
    'percofn': {'driverattr': 'percofn', 'clarg': '--perco', 'help': 'Percolator '
        'XML output file', 'type': 'file'},
    'mzidfns': {'driverattr': 'mzidfns', 'clarg': '--mzids', 'help': 'MzIdentML '
        ' output files belonging to PSM table TSV files, use same order as for TSVs', 
        'type': 'file', 'nargs': '+'},
    'splitcol': {'driverattr': 'splitcol', 'clarg': '--splitcol',
                 'help': 'Either a column number to split a PSM table on, or '
                 '"TD", "bioset" for splitting on target/decoy or biological '
                 'sample set columns (resulting from msstitch perco2psm or '
                 'msstitch psmtable. First column is number 1.'
                 },
}
psmtable_options['quantcolpattern'] = {
    k: v for k, v in shared_options['quantcolpattern'].items()}
psmtable_options['quantcolpattern']['required'] = True
psmtable_options['featcol'] = {k: v for k, v 
        in shared_options['featcol'].items()}
psmtable_options['featcol'].update(
        {'help': 'Column number in table in which desired '
        'accessions to summarize are stored. First column number is 1. If not '
        'specified this will summarize to PSMs themselves, i.e. only calculate '
        'the ratios and add those to the PSM table.'})

pepprottable_options = {
    # a mock infile to make sure we don't show or need an infile, e.g. in
    # case of building something from lookup
    'mock_infn': {'driverattr': 'fn', 'clarg': '-i', 'required': False,
                  'help': argparse.SUPPRESS, 'type': 'file'},
    'fdr': {'driverattr': 'fdr', 'clarg': '--fdr', 'action': 'store_const',
            'default': False, 'const': True, 'required': False,
            'help': 'Output FDR data to table'},
    'quantfile': {'driverattr': 'quantfile', 'clarg': '--quantfile',
                  'type': 'file', 'help': 'File containing isobaric quant '
                  'data to add to table.'},
    'scorecolpattern': {'driverattr': 'scorecolpattern', 'type': str,
                        'clarg': '--scorecolpattern', 'help': 'Regular '
                        'expression pattern to find column where score '
                        'to use (e.g. percolator svm-score) is written.'},

}


prottable_options = {k: v for k, v in pepprottable_options.items()}
prottable_options.update({
    'psmfile': {'driverattr': 'psmfile', 'clarg': '--psmtable', 'type': 'file',
                'help': 'PSM table file containing isobaric quant data to '
                'add to table.', 'conditional_required': ['quantcolpattern']},
    'decoyfn': {'driverattr': 'decoyfn',
                'help': 'Decoy peptide table input file',
                'type': 'file', 'clarg': '--decoyfn'},
    'minlogscore': {'driverattr': 'minlogscore', 'clarg': '--logscore',
                    'action': 'store_const', 'default': False, 'const': True,
                    'required': False, 'help': 'Score, e.g. q-values will '
                    'be converted to -log10 values.'},
    't_fasta': {'driverattr': 't_fasta', 'clarg': '--targetfasta',
                'type': 'file', 'help': 'FASTA file with target proteins '
                'to determine best scoring proteins of target/decoy pairs '
                'for picked FDR. In case using --picktype ensg/genename',},
    'd_fasta': {'driverattr': 'd_fasta', 'clarg': '--decoyfasta',
                'type': 'file', 'help': 'FASTA file with decoy proteins '
                'to determine best scoring proteins of target/decoy pairs '
                'for picked FDR. In case using --picktype ensg/genename',},
    'picktype': {'driverattr': 'picktype', 'clarg': '--picktype', 'required': False,
                'type': str, 'choices': ['fasta', 'result'], 'default': 'fasta',
                 'help': 'Feature type to use for determining picked FDR. Can '
                 'be one of [ensg, genename, result]. "result" will infer T/D pairs '
                 'from the protein table and may not be as cleanly matched as '
                 'using ensg or genename options, which are based on being able to'
                 ' find genes in the fasta. Use result when genes are coming '
                 'from another source than your fasta input, or when you dont '
                 'have exact matching target/decoy paired fasta.'},
    'mergecutoff': {'driverattr': 'mergecutoff', 'clarg': '--mergecutoff',
                    'type': float, 'default': False, 'help': 'FDR cutoff when '
                    'building merged protein table, to use when a cutoff has '
                    'been used before storing the table to lookup. FDR values '
                    'need to be stored in the lookup', 'required': False},
})


peptable_options = {k: v for k, v in pepprottable_options.items()}
peptable_options.update({
    'spectracol': {'driverattr': 'spectracol',
                   'type': int, 'clarg': '--spectracol', 'help':
                   'Specify this column number (first col. is 1) '
                   'containing PSM table spectrafiles (e.g. mzML) '
                   'if you want to track PSMs when creating peptide '
                   'tables', 'required': False},
    'nogroup': {'driverattr': 'nogroup',
                   'clarg': '--no-group-annotation', 'action': 'store_const',
                   'const': True, 'default': False, 'required': False,
                   'help': 'For peptide table merging. Do not include protein group '
                   'or gene data in output, just use accessions. '},
    'genecentric': {'driverattr': 'genecentric',
                    'clarg': '--genecentric', 'action': 'store_const',
                    'const': True, 'default': False, 'required': False,
                    'help': 'For peptide table merging. Do not include protein group '
                    'data in output, but use gene names instead to count peptides '
                    'per feature, determine peptide-uniqueness.',
                    },
    'modelqvals': {'driverattr': 'modelqvals',
                    'clarg': '--modelqvals', 'action': 'store_const',
                    'const': True, 'default': False, 'required': False,
                    'help': 'Create linear-modeled q-vals for peptides, to avoid '
                    'overlapping stepped low-qvalue data of peptides with '
                    'different scores', },
    'qvalthreshold': {'driverattr': 'qvalthreshold',
        'type': float, 'clarg': '--qvalthreshold', 'help': 'Specifies the '
        'inclusion threshold for q-values to fit a linear model to. Any scores/'
        'q-values below this threshold will not be used.', 'default': 10e-4,
        'required': False},
    'minpeptidenr': {'driverattr': 'minpeptidenr',
        'type': int, 'clarg': '--minpepnr', 'default': 10, 'help': 'Specifies '
        'the minimal amount of peptides (passing the --qvalthreshold) needed '
        'to fit a linear model, default is 10.', 'required': False},
    'totalprotfn': {'driverattr': 'totalprotfn', 'type': 'file', 
        'clarg': '--totalproteome', 'required': False,
        'help': 'File containing total proteome quantification to normalize PTM '
        'peptide quantification against, i.e. Phospho peptides isobaric quant '
        'ratios are divided by their proteins to distinguish differential '
        'phosphorylation from the protein expression variation in the sample. '
        'This file can also be a gene names or ENSG table. Accession should be '
        'in the first column. The file is preferably generated from a search '
        'without the relevant PTM, and should not be normalized to channel '
        'medians'},
})
