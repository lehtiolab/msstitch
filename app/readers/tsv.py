def get_tsv_header(tsvfn):
    with open(tsvfn) as fp:
        return next(fp).strip().split('\t')


def get_header_index(header, index_options):
    """Returns index position in list header of the first found
    in the list index_options. If not found throws a ValueError"""
    for option in index_options:
        try:
            index = header.index(option)
        except ValueError:
            continue
        else:
            return index
    raise ValueError('Cannot find any of proposed header indices [ {0} ] in'
                     ' MzIdentML tsv header'.format(', '.join(index_options)))


def get_mzidtsv_lines_scannr_specfn(fn):
    """Returns generator of lines of tsv, skipping header, as split lists,
    and a tuple containing (spectra file, scan nr)."""
    with open(fn) as fp:
        header = next(fp).strip().split('\t')
        scan_ix = get_header_index(header, ['ScanNum', 'scannr',
                                            'scan_nr', 'Scan number'])
        fn_ix = get_header_index(header, ['#SpecFile', 'spectra_file',
                                          'specfile'])

        for line in fp:
            line = line.strip().split('\t')
            yield line, (line[fn_ix], line[scan_ix])
