from . import tsv


def write_mzid_tsv(header, psms, outfn):
    tsv.write_tsv(header, psms, outfn)


def write_multi_mzidtsv(baseheader, oldheader, psms, base_outfile):
    outfile_handles = {}
    for psm in psms:
        line = psm['psm']
        split_pool = psm['split_pool']
        if split_pool not in outfile_handles:
            outfile_handles[split_pool] = open(base_outfile.format(split_pool),
                                               'w')
            header = format_split_header(baseheader, split_pool)
            tsv.write_tsv_line_from_list(header, outfile_handles[split_pool])
        tsv.write_tsv_line_from_list([line[field] for field in oldheader],
                                     outfile_handles[split_pool])


def format_split_header(baseheader, split_pool):
    """Returns header with split pool name in header fields where
    applicable."""
    return [x.format(split_pool) for x in baseheader]
