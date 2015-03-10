from . import tsv


def write_mzid_tsv(header, psms, outfn):
    tsv.write_tsv(header, psms, outfn)


def write_multiple_mzidtsv(baseheader, psms, base_outfile):
    outfile_handles = {}
    for psm in psms:
        line = psm['psm']
        split_pool = psm['split_pool']
        header = baseheader.format(split_pool)
        if split_pool not in outfile_handles:
            outfile_handles[split_pool] = open(base_outfile.format(split_pool),
                                               'w')
            tsv.write_tsv_line_from_list(header, outfile_handles[split_pool])
        tsv.write_tsv_line_from_list([line[field] for field in header],
                                     outfile_handles[split_pool])
