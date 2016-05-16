from app.writers import tsv


def write_mzid_tsv(header, psms, outfn):
    tsv.write_tsv(header, psms, outfn)


def write_multi_mzidtsv(header, oldheader, psms, base_outfile):
    outfile_handles = {}
    for psm in psms:
        line = psm['psm']
        split_pool = psm['split_pool']
        if split_pool not in outfile_handles:
            outfile_handles[split_pool] = open(base_outfile.format(split_pool),
                                               'w')
            tsv.write_tsv_line_from_list(header, outfile_handles[split_pool])
        tsv.write_tsv_line_from_list([line[field] for field in oldheader],
                                     outfile_handles[split_pool])
    [handle.close() for handle in outfile_handles.values()]
