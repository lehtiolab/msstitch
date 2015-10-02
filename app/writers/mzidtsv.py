from . import tsv


def write_mzid_tsv(header, psms, outfn):
    tsv.write_tsv(header, psms, outfn)


def write_multi_mzidtsv(header, oldheader, psms, setnames, base_outfile):
    outfile_handles = {}
    for i, setname in enumerate(setnames):
        outfile_handles[setname] = open(base_outfile.format(str(i)), 'w')
        tsv.write_tsv_line_from_list(header, outfile_handles[setname])
    for psm in psms:
        line = psm['psm']
        split_pool = psm['split_pool']
        tsv.write_tsv_line_from_list([line[field] for field in oldheader],
                                     outfile_handles[split_pool])
