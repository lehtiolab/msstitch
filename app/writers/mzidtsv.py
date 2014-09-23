from . import tsv


def write_mzid_tsv(header, psms, outfn):
    tsv.write_tsv(header, psms, outfn)
