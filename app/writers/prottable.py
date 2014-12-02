from app.writers import tsv


def write_prottable(header, psms, outfn):
    tsv.write_tsv(header, psms, outfn)
