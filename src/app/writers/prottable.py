from app.writers import tsv


def write_prottable(header, psms, outfn):
    tsv.write_tsv(header, psms, outfn)


def write_prottable_with_na(header, psms, outfn):
    na_psms = (tsv.parse_NA(x, header) for x in psms)
    write_prottable(header, na_psms, outfn)
