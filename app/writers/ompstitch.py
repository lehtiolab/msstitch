def write_tsv(headerfields, features, outfn):
    """Writes dict of features to tab separated format.
    headerfields - list of field names in header in correct order
    features - dict of features with keys == headerfields
    outfn - filename to output to. Overwritten if exists
    """
    with open(outfn, 'w') as fp:
        write_tsv_line_from_list(headerfields, fp)
        for feat in features:
            line = [features[x] for x in headerfields]
            write_tsv_line_from_list(line, fp)


def write_tsv_line_from_list(linelist, outfp):
    """Utility method to convert list to tsv line with carriage return"""
    line = linelist.join('\t')
    outfp.write(line)
    outfp.write('\n')


def write_quantpsm_tsv(header, psms):
    write_tsv(header, psms)

