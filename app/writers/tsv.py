def write_tsv(headerfields, features, outfn):
    """Writes header and generator of lines to tab separated file.

    headerfields - list of field names in header in correct order
    features - generates 1 list per line that belong to header
    outfn - filename to output to. Overwritten if exists
    """
    with open(outfn, 'w') as fp:
        write_tsv_line_from_list(headerfields, fp)
        for line in features:
            write_tsv_line_from_list(line, fp)


def write_tsv_line_from_list(linelist, outfp):
    """Utility method to convert list to tsv line with carriage return"""
    line = '\t'.join(linelist)
    outfp.write(line)
    outfp.write('\n')


def write_quantpsm_tsv(header, psms):
    write_tsv(header, psms)
