from . import tsv


def write_mzid_tsv(header, scanresults, outfn):
    """Writes psms from generator to a tsv file. PSMs are dicts containing
    a 'line' key with the line to write, and an optional 'rank'. Rank is
    added to the line in case it is not None, which it would be by including
    multiple PSMs per scan.
    """
    with open(outfn, 'w') as fp:
        tsv.write_tsv_line_from_list(header, fp)
        for scanresult in scanresults:
            for psm in scanresult:
                if psm['rank'] is not None:
                    line = [psm['rank']]
                else:
                    line = []
                line.extend(psm['line'])
                tsv.write_tsv_line_from_list(line, fp)
