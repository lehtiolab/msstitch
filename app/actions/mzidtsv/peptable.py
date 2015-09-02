from statistics import median

from app.dataformats import mzidtsv as mzidtsvdata
from app.readers import tsv as reader


def get_quantcols(pattern, oldheader, coltype):
    """Searches for quantification columns using pattern and header list.
    Calls reader function to do regexp. Returns either a single column for
    precursor quant, or a list of columns for isobaric quant."""
    if coltype == 'isob':
        oldcols = reader.get_cols_in_file(pattern, oldheader)
        return {col: 'PSM median {0}'.format(col) for col in oldcols}
    elif coltype == 'precur':
        return reader.get_cols_in_file(pattern, oldheader, single_col=True)


def get_peptable_header(oldheader, isobqfieldmap=False, precurqfield=False):
    header = oldheader[:]
    if isobqfieldmap:
        for field, medianfield in isobqfieldmap.items():
            header = [medianfield if x==field else x for x in header]
    if precurqfield:
        header = [mzidtsvdata.HEADER_PRECURSOR_QUANT_MEDIAN if x==precurqfield else x for x in header]
    peptable_header = [mzidtsvdata.HEADER_LINKED_PSMS]
    ix = header.index(mzidtsvdata.HEADER_PEPTIDE)
    return header[:ix] + peptable_header + header[ix:]


def generate_peptides(tsvfn, oldheader, scorecol, isofieldmap,
                      precurquantcol, fncol=None, higherbetter=True):
    if fncol is None:
        fncol = mzidtsvdata.HEADER_SPECFILE
    peptides = {}
    for psm in reader.generate_tsv_psms(tsvfn, oldheader):
        try:
            existing_score = peptides[psm[mzidtsvdata.HEADER_PEPTIDE]]['score']
        except KeyError:
            add_peptide(peptides, psm, fncol, scorecol, new=True)
        else:
            if higherbetter and psm[scorecol] > existing_score:
                add_peptide(peptides, psm, fncol, scorecol)
            elif not higherbetter and psm[scorecol] < existing_score:
                add_peptide(peptides, psm, fncol, scorecol)
        finally:
            add_quant_values(peptides, psm, isofieldmap, precurquantcol)
    for peptide in peptides.values():
        peptide['line'][mzidtsvdata.HEADER_LINKED_PSMS] = '; '.join(
            peptide['psms'])
        for qtype, pepquant in peptide['quant'].items():
            peptide['line'].update(parse_quant_data(qtype, pepquant,
                                                    isofieldmap))
        yield peptide['line']


def parse_quant_data(qtype, pepquant, fieldmap=None):
    if qtype == 'isob':
        quants = {fieldmap[x]: get_median(pepquant[x]) for x in fieldmap}
    elif qtype == 'precur':
        quants = {mzidtsvdata.HEADER_PRECURSOR_QUANT_MEDIAN: get_median(
            pepquant)}
    return quants


def get_median(quantdata):
    """Parses lists of quantdata and gets the median from them. Stirps NA"""
    quantfloats = []
    for q in quantdata:
        try:
            quantfloats.append(float(q))
        except(TypeError, ValueError):
            pass
    if not quantfloats:
        return 'NA'
    return str(median(quantfloats))


def add_peptide(allpeps, psm, fncol, scorecol=False, new=False):
    peptide = {'score': psm[scorecol],
               'line': psm,
               'psms': []
               }
    if not new:
        peptide['psms'] = allpeps[psm[mzidtsvdata.HEADER_PEPTIDE]]['psms']
    peptide['psms'].append('{0}_{1}'.format(psm[fncol],
                                            psm[mzidtsvdata.HEADER_SCANNR]))
    allpeps[psm[mzidtsvdata.HEADER_PEPTIDE]] = peptide


def add_quant_values(allpeps, psm, isobq_fieldmap, precurq_field):
    try:
        quants = allpeps[psm[mzidtsvdata.HEADER_PEPTIDE]]['quant']
    except KeyError:
        quants = {}
    if isobq_fieldmap:
        if 'isob' not in quants:
            quants['isob'] = {x: [psm[x]] for x in isobq_fieldmap}
        else:
            for field in isobq_fieldmap:
                quants['isob'][field].append(psm[field])
    if precurq_field:
        try:
            quants['precur'].append(psm[precurq_field])
        except KeyError:
            quants['precur'] = [psm[precurq_field]]
    allpeps[psm[mzidtsvdata.HEADER_PEPTIDE]]['quant'] = quants
