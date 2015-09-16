from statistics import median

from app.dataformats import mzidtsv as mzidtsvdata
from app.dataformats import peptable as peptabledata
from app.readers import tsv as reader
from app.actions.peptable.base import evaluate_peptide
from app.actions.headers.peptable import switch_psm_to_peptable_fields


def get_quantcols(pattern, oldheader, coltype):
    """Searches for quantification columns using pattern and header list.
    Calls reader function to do regexp. Returns either a single column for
    precursor quant, or a list of columns for isobaric quant."""
    if pattern is None:
       return False
    if coltype == 'isob':
        oldcols = reader.get_cols_in_file(pattern, oldheader)
        return {col: 'PSM median {0}'.format(col) for col in oldcols}
    elif coltype == 'precur':
        return reader.get_cols_in_file(pattern, oldheader, single_col=True)


def generate_peptides(tsvfn, oldheader, scorecol, isofieldmap,
                      precurquantcol, fncol=None, higherbetter=True):
    if fncol is None:
        fncol = mzidtsvdata.HEADER_SPECFILE
    peptides = {}
    switch_map = switch_psm_to_peptable_fields(oldheader)
    for psm in reader.generate_tsv_psms(tsvfn, oldheader):
        for oldkey, newkey in switch_map.items():
            try:
                psm[newkey] = psm.pop(oldkey)
            except KeyError:
                pass
        pepseq = psm[peptabledata.HEADER_PEPTIDE]
        peptides = evaluate_peptide(peptides, psm, pepseq, higherbetter,
                                    scorecol, fncol)
        add_quant_values(peptides, psm, isofieldmap, precurquantcol)
    for peptide in peptides.values():
        peptide['line'][peptabledata.HEADER_LINKED_PSMS] = '; '.join(
            peptide['psms'])
        for qtype, pepquant in peptide['quant'].items():
            peptide['line'].update(parse_quant_data(qtype, pepquant,
                                                    isofieldmap))
        yield peptide['line']


def parse_quant_data(qtype, pepquant, fieldmap=None):
    if qtype == 'isob':
        quants = {fieldmap[x]: get_median(pepquant[x]) for x in fieldmap}
    elif qtype == 'precur':
        quants = {peptabledata.HEADER_AREA: max([float(x) for x in pepquant])}
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


def add_quant_values(allpeps, psm, isobq_fieldmap, precurq_field):
    try:
        quants = allpeps[psm[peptabledata.HEADER_PEPTIDE]]['quant']
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
    allpeps[psm[peptabledata.HEADER_PEPTIDE]]['quant'] = quants
