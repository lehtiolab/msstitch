from statistics import median

from app.dataformats import mzidtsv as mzidtsvdata
from app.dataformats import peptable as peptabledata
from app.readers import tsv as reader
from app.actions.peptable.base import evaluate_peptide
from app.actions.headers.peptable import switch_psm_to_peptable_fields


def get_quantcols(pattern, oldheader, coltype):
    """Searches for quantification columns using pattern and header list.
    Calls reader function to do regexp. Returns a single column for
    precursor quant."""
    if pattern is None:
       return False
    if coltype == 'precur':
        return reader.get_cols_in_file(pattern, oldheader, single_col=True)


def generate_peptides(tsvfn, oldheader, scorecol, precurquantcol, fncol=None,
                      higherbetter=True):
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
        add_quant_values(peptides, psm, precurquantcol)
    for peptide in peptides.values():
        peptide['line'][peptabledata.HEADER_LINKED_PSMS] = '; '.join(
            peptide['psms'])
        for qtype, pepquant in peptide['quant'].items():
            peptide['line'].update(parse_quant_data(qtype, pepquant))
        yield peptide['line']


def parse_quant_data(qtype, pepquant):
    if qtype == 'precur':
        quants = {peptabledata.HEADER_AREA: get_peptide_quant(pepquant, 'precur')}
    return quants


def get_peptide_quant(quantdata, quanttype):
    """Parses lists of quantdata and returns maxvalue from them. Strips NA"""
    parsefnx = {'precur': max}
    quantfloats = []
    for q in quantdata:
        try:
            quantfloats.append(float(q))
        except(TypeError, ValueError):
            pass
    if not quantfloats:
        return 'NA'
    return str(parsefnx[quanttype](quantfloats))


def add_quant_values(allpeps, psm, precurq_field):
    try:
        quants = allpeps[psm[peptabledata.HEADER_PEPTIDE]]['quant']
    except KeyError:
        quants = {}
    if precurq_field:
        try:
            quants['precur'].append(psm[precurq_field])
        except KeyError:
            quants['precur'] = [psm[precurq_field]]
    try:
        allpeps[psm[peptabledata.HEADER_PEPTIDE]]['quant'] = quants
    except KeyError:
        # PSM sequence does not exist in allpeps because has no good-scoring peptides
        pass
