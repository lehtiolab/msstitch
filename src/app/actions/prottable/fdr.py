"""Assigns FDR to proteins from qvality"""
from app.actions.pycolator import reassign as pyreassign
from app.dataformats import prottable as prottabledata


def assign_protein_fdr(qvalityfn, proteins, headerfields, scorefield):
    qvalityout = pyreassign.parse_qvality_output(qvalityfn)
    fdrheader = headerfields['proteinfdr'][prottabledata.HEADER_QVAL][None]
    pepheader = headerfields['proteinpep'][prottabledata.HEADER_PEP][None]
    for protein in proteins:
        outprotein = {k: v for k, v in protein.items()}
        try:
            score = round(float(outprotein[scorefield]), 5)
        except ValueError:
            qval, pep = 'NA', 'NA'
        else:
            qval, pep, warning = pyreassign.lookup_statistic(score, qvalityout)
        outprotein.update({fdrheader: qval, pepheader: pep})
        yield outprotein
