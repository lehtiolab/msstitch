from app.dataformats import mzidtsv as mzidtsvdata
from app.dataformats import prottable as prottabledata


def generate_master_proteins(psms, protcol):
    """Fed with a psms generator, this returns the master proteins present
    in the PSM table. PSMs with multiple master proteins are excluded."""
    master_proteins = {}
    if protcol is None:
        protcol = mzidtsvdata.HEADER_MASTER_PROT
    for psm in psms:
        protacc = psm[protcol]
        if ';' in protacc:
            continue
        master_proteins[protacc] = 1
    for protacc in master_proteins:
        yield {prottabledata.HEADER_PROTEIN: protacc}
