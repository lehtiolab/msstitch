from app.dataformats import mzidtsv as mzidtsvdata


def generate_master_proteins(psms):
    """Fed with a psms generator, this returns the master proteins present
    in the PSM table. PSMs with multiple master proteins are excluded."""
    master_proteins = {}
    for psm in psms:
        protacc = psm[mzidtsvdata.HEADER_MASTER_PROT]
        if ';' in protacc:
            continue
        master_proteins[protacc] = 1
    master_proteins = (x for x in master_proteins)
    return master_proteins
