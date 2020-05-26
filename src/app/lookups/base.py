from app.lookups.sqlite import (proteingroups, quant, searchspace,
                                biosets, spectra, prottable, psms,
                                peptable)


def get_lookup(fn, lookuptype):
    lookupmap = {'spectra': spectra.SpectraDB,
                 'specquant': quant.QuantDB,
                 'psm': psms.PSMDB,
                 'proteingroups': proteingroups.ProteinGroupDB,
                 'searchspace': searchspace.SearchSpaceDB,
                 'peptidetable': peptable.PepTableProteinCentricDB,
                 'peptidegenecentrictable': peptable.PepTableGeneCentricDB,
                 'peptidetableplain': peptable.PepTablePlainDB,
                 'prottable': prottable.ProtTableDB,
                 'genetable': prottable.GeneTableDB,
                 'associdtable': prottable.GeneTableAssocIDsDB,
                 }
    return lookupmap[lookuptype](fn)


def create_new_lookup(fn, lookuptype):
    with open(fn, 'w'):
        pass
    return get_lookup(fn, lookuptype)
