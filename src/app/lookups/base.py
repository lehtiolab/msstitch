from app.lookups.sqlite import (proteingroups, quant, searchspace,
                                biosets, spectra, prottable, psms,
                                peptable)


def get_lookup(fn, lookuptype):
    lookupmap = {'biosets': biosets.BioSetDB,
                 'spectra': spectra.SpectraDB,
                 'psm': psms.PSMDB,
                 'proteingroups': proteingroups.ProteinGroupDB,
                 'quant': quant.QuantDB,
                 'isoquant': quant.IsobaricQuantDB,
                 'ms1quant': quant.PrecursorQuantDB,
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
