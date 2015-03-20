from app.lookups.sqlite import (proteingroups, quant, searchspace,
                                biosets, spectra, proteinquant)


def get_lookup(fn, lookuptype):
    lookupmap = {'biosets': biosets.BioSetDB,
                 'spectra': spectra.SpectraDB,
                 'proteingroups': proteingroups.ProteinGroupDB,
                 'quant': quant.QuantDB,
                 'searchspace': searchspace.SearchSpaceDB,
                 'proteinquant': proteinquant.ProtQuantDB,
                 'protgroupprottable': proteingroups.ProteinGroupProteinTableDB,
                 }
    return lookupmap[lookuptype](fn)


def create_new_lookup(fn, lookuptype):
    with open(fn, 'w'):
        pass
    return get_lookup(fn, lookuptype)
