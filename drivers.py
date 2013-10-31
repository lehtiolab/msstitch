import readers
import filtering
import writers


def merge_multiple_fractions(fns):
    """Performs the work to merge parallelized percolator fractions.
    Target/decoy split, filtering unique peptides, running qvality on resulting
    score distributions for psms and peptides and setting values."""
    pass


def split_target_decoy(fn, targetfn='target.xml', decoyfn='decoy.xml'):
    """ Calls splitter to split percolator output into target/decoy elements.
        Writes two new xml files with features. Currently only psms and
        peptides. Proteins not here, since one cannot do protein inference
        before having merged and remapped multifraction data anyway.
    """
    namespace = readers.get_namespace(fn)
    static_xml = readers.get_percolator_static_xml(fn, namespace)
    split_elements = filtering.split_target_decoy(fn, namespace)
    writers.write_percolator_xml(static_xml, split_elements['target'], targetfn)
    writers.write_percolator_xml(static_xml, split_elements['decoy'], decoyfn)


def merge_filter_unique_peptides(fns, score):
    """Make sure fractions are from same
    percolator run."""
    psm_generators = []
    namespace = readers.get_namespace(fns[0])
    for fn in fns:
        psm_generators.append(readers.get_psms(fn, namespace))
    filtering.filter_unique_peptides(fns, score, namespace)

