import os

import readers
import filtering
import writers

def create_outfilepath(fn, outdir, suffix=None):
    basefn = os.path.basename(fn)
    outfn = basefn + suffix
    return os.path.join(outdir, outfn)


def merge_multiple_fractions(fns):
    """Performs the work to merge parallelized percolator fractions.
    Target/decoy split, filtering unique peptides, running qvality on resulting
    score distributions for psms and peptides and setting values."""
    pass


def split_target_decoy(fns, outdir, targetsuffix='_target.xml', decoysuffix='_decoy.xml'):
    """ Calls splitter to split percolator output into target/decoy elements.
        Writes two new xml files with features. Currently only psms and
        peptides. Proteins not here, since one cannot do protein inference
        before having merged and remapped multifraction data anyway.
    """
    for fn in fns:
        namespace = readers.get_namespace(fn)
        static_xml = readers.get_percolator_static_xml(fn, namespace)
        split_elements = filtering.split_target_decoy(fn, namespace)
        targetfn = create_outfilepath(fn, outdir, targetsuffix)
        decoyfn = create_outfilepath(fn, outdir, decoysuffix)
        writers.write_percolator_xml(static_xml, split_elements['target'], 
                                        targetfn)
        writers.write_percolator_xml(static_xml, split_elements['decoy'],
                                        decoyfn)




def merge_unique_best_scoring_peptides(fns, outdir, score='svm',
                outsuffix='_merged.xml'):
    """This function processes multiple percolator runs from fractions and
    filters out the best scoring peptides. It writes a single fraction with
    those peptides and ALL psms from all fractions. 
    
    Namespace and static xml come from first percolator file. 
    Make sure fractions are from same percolator run."""
    namespace = readers.get_namespace(fns[0])
    static_xml = readers.get_percolator_static_xml(fns[0], namespace)
    allpsms = readers.generate_psms_multiple_fractions(fns, namespace)
    allpeps = readers.generate_psms_multiple_fractions(fns, namespace)
    uniquepeps = filtering.filter_unique_peptides(allpeps, score, namespace)
    features = {'psm': allpsms, 'peptide': uniquepeps}
    merged_fn = create_outfilepath(fns[0], outdir, outsuffix)
    writers.write_percolator_xml(static_xml, features, merged_fn)
