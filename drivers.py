# need to install:
#genshi
#libhdf5-dev 
#python-dev
#NumPy 1.5 or newer (1.6 recommended)
#


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


def merge_unique_best_scoring_peptides(fns, score, mergedfn):
    """This function processes multiple percolator runs from fractions and
    filters out the best scoring peptides. It writes a single fraction with
    those peptides and ALL psms from all fractions. 
    
    Namespace and static xml come from first percolator file. 
    Make sure fractions are from same percolator run."""
    namespace = readers.get_namespace(fns[0])
    static_xml = readers.get_percolator_static_xml(fns[0], namespace)
    allpsms = readers.generate_psms_multiple_fractions(fns, namespace)
    uniquepeps = filtering.filter_unique_peptides(fns, score, namespace)
    features = {'psm': allpsms, 'peptide': uniquepeps}
    writers.write_percolator_xml(static_xml, features, mergedfn)


if __name__ == '__main__':
    # just for testing
    import time
    t = time.time()
    print 'lets roll'
    fns = [
    '/mnt/kalevalatmp/test_galaxydb/files/000/dataset_108_files/dataset_108.dat_task_0',
    '/mnt/kalevalatmp/test_galaxydb/files/000/dataset_108_files/dataset_108.dat_task_1'
    ]
    print merge_unique_best_scoring_peptides(fns, 'svm',
                    'merged.xml')
    print 'Took %s seconds' % str(time.time() - t)

