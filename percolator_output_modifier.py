# dependencies: pip install genshi
import sys
import os
import tarfile
from genshi.template import MarkupTemplate
from genshi.core import Markup

""" - Splits Percolator output into decoy and target files.
    - Extracts unique PSM/peptides/proteins out of a Percolator output file.
    - Merges Percolator output files
    Usage: python percolator_output_modifier.py command psm/peptides/proteins [score] infile outfile [outfile2]
    """

try:
    from lxml import etree
except Exception:
    sys.stderr.write('Failed to import lxml module.')


def readPercout(fname):
    doc = None
    try:
        doc = etree.parse(fname)
    except Exception:
        sys.stderr.write('Could not parse XML provided in %s or error reading file. \n' % (fname))
    return doc


def outputTabSep(fn, to_process, outputfn, ns):
    features = etree.iterparse(fn, tag='{%s}%s' % (ns['xmlns'], to_process) )
    header = {'psm':['ID', 'score', 'q-value', 'PEP', 'peptide', 'protein IDs\n'],
          'peptide':['ID', 'score', 'q-value', 'PEP', 'psm_ids', 'protein IDs\n']
          }
    with open(outputfn, 'w') as fp:
        fp.write('\t'.join(header[to_process]))
        for ac,el in features:
            p_id = el.attrib['{%s}%s_id' % (ns['xmlns'], to_process)]
            score = el.xpath('xmlns:svm_score', namespaces=ns)[0].text
            q = el.xpath('xmlns:q_value', namespaces=ns)[0].text
            PEP = el.xpath('xmlns:pep', namespaces=ns)[0].text
            proteins = '; '.join( [x.text for x in el.xpath('xmlns:protein_id',
                                    namespaces=ns) ] )

            if to_process == 'peptide':
                psm_ids = '; '.join([x.text for x in el.xpath('xmlns:psm_id',
                                    namespaces=ns) ])
                fp.write('\t'.join([p_id, score, q, PEP, psm_ids, proteins]))

            elif to_process == 'psm':
                pepseq = el.xpath('xmlns:peptide_seq',
                        namespaces=ns)[0].attrib['seq']
                fp.write('\t'.join([p_id, score, q, PEP, pepseq, proteins]))
            fp.write('\n')

def target_decoy_generator(element_generator, decoy, ns):
    for ev,el in element_generator:
        if el.attrib['{%s}decoy' % ns['xmlns']] == decoy:
            strxml = etree.tostring(el)
            strxml = strxml.replace('xmlns="{0}" '.format(ns['xmlns']), '')
            strxml = strxml.replace('xmlns:p="{0}" '.format(ns['xmlns:p']), '')
            strxml = strxml.replace('xmlns:xsi="{0}" '.format(ns['xmlns:xsi']), '')
            strxml = Markup(strxml)
        else:
            continue
        el.clear()
        while el.getprevious() is not None:
            del el.getparent()[0]
        yield strxml

def splitTargetDecoy(fn, ns):
    """ Splits XML into target/decoy/notspecified elements.
        Usage: splitTargetDecoy('test.xml', namespace)
    """
    split_elements = {'target': {}, 'decoy': {}}
    feats_to_process = ['psm', 'peptide']
    for feat in feats_to_process:
        targetfeatgen = etree.iterparse(fn, tag='{%s}%s' % (ns['xmlns'], feat))
        decoyfeatgen = etree.iterparse(fn, tag='{%s}%s' % (ns['xmlns'], feat))
        split_elements['target'][feat] = target_decoy_generator(targetfeatgen, 'false', ns)
        split_elements['decoy'][feat] = target_decoy_generator(decoyfeatgen, 'true', ns)
        
    return split_elements

def get_percolator_static_xml(fn, ns):
    rootgen = etree.iterparse(fn, tag='{%s}percolator_output' % ns['xmlns'], events=('start',))
    root = rootgen.next()[1]
    for child in root.getchildren():
        root.remove(child)
    process = etree.iterparse(fn, tag='{%s}process_info' % ns['xmlns'], events=('start',))
    root.append(process.next()[1])
    return root

def write_percolator_xml(fn, staticxml, feats, ns):
    """Given the static percolator xml root and process info nodes, and all
    psms and peptides as iterators in a dict, this generates percolator out data
    into a Genshi template."""

    # First prepare xml template for Genshi containing static percolator xml 
    # and a for loop for Genshi to fill in
    etree.register_namespace('py', 'http://genshi.edgewall.org/')
    peptides = etree.SubElement(staticxml, 'peptides', nsmap={'py': 'http://genshi.edgewall.org/'})
    psms = etree.SubElement(staticxml, 'psms', nsmap={'py': 'http://genshi.edgewall.org/'})
    etree.SubElement(peptides, '{http://genshi.edgewall.org/}for', each='pep in peps').text = '${pep}'
    etree.SubElement(psms, '{http://genshi.edgewall.org/}for', each='psm in psms').text = '${psm}'
    root = etree.tostring(staticxml, pretty_print=True, xml_declaration=True)
    
    # then insert peptides/psms into new formed template
    tpl = MarkupTemplate(root, encoding='utf-8')
    stream = tpl.generate(peps=feats['peptide'], psms=feats['psm'])
    with open(fn, 'w') as fp:
        for chunk in stream.serialize():
            fp.write(chunk)
    del(stream)
        

def refillTree(fn, ns, elements, tdb, iterfind_used=False):
    """elements: {'psm': {'target': generator, 'decoy':generator, 'both':generator}, etc } tdb='target' or 'decoy' or 'both' """
    # generate percolator_output root element with process_info element from existing file fn
    topgen = etree.iterparse(fn, tag='{%s}percolator_output' % ns['xmlns'])
    for ac, outtree in topgen:
        outtree.clear()
    process_gen = etree.iterparse(fn, tag='{%s}process_info' % ns['xmlns'])
    for ac, proc in process_gen:
        outtree.append(proc)
    
    # loop through generators and attach to tree
    for feat in elements:
        top_element = etree.Element(feat+'s')
        outtree.append(top_element)
        if iterfind_used:
            for el in elements[feat][tdb]:
                top_element.append(el)
        else:
            for ac,el in elements[feat][tdb]:
                top_element.append(el)
    
    return etree.ElementTree(outtree)



def filterUniques(input_file, files_dir, score, ns):
    """ Filters unique peptides from (multiple) Percolator output XML files.
        Takes a dir with a set of XML files, a score to filter on and a namespace.
        Outputs an ElementTree.
    """
    # lookup dict/initializing
    scores = {'q':'q_value', 'pep':'pep', 'p':'p_value', 'svm':'svm_score'} 
    filtered_peptides = {}
    try:
        members = [os.path.join(files_dir, fn) for fn in os.listdir(files_dir)]
    except OSError: # no extra_files_dir
        members = [input_file]
    if len(members) == 1:
        print >> sys.stdout, 'WARNING: You are attempting to filter unique peptides from a non-fractioned dataset. This may be unneccessary.'
    
    # parse XML
    peptide_generators_fractions = []
    for fn in members:
        peptide_generators_fractions.append(etree.iterparse(fn, tag='{%s}peptide' % ns['xmlns']))
    
    # get results
    for peptide_generator in peptide_generators_fractions:
        for ac, peptide_el in peptide_generator:
            # It's actually faster to loop through the feat's children, 
            # but this is 2-line code and still readable.
            featscore = float(peptide_el.xpath('xmlns:%s' % scores[score], namespaces=ns)[0].text)
            seq = str(peptide_el.xpath('@xmlns:peptide_id', namespaces=ns)[0])
            
            if seq not in filtered_peptides:
                filtered_peptides[seq] = peptide_el
            elif score == 'svm': # greater than score is accepted
                if featscore > filtered_peptides[seq]:
                    filtered_peptides[seq] = peptide_el
            else: # lower than score is accepted
                if featscore < filtered_peptides[seq]:
                    filtered_peptides[seq] = peptide_el
        
    # create top element and append the found unique high scorers
    top_el = etree.iterparse(fn, tag='{%s}peptides' % ns['xmlns'] ).next()[1]
    top_el.clear()
    
    for seq in filtered_peptides:
        top_el.append(filtered_peptides[seq])
    # convert to generator
    filtered_peptides = {'peptide': {} }
    filtered_peptides['peptide']['both'] = etree.iterwalk(top_el, tag='{%s}peptide' % ns['xmlns'] )
    
    # output full percolator tree
    return refillTree(fn, ns, filtered_peptides, 'both')


def mergeXML(datasets, ns):
    """ Merges multiple Percolator XML output documents into a single doc.
    Useful for example when a target and decoy set have been created that have to be
    remerged.
    """
    outdoc = readPercout(datasets[0])
    root = outdoc.getroot()
    
    for parentfeat in ['psms', 'peptides', 'proteins']:
        featname = parentfeat[:-1]
        parent_el = root.xpath('//xmlns:%s' % parentfeat, namespaces=ns)
        print >>sys.stdout, len(parent_el)
        if parent_el: # empty lists are skipped
            parent_el = parent_el[0]
        for fn in datasets[1:]:
            elements = etree.iterparse(fn, tag='{%s}%s' % (ns['xmlns'], featname) )
            for ac,el in elements:
                parent_el.append(el)
        print >>sys.stdout, len(parent_el)
    
    return etree.ElementTree(root)


def untar(tar):
    try:
        with tarfile.open(tar, 'r') as f:
            members = f.getmembers()
            f.extractall()
            return [x.name for x in members]
    except:
        sys.stdout.write('Could not extract Percolator files from dataset: %s \n' % tar)
        return 1

def writeXML(*args):
    """Takes a filename and _ElementTree arguments in tuples ((fn, tree), (fn2,tree2),...) 
    and writes an XML file.
    """
    try:
    # hello! shouldn't these be exceptions instead?
        assert False not in [isinstance(x, tuple) for x in args] # arguments are tuples
        assert False not in [len(x)>2 for x in args] # tuples are 2-length
        assert False not in [isinstance(x[0], str) for x in args] # first arg from tuple is string for fn.
        assert False not in [isinstance(x[1], etree._ElementTree) for x in args] # second arg from tuple is ElementTree.
    except:
        Exception('function writeXML takes arguments in form of tuples (filename, XMLtree)')
        
    for fn, doc in args:
        try:
            outfile = open(fn, 'w')
        except Exception:
            sys.stderr('Unable to write XML output to file %s \n' % fn)
        doc.write('<?xml version="1.0" encoding="UTF-8"?>')
        doc.write(outfile)

def parseOptions(args):
    ns = {'xmlns':'http://per-colator.com/percolator_out/14',
    'xmlns:p' : 'http://per-colator.com/percolator_out/14',
    'xmlns:xsi':'http://www.w3.org/2001/XMLSchema-instance',
    'xsi:schemaLocation':'http://per-colator.com/percolator_out/14,https://github.com/percolator/percolator/raw/pout-1-4/src/xml/percolator_out.xsd',
    'p:majorVersion':'0',
    'p:minorVersion':'00',
    'p:percolator_version':'Percolator version 0.00'
    }

    if args[0] == 'splittd':
        fname = args[1]
        staticxml = get_percolator_static_xml(fname, ns)
        split_els = splitTargetDecoy(fname, ns)
        write_percolator_xml(args[2], staticxml, split_els['target'], ns)
        write_percolator_xml(args[3], staticxml, split_els['decoy'], ns)
    
    elif args[0] == 'filter_uni':
        uniques = filterUniques(args[1], args[2], args[3], ns)
        writeXML((args[4], uniques))
    
    elif args[0] == 'merge':
        with open(args[1]) as fp:
            datasets = fp.read()
        datasets = [opt.strip() for opt in datasets.strip().split('\n')]
        merged = mergeXML(datasets, ns)
        writeXML((args[2], merged))
    
    elif args[0] == 'untar_merge':
        files = os.listdir(args[1])
        merged = mergeXML(files, ns)
        writeXML((args[2], merged))

    elif args[0] == 'tab':
        outputTabSep(args[1], args[2], args[3], ns)
        
    else:
        sys.stderr.write('Argument %s not recognized.\n' % args[0])
        return 1


def main():
    parseOptions(sys.argv[1:])


if __name__ == '__main__':
    main()
