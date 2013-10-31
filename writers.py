from lxml import etree
from genshi.template import MarkupTemplate


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


