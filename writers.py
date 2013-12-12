from lxml import etree

def write_percolator_xml(staticxml, feats, fn):
    """Given the static percolator xml root and process info nodes, and all
    psms and peptides as iterators in a dict {'peptide': pep_iterator, 'psm':
    psm_iterator}, this generates percolator out data into a file."""

    etree.SubElement(staticxml, 'psms').text='{$psms}'
    root = etree.tostring(staticxml, pretty_print=True, xml_declaration=True)
    root = root[:root.find('{$psms}')]
    with open(fn, 'w') as fp:
        fp.write(root)
        fp.write('\n')
    with open(fn, 'a') as fp:
        psmcount = 0
        for psm in feats['psm']:
            psmcount += 1
            fp.write(psm)
            fp.write('\n')
        fp.write('</psms><peptides>\n')

        peptidecount = 0
        for pep in feats['peptide']:
            peptidecount += 1
            fp.write(pep)
            fp.write('\n')
        fp.write('</peptides></percolator_output>')
    print 'Wrote {0} psms, {1} peptides to file {2}'.format(psmcount,
                                            peptidecount, fn)


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


