from lxml import etree


def write_percolator_xml(staticxml, feats, fn):
    """Given the static percolator xml root and process info nodes, and all
    psms and peptides as iterators in a dict {'peptide': pep_iterator, 'psm':
    psm_iterator}, this generates percolator out data into a file."""

    # First get xml until psms opening element is found.
    etree.SubElement(staticxml, 'psms').text = '***psms***'
    root = etree.tostring(staticxml, pretty_print=True,
                          xml_declaration=True, encoding='UTF-8')
    root = root.decode('utf-8')
    root = root[:root.find('***psms***')]

    # Write opening xml
    with open(fn, 'w') as fp:
        fp.write(root)
        fp.write('\n')

    # Then write features
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
    print('Wrote {0} psms, {1} peptides to file {2}'.format(psmcount,
                                                            peptidecount, fn))


def write_qvality_input(scores, fn):
    with open(fn, 'w') as fp:
        for score in scores:
            fp.write('{0}\n'.format(score))
