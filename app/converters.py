def peptide_to_tsv(features, ns):
    for feature in features:
        line = percolator_feature_to_tsv(feature, ns, 'peptide')
        line['psm_ids'] = '; '.join([x.text for x in feature.xpath('xmlns:psm_id',
                                                                   namespaces=ns)])
        yield line


def psm_to_tsv(features, ns):
    for feature in features:
        line = percolator_feature_to_tsv(feature, ns, 'psm')
        line['pepseq'] = feature.xpath(
            'xmlns:peptide_seq', namespaces=ns)[0].attrib['seq']
        yield line


def percolator_feature_to_tsv(feature, ns, feattype):
    line = {'p_id': feature.attrib['{%s}%s_id' % (ns['xmlns'], feattype)],
            'score': feature.xpath('xmlns:svm_score', namespaces=ns)[0].text,
            'q': feature.xpath('xmlns:q_value', namespaces=ns)[0].text,
            'PEP': feature.xpath('xmlns:pep', namespaces=ns)[0].text,
            'proteins': '; '.join([x.text for x in feature.xpath(
                'xmlns:protein_id', namespaces=ns)]),
            }
    return line
