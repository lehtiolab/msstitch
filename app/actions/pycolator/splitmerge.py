from lxml import etree

from app.readers import xmlformatting as formatting
from app.readers import pycolator as reader


def create_merge_psm_map(peptides, ns):
    """Loops through peptides, stores sequences mapped to PSM ids."""
    psmmap = {}
    for peptide in peptides:
        seq = reader.get_peptide_seq(peptide, ns)
        psm_ids = reader.get_psm_ids_from_peptide(peptide, ns)
        for psm_id in psm_ids:
            try:
                psmmap[seq][psm_id.text] = 1
            except KeyError:
                psmmap[seq] = {psm_id.text: 2}
    for seq, psm_id_dict in psmmap.items():
        psmmap[seq] = [x for x in psm_id_dict]
    return psmmap


def merge_peptides(fns, ns):
    """Loops peptides from multiple files, fetches PSMs from
    sequence:PSM map, outputs correctly PSM mapped peptides"""
    peptides_to_map = reader.generate_peptides_multiple_fractions(fns, ns)
    psmmap = create_merge_psm_map(peptides_to_map, ns)
    peptides = reader.generate_peptides_multiple_fractions(fns, ns)
    for peptide in peptides:
        seq = reader.get_peptide_seq(peptide, ns)
        psm_ids = reader.get_psm_ids_from_peptide(peptide, ns)
        # remove current psm ids, repopulate with stored ones
        psm_ids.clear()
        for new_psm_id in psmmap[seq]:
            etree.SubElement(psm_ids, 'psm_id').text = new_psm_id
        yield formatting.string_and_clear(peptide, ns)


def target_decoy_generator(element_generator, decoy, ns):
    for el in element_generator:
        if el.attrib['{%s}decoy' % ns['xmlns']] == decoy:
            yield formatting.string_and_clear(el, ns)
        else:
            formatting.clear_el(el)


def split_target_decoy(elements, ns, filter_type):
    td = {'target': 'false', 'decoy': 'true'}
    feats_to_process = {'psm': None, 'peptide': None}
    for feat in feats_to_process:
        feats_to_process[feat] = target_decoy_generator(
            elements[feat], td[filter_type], ns)
    return feats_to_process
