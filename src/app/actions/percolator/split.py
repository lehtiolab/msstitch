from lxml import etree
import re

from app.readers import xmlformatting as formatting
from app.readers import percolator as reader


def protein_header_split_generator(elements, ns, can_headers, headers):
    """Loop through proteins of each PSM/peptide. If a protein does not
    match any of headers, discard PSM/peptide immediately"""
    for el in elements:
        header_matching = False
        can = False
        for protein in el.findall('{%s}protein_id' % ns['xmlns']):
            if any(re.search(h, protein.text) for h in can_headers):
                can = True
                break #as soon as a canonical match was found break
            """for classes other than known,
               check if there is at least one protein matching the specified header
               and those with matches to the canonical proteins will not be used"""
            if any(re.search(h, protein.text) for h in headers):
                    header_matching = True
        if (header_matching and not can) or ((headers == can_headers) and can):
            yield formatting.string_and_clear(el, ns)
        else:
            formatting.clear_el(el)


def split_protein_header_id_type(elements, ns, protheaders):
    can_headers = protheaders.split('|')[0].split(':')[-1].strip(';').split(';')
    other_headers = protheaders.split('|')[-1].split(':')[-1].strip(';').split(';')
    return {x: protein_header_split_generator(elements[x], ns, can_headers, other_headers)
            for x in ['psm', 'peptide']}

