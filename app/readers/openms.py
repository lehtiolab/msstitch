from . import basereader


def quant_generator(consfiles):
    return basereader.generate_tags_multiple_files(
        consfiles,
        'consensusElement',
        ['consensusElementList'],
    )


def get_consxml_rt(cons_el):
    """Returns consensusXML in minutes"""
    rt = cons_el.find('centroid').attrib['rt']
    return float(rt)
