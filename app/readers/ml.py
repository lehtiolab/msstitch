def get_scan_nr(element, attribname):
    """General method to get a scan nr from xml element of mzML/mzIdentML"""
    info = element.attrib[attribname]
    infomap = {y[0]: y[1] for x in info.split() for y in x.split('=')}
    return infomap['scan']
