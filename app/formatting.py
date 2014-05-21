from lxml import etree

def string_and_clear(el, ns):
    str_el = stringify_strip_namespace_declaration(el, ns)
    clear_el(el)
    return str_el

def stringify_strip_namespace_declaration(el, ns):
    strxml = etree.tostring(el, encoding='utf-8').decode('utf-8')
    strxml = strxml.replace('xmlns="{0}" '.format(ns['xmlns']), '')
    strxml = strxml.replace('xmlns:p="{0}" '.format(ns['xmlns:p']), '')
    strxml = strxml.replace('xmlns:xsi="{0}" '.format(ns['xmlns:xsi']), '')
    return strxml

def clear_el(el):
    el.clear()
    while el.getprevious() is not None:
        del(el.getparent()[0])

