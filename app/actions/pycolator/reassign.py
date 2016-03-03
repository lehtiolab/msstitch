import sys

from app.readers import xmlformatting as formatting


def parse_qvality_output(fn):
    statistics = {}
    with open(fn) as fp:
        for line in fp:
            line = line.strip()
            line = line.split('\t')
            if len(line) == 3:
                try:
                    score = float(line[0])
                except ValueError:
                    continue
                statistics[score] = {'PEP': line[1], 'q': line[2]}
    return statistics


def reassign_elements(elements, stats, ns):
    for el in elements:
        score = round(float(el.xpath('xmlns:svm_score',
                                     namespaces=ns)[0].text), 5)
        oldq = el.xpath('xmlns:q_value', namespaces=ns)[0]
        oldpep = el.xpath('xmlns:pep', namespaces=ns)[0]
        newq, newpep, warning = lookup_statistic(score, stats)
        if warning is not None:
            sys.stdout.write(warning)
        oldq.text, oldpep.text = newq, newpep
        yield formatting.string_and_clear(el, ns)


def lookup_statistic(score, stats):
    """ Finds statistics that correspond to PSM/peptide/protein feature's
    score. Loops through list of qvality generated scores until it finds
    values closest to the feature's svm_score."""
    if score in stats:
        return stats[score]['q'], stats[score]['PEP'], None

    else:
        lower, warning = None, None
        for stat_score in sorted(stats.keys()):
            if score < stat_score:
                break
            lower = stat_score
        if score > stat_score:
            warning = ('WARNING! Values found with higher svm_score than in '
                       'qvality recalculation!')
        if not lower:
            warning = ('WARNING! Values found with lower svm_score than in '
                       'qvality recalculation were set to PEP=1, qval=1.')
            return '1', '1', warning
        qval = (float(stats[stat_score]['q']) + float(stats[lower]['q'])) / 2
        pep = (float(stats[stat_score]['PEP']) + float(stats[lower]['PEP']))/2

    return str(qval), str(pep), warning


def fix_psmid_with_incomplete_filename(psm_id, fn):
    """It looks like that filenames are cut by msgf2pin for inclusion
    in the psm id. Either by removing extension (may contain fraction or task),
    or by truncating a certain length."""
    if psm_id.startswith(fn):
        return psm_id

    charindex = 0
    while True:
        charindex += 1
        # cut incomplete filename from psm id
        if psm_id[:charindex] != fn[:charindex]:
            new_psm_id = psm_id[charindex - 1:]
            break
    # add underscore if needed, then filename
    if not new_psm_id[0] == '_':
        new_psm_id = '_{0}'.format(new_psm_id)
    new_psm_id = '{0}{1}'.format(fn, new_psm_id)
    return new_psm_id
