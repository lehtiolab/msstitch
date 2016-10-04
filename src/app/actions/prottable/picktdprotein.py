from app.readers import tsv as reader
from app.readers import fasta
from app.dataformats import prottable as prottabledata
from app.dataformats import mzidtsv as mzidtsvdata


TARGET = 't'
DECOY = 'd'


def write_pick_td_tables(target, decoy, theader, dheader,
                         targetfasta, decoyfasta, inferencetype,
                         fastadelim, genefield):
    tfile, dfile = 'target.txt', 'decoy.txt'
    tdmap = {}
    with open(tfile, 'w') as tdmap[TARGET], open(dfile, 'w') as tdmap[DECOY]:
        for ptype, score in generate_pick_fdr(target, decoy, theader, dheader,
                                              targetfasta, decoyfasta,
                                              inferencetype,
                                              fastadelim, genefield):
            tdmap[ptype].write('\n{}'.format(score))
    return tfile, dfile


def generate_pick_fdr(targetprot, decoyprot, theader, dheader,
                      targetfasta, decoyfasta, picktype,
                      fastadelim, genefield):
    t_scores, d_scores = generate_scoremaps(targetprot, decoyprot,
                                            theader, dheader)
    if picktype == 'fasta':
        tdmap = create_td_gene_map(targetfasta, decoyfasta,
                                   fastadelim, genefield)
    elif picktype == 'result':
        tdmap = create_td_assoc_map(targetprot, decoyprot, theader, dheader)
    td_picked = generate_picked_genes(tdmap, t_scores, d_scores)
    for ptype, score in td_picked:
        yield ptype, score


def generate_picked_genes(tdmap, tscores, dscores):
    for tgene, dgene in tdmap.items():
        picked = pick_target_decoy(tscores.get(tgene), dscores.get(dgene))
        if picked:
            yield picked


def get_score(protein):
    return protein[prottabledata.HEADER_QSCORE]


def create_td_gene_map(tfastafn, dfastafn, fastadelim, genefield):
    tfasta = (x[1].split('.')[0] for x in
              fasta.get_proteins_genes(tfastafn, fastadelim, genefield))
    dfasta = (x[1].split('.')[0] for x in
              fasta.get_proteins_genes(dfastafn, fastadelim, genefield))
    tdmap = {}
    for target, decoy in zip(tfasta, dfasta):
        tdmap[target] = decoy
    return tdmap


def create_td_assoc_map(tproteins, dproteins, theader, dheader):
    tdmap = {}
    for dprot in reader.generate_tsv_proteins(dproteins, dheader):
        dprot = dprot[prottabledata.HEADER_PROTEIN]
        faketprot = dprot.replace(mzidtsvdata.DECOY_PREFIX, '')
        tdmap[faketprot] = dprot
    for tprot in reader.generate_tsv_proteins(tproteins, theader):
        tprot = tprot[prottabledata.HEADER_PROTEIN]
        if not tdmap.get(tprot, False):
            tdmap[tprot] = None
    return tdmap


def generate_scoremaps(targetprot, decoyprot, theader, dheader):
    t_scores, d_scores = {}, {}
    for protein in reader.generate_tsv_proteins(targetprot, theader):
        t_scores[protein[prottabledata.HEADER_PROTEIN]] = get_score(protein)
    for protein in reader.generate_tsv_proteins(decoyprot, dheader):
        d_scores[protein[prottabledata.HEADER_PROTEIN]] = get_score(protein)
    return t_scores, d_scores


def pick_target_decoy(tscore, dscore):
    """Feed it with a target and decoy score and the protein/gene/id names,
    and this will return target/decoy type, the winning ID and the score"""
    try:
        tscore = float(tscore)
    except (ValueError, TypeError):
        tscore = False
    try:
        dscore = float(dscore)
    except (ValueError, TypeError):
        dscore = False
    falsecheck = {score: (ptype, score) for score, ptype
                  in zip([tscore, dscore], [TARGET, DECOY])}
    if len(falsecheck) == 1:
        return False
    elif False in falsecheck:
        falsecheck.pop(False)
        return next(iter(falsecheck.values()))
    elif tscore > dscore:
        return TARGET, tscore
    elif tscore < dscore:
        return DECOY, dscore
    else:
        # in case uncaught edgecase occurs
        print('WARNING, target score {} and decoy score {} could not be '
              'compared'.format(tscore, dscore))
        return False
