from app.readers import tsv as reader
from app.readers import fasta
from app.dataformats import prottable as prottabledata
from app.dataformats import mzidtsv as mzidtsvdata


def generate_protein_fdr(target, decoy, theader, dheader, headerfields):
    tproteins = [x for x in reader.generate_tsv_proteins(target, theader)]
    dproteins = [x for x in reader.generate_tsv_proteins(decoy, dheader)]
    [x.update({'target_decoy': 'target'}) for x in tproteins]
    [x.update({'target_decoy': 'decoy'}) for x in dproteins]
    proteins = sorted(tproteins + dproteins, key=lambda x: get_score(x),
                      reverse=True)
    fdrheader = headerfields['proteinfdr'][prottabledata.HEADER_QVAL][None]
    for protein in qvalue_generator(fdrheader, proteins):
        yield protein


def generate_pick_fdr(target, decoy, theader, dheader, headerfields,
                      targetfasta, decoyfasta, inferencetype,
                      fastadelim, genefield):
    picked_proteins = get_picked_proteins(target, decoy, theader, dheader,
                                          targetfasta, decoyfasta,
                                          inferencetype, fastadelim, genefield)
    sorted_proteins = sorted(picked_proteins, key=lambda x: get_score(x),
                             reverse=True)
    fdrheader = headerfields['proteinfdr'][prottabledata.HEADER_QVAL][None]
    for protein in qvalue_generator(fdrheader, sorted_proteins):
        yield protein


def qvalue_generator(fdrheader, sorted_features):
    tdcounter = {'target': 0, 'decoy': 0}
    for feat in sorted_features:
        outfeat = {k: v for k, v in feat.items()}
        tdcounter[feat['target_decoy']] += 1
        fdr = tdcounter['decoy'] / float(tdcounter['decoy'] +
                                         tdcounter['target'])
        if outfeat['target_decoy'] == 'target':
            outfeat[fdrheader] = fdr
            yield outfeat


def get_picked_proteins(targetprot, decoyprot, theader, dheader, targetfasta,
                        decoyfasta, picktype, fastadelim, genefield):
    t_feats, d_feats = generate_accessionmaps(targetprot, decoyprot, theader,
                                              dheader)
    if picktype == 'fasta':
        tdmap = create_td_gene_map(targetfasta, decoyfasta,
                                   fastadelim, genefield)
    elif picktype == 'result':
        tdmap = create_td_assoc_map(targetprot, decoyprot, theader, dheader)
    picked_proteins = []
    for tgene, dgene in tdmap.items():
        picked = pick_target_decoy(t_feats.get(tgene), d_feats.get(dgene))
        if picked:
            picked_proteins.append(picked)
    return picked_proteins


def get_score(protein):
    try:
        return float(protein[prottabledata.HEADER_QSCORE])
    except (TypeError, ValueError):
        return False


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


def generate_accessionmaps(targetprot, decoyprot, theader, dheader):
    t_scores, d_scores = {}, {}
    for protein in reader.generate_tsv_proteins(targetprot, theader):
        acc = protein[prottabledata.HEADER_PROTEIN]
        t_scores[acc] = protein
        t_scores[acc]['target_decoy'] = 'target'
    for protein in reader.generate_tsv_proteins(decoyprot, dheader):
        acc = protein[prottabledata.HEADER_PROTEIN]
        d_scores[acc] = protein
        d_scores[acc]['target_decoy'] = 'decoy'
    return t_scores, d_scores


def pick_target_decoy(tfeature, dfeature):
    """Feed it with a target and decoy score and the protein/gene/id names,
    and this will return target/decoy type, the winning ID and the score"""
    tscore, dscore = get_score(tfeature), get_score(dfeature)
    if tscore == dscore:
        # same score or both False
        return False
    elif False in [tscore, dscore]:
        # return the non-False feature
        return [v for k, v in {tscore: tfeature, dscore: dfeature}.items()
                if k is not False][0]
    elif tscore > dscore:
        return tfeature
    elif tscore < dscore:
        return dfeature
    else:
        # in case uncaught edgecase occurs
        print('WARNING, target score {} and decoy score {} could not be '
              'compared'.format(tscore, dscore))
        return False
