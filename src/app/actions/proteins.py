from math import log

from app.dataformats import peptable as peptabledata
from app.dataformats import prottable as prottabledata
from app.dataformats import mzidtsv as mzidtsvdata

from app.actions.psmtopeptable import evaluate_peptide
from app.readers import fasta


def generate_bestpep_proteins(peptides, scorecol, minlog, outputaccfield, 
        protcol, higherbetter=True):
    """Best peptide for each protein in a table"""
    protein_peptides = {}
    if minlog:
        higherbetter = False
    if not protcol:
        protcol = peptabledata.HEADER_MASTERPROTEINS
    for peptide in peptides:
        p_acc = peptide[protcol]
        if ';' in p_acc or p_acc == 'NA':
            continue
        protein_peptides = evaluate_peptide(protein_peptides, peptide, p_acc,
                higherbetter, scorecol, fncol=False)
    if minlog:
        try:
            nextbestscore = min([pep['score'] for pep in
                                 protein_peptides.values()
                                 if pep['score'] > 0])
        except ValueError:
            print('WARNING: Cannot find score of type {} which is above 0. '
                    'Only scores above zero can have a -log value.'.format(scorecol))
            nextbestscore = -log(1e-06, 10) # 1 decoy in a million, fake value
        else:
            nextbestscore = -log(nextbestscore, 10)
    for protein, peptide in protein_peptides.items():
        if minlog and peptide['score'] != 'NA':
            peptide['score'] = log_score(peptide['score'], nextbestscore)
        yield {
                outputaccfield: protein, 
                prottabledata.HEADER_QSCORE: str(peptide['score'])
                }


def log_score(score, nextbestscore):
    try:
        logged_score = -log(score, 10)
    except ValueError:
        logged_score = nextbestscore + 2
    return logged_score


def generate_protein_fdr(target, decoy, featfield):
    tproteins = [x for x in target if get_score(x) is not None]
    dproteins = [x for x in decoy if get_score(x) is not None]
    [x.update({'target_decoy': 'target'}) for x in tproteins]
    [x.update({'target_decoy': 'decoy'}) for x in dproteins]
    proteins = sorted(tproteins + dproteins, key=lambda x: get_score(x),
                      reverse=True)
    #fdrheader = headerfields['proteinfdr'][prottabledata.HEADER_QVAL][None]
    # FIXME complex stupid headerfields system, remove it?
    fdrheader = prottabledata.HEADER_QVAL
    for protein in qvalue_generator(fdrheader, proteins, featfield):
        yield protein


def generate_pick_fdr(target, decoy, tfastafn, dfastafn, picktype, featfield,
        fastadelim, genefield):
    t_scores, d_scores = {}, {}
    for protein in target:
        acc = protein[featfield]
        t_scores[acc] = protein
        t_scores[acc]['target_decoy'] = 'target'
    for protein in decoy:
        acc = protein[featfield]
        d_scores[acc] = protein
        d_scores[acc]['target_decoy'] = 'decoy'
    if picktype == 'fasta':
        prefixlen = len(mzidtsvdata.DECOY_PREFIX)
        genetype = 'genename' if featfield == mzidtsvdata.HEADER_SYMBOL else 'ensg'
        tfasta = fasta.get_genes_pickfdr(tfastafn, genetype, fastadelim, genefield)
        dfasta = fasta.get_genes_pickfdr(dfastafn, genetype, fastadelim, genefield)
        tdmap = {}
        for target, decoy in zip(tfasta, dfasta):
            tdmap[target] = decoy
    elif picktype == 'result':
        # FIXME this code path is not accessed but kept  in case we ever need
        # to get back to result picktype. See commit 58080fdc11d0800117a09122d66173e7ac54292a
        tdmap = {}
        for dprot in d_scores:
            faketprot = dprot.replace(mzidtsvdata.DECOY_PREFIX, '')
            tdmap[faketprot] = dprot
        for tprot in t_scores:
            if not tdmap.get(tprot, False):
                tdmap[tprot] = None
    picked_proteins = []
    for tgene, dgene in tdmap.items():
        picked = pick_target_decoy(t_scores.get(tgene), d_scores.get(dgene))
        if picked:
            picked_proteins.append(picked)
    sorted_proteins = sorted(picked_proteins, key=lambda x: get_score(x),
                             reverse=True)
    #fdrheader = headerfields['proteinfdr'][prottabledata.HEADER_QVAL][None]
    fdrheader = prottabledata.HEADER_QVAL
    for protein in qvalue_generator(fdrheader, sorted_proteins, featfield):
        yield protein


def qvalue_generator(fdrheader, sorted_features, featfield):
    tdcounter = {'target': 0, 'decoy': 0}
    previousscore = get_score(sorted_features[0])
    outfeats = []
    for feat in sorted_features:
        outfeat = {k: v for k, v in feat.items()}
        score = get_score(outfeat)
        if score != previousscore:
            # new score, all proteins with previous score get same fdr
            previousscore = score
            try:
                fdr = tdcounter['decoy'] / float(tdcounter['target'])
            except ZeroDivisionError:
                fdr = 1
            for feat in sorted(outfeats, key=lambda x: x[featfield]):
                feat[fdrheader] = fdr
                yield feat
            outfeats = []
        tdcounter[outfeat['target_decoy']] += 1
        # Only report target hits so FDR=D/T
        if outfeat['target_decoy'] == 'target':
            outfeats.append(outfeat)
    # All proteins from bottom of list (no new score) get FDR as well
    try:
        fdr = tdcounter['decoy'] / float(tdcounter['target'])
    except ZeroDivisionError:
        fdr = 1
    for feat in sorted(outfeats, key=lambda x: x[featfield]):
        feat[fdrheader] = fdr
        yield feat


def get_score(protein):
    try:
        return float(protein[prottabledata.HEADER_QSCORE])
    except (TypeError, ValueError):
        return None


def pick_target_decoy(tfeature, dfeature):
    """Feed it with a target and decoy score and the protein/gene/id names,
    and this will return target/decoy type, the winning ID and the score"""
    tscore, dscore = get_score(tfeature), get_score(dfeature)
    if tscore is None and dscore is None:
        return False
    elif None in [tscore, dscore]:
        # return the non-False feature
        return [v for k, v in {tscore: tfeature, dscore: dfeature}.items()
                if k is not False][0]
    elif tscore == dscore:
        # same score or both False
        return False
    elif tscore > dscore:
        return tfeature
    elif tscore < dscore:
        return dfeature
    else:
        # in case uncaught edgecase occurs
        print('WARNING, target score {} and decoy score {} could not be '
              'compared'.format(tscore, dscore))
        return False


def generate_top_ms1_peptides(peptides, protcol):
    """Fed with a peptides generator, this returns the 3 PSMs with
    the highest precursor intensities (or areas, or whatever is
    given in the HEADER_PRECURSOR_QUANT"""
    top_ms1_peps = {}
    for peptide in peptides:
        protacc = peptide[protcol]
        precursor_amount = peptide[peptabledata.HEADER_AREA]
        if ';' in protacc or precursor_amount == 'NA':
            continue
        precursor_amount = float(precursor_amount)
        pepseq = peptide[peptabledata.HEADER_PEPTIDE]
        try:
            peptide_area = top_ms1_peps[protacc][pepseq]
        except KeyError:
            try:
                top_ms1_peps[protacc][pepseq] = precursor_amount
            except KeyError:
                top_ms1_peps[protacc] = {pepseq: precursor_amount}
        else:
            if precursor_amount > peptide_area:
                top_ms1_peps[protacc][pepseq] = precursor_amount
    return top_ms1_peps


def calculate_protein_precursor_quant(top_ms1_peps, prot_acc):
    try:
        amounts = top_ms1_peps[prot_acc].values()
    except KeyError:
        return 'NA'
    else:
        amounts = sorted([x for x in amounts if x > 0], reverse=True)[:3]
        return sum(amounts) / len(amounts)


def add_ms1_quant_from_top3_mzidtsv(features, peptides, outputaccfield, featcol):
    """Collects peptides with the highest precursor quant values,
    adds sum of the top 3 of these to a protein table"""
    top_ms1_peps = generate_top_ms1_peptides(peptides, featcol)
    for feat in features:
        acc = feat[outputaccfield]
        prec_area = calculate_protein_precursor_quant(top_ms1_peps, acc)
        outfeat = {k: v for k, v in feat.items()}
        outfeat[prottabledata.HEADER_AREA] = str(prec_area)
#        outprotein[headerfields['precursorquant'][
#            prottabledata.HEADER_AREA][None]] = str(prec_area)
        yield outfeat 
