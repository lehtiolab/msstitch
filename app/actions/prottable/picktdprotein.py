from app.readers import tsv as reader
from app.readers import fasta
from app.dataformats import prottable as prottabledata


TARGET = 't'
DECOY = 'd'


def write_pick_td_tables(target, decoy, theader, dheader,
                         targetfasta, decoyfasta):
    tfile, dfile = 'target.txt', 'decoy.txt'
    tdmap = {}
    header = '{}\t{}'.format(prottabledata.HEADER_PROTEIN,
                             prottabledata.HEADER_QSCORE)
    with open(tfile, 'w') as tdmap[TARGET], open(dfile, 'w') as tdmap[DECOY]:
        tdmap[TARGET].write(header)
        tdmap[DECOY].write(header)
        for pdata in generate_pick_fdr(target, decoy, theader, dheader,
                                       targetfasta, decoyfasta):
            ptype, protein, score = pdata
            tdmap[ptype].write('\n{}\t{}'.format(protein, score))
    return tfile, dfile


def generate_pick_fdr(targetprot, decoyprot, theader, dheader, 
                      targetfasta, decoyfasta):
    t_scores, d_scores = generate_scoremaps(targetprot, decoyprot, 
                                            theader, dheader)
    tfasta = fasta.generate_proteins_id(targetfasta)
    dfasta = fasta.generate_proteins_id(decoyfasta)
    for target, decoy in zip(tfasta, dfasta):
        picked = pick_target_decoy(t_scores, d_scores, target, decoy)
        if picked:
            ptype, pickedprotein, score = picked
            yield ptype, pickedprotein, score


def get_score(protein):
    return protein[prottabledata.HEADER_QSCORE]


def generate_scoremaps(targetprot, decoyprot, theader, dheader):
    t_scores, d_scores = {}, {}
    for protein in reader.generate_tsv_proteins(targetprot, theader):
        t_scores[protein[prottabledata.HEADER_PROTEIN]] = get_score(protein)
    for protein in reader.generate_tsv_proteins(decoyprot, dheader):
        d_scores[protein[prottabledata.HEADER_PROTEIN]] = get_score(protein)
    return t_scores, d_scores


def pick_target_decoy(targetscores, decoyscores, target, decoy):
    # target not, decoy not, both target and decoy, target or decoy
    try:
        tscore = float(targetscores[target])
    except KeyError:
        tscore = False
    try:
        dscore = float(decoyscores[decoy])
    except KeyError:
        dscore = False
    if tscore > dscore:
        return TARGET, target, tscore
    elif tscore < dscore:
        return DECOY, decoy, dscore
    else:
        return False
