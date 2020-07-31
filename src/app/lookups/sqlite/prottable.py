from app.lookups.sqlite.protpeptable import ProtPepTable
from app.dataformats import prottable as ph


class ProtGeneTableBase(ProtPepTable):
    stdheaderfields = [
            ph.HEADER_NO_PSM,
            ph.HEADER_NO_PEPTIDE,
            ph.HEADER_NO_UNIPEP,
            ]
    singlefields = stdheaderfields + [ph.HEADER_QVAL, ph.HEADER_AREA, ph.HEADER_NO_FULLQ_PSMS]
    

class ProtTableDB(ProtGeneTableBase):
    datatype = 'protein'
    colmap = {'protein_group_master': ['master_id', 'protein_acc'],
              'proteins': ['pacc_id', 'protein_acc'],
              'protein_precur_quanted': ['pacc_id', 'prottable_id', 'quant'],
              'protein_fdr': ['pacc_id', 'prottable_id', 'fdr'],
              'protein_iso_fullpsms': ['pacc_id', 'prottable_id', 'amount_psms'],
              'protquant_channels': ['channel_id', 'prottable_id',
                                     'channel_name', 'amount_psms_name'],
              'protein_iso_quanted': ['proteinquant_id', 'pacc_id',
                                      'channel_id', 'quantvalue',
                                      'amount_psms'],
              }

    def create_pdata_map(self):
        """This runs only once, returns the data which is not dependent on sets,
        in a dict with accessions as keys"""
        sql = """
SELECT pgm.master_id, p.protein_acc, IFNULL(g.gene_acc, 'NA'),
    IFNULL(aid.assoc_id, 'NA'), cov.coverage, sub.pgc, 
    sub.pgcnr, IFNULL(pd.description, 'NA') FROM protein_group_master AS pgm
    LEFT OUTER JOIN (
            SELECT master_id, GROUP_CONCAT(protein_acc) AS pgc, 
            COUNT(protein_acc) AS pgcnr FROM protein_group_content 
            GROUP BY master_id
            ) AS sub ON sub.master_id=pgm.master_id 
    INNER JOIN proteins AS p ON pgm.pacc_id=p.pacc_id
    LEFT OUTER JOIN protein_coverage AS cov ON 
    p.protein_acc=cov.protein_acc 
    LEFT OUTER JOIN ensg_proteins AS egp ON pgm.pacc_id=egp.pacc_id
    LEFT OUTER JOIN genes AS g ON egp.gene_id=g.gene_id
    LEFT OUTER JOIN genename_proteins AS gnp ON pgm.pacc_id=gnp.pacc_id
    LEFT OUTER JOIN associated_ids AS aid ON aid.gn_id=gnp.gn_id
    LEFT OUTER JOIN prot_desc AS pd ON pd.pacc_id=p.pacc_id
    """
        cursor = self.get_cursor()
        pgdata = {}
        for mid, macc, gacc, aid, cov, pgc, pgnr, desc in cursor.execute(sql):
            pgdata[mid] = {
                    ph.HEADER_PROTEIN: macc,
                    ph.HEADER_GENEID: gacc,
                    ph.HEADER_GENENAME: aid,
                    ph.HEADER_COVERAGE: cov,
                    ph.HEADER_CONTENTPROT: pgc,
                    ph.HEADER_NO_PROTEIN: pgnr,
                    ph.HEADER_DESCRIPTION: desc,
                    }
        return pgdata

    def merge_features(self):
        sql = """
    SELECT bs.set_name, pgm.master_id, COUNT(DISTINCT ppg.psm_id), 
    COUNT (DISTINCT ps.pep_id), COUNT(DISTINCT uni.pep_id), 
    IFNULL(pf.fdr, 'NA'), ppq.quant, fqpsm.amount_psms, GROUP_CONCAT(pqc.channel_name), GROUP_CONCAT(piq.quantvalue),
    GROUP_CONCAT(piq.amount_psms)
        FROM psm_protein_groups AS ppg 
        INNER JOIN psms ON ppg.psm_id=psms.psm_id 
        INNER JOIN mzml ON psms.spectra_id=mzml.spectra_id
        INNER JOIN mzmlfiles ON mzml.mzmlfile_id=mzmlfiles.mzmlfile_id
        INNER JOIN biosets AS bs ON mzmlfiles.set_id=bs.set_id
        INNER JOIN peptide_sequences AS ps ON psms.pep_id=ps.pep_id
        INNER JOIN protein_tables AS pt ON pt.set_id=bs.set_id
        INNER JOIN protein_group_master AS pgm ON pgm.master_id=ppg.master_id
        INNER JOIN proteins AS prots ON prots.pacc_id=pgm.pacc_id
        INNER JOIN protein_fdr AS pf ON pf.prottable_id=pt.prottable_id AND 
            pf.pacc_id=prots.pacc_id
        LEFT OUTER JOIN protein_precur_quanted AS ppq ON ppq.prottable_id=pt.prottable_id AND 
            ppq.pacc_id=prots.pacc_id
        LEFT OUTER JOIN protein_iso_fullpsms AS fqpsm ON fqpsm.prottable_id=pt.prottable_id AND 
            fqpsm.pacc_id=prots.pacc_id
        LEFT OUTER JOIN protquant_channels AS pqc ON pqc.prottable_id=pt.prottable_id
        LEFT OUTER JOIN protein_iso_quanted AS piq ON piq.channel_id=pqc.channel_id AND
            piq.pacc_id=prots.pacc_id
        LEFT OUTER JOIN (
                SELECT ppg.pep_id AS pep_id FROM (
                    SELECT psms.pep_id AS pep_id, COUNT (DISTINCT ppg.master_id) AS nrpep 
                        FROM psm_protein_groups AS ppg INNER JOIN psms USING(psm_id)
                        GROUP BY psms.pep_id
                    ) AS ppg WHERE ppg.nrpep==1
                ) AS uni ON uni.pep_id=ps.pep_id
        GROUP BY pgm.master_id, bs.set_id
        """
        cursor = self.get_cursor()
        cursor.execute(sql)
        return cursor



class GeneTableDB(ProtGeneTableBase):
    datatype = 'gene'
    colmap = {'genes': ['gene_id', 'gene_acc'],
              'gene_precur_quanted': ['gene_id', 'genetable_id', 'quant'],
              'gene_fdr': ['gene_id', 'genetable_id', 'fdr'],
              'gene_iso_fullpsms': ['gene_id', 'genetable_id', 'amount_psms'],
              'genequant_channels': ['channel_id', 'genetable_id',
                                     'channel_name', 'amount_psms_name'],
              'gene_iso_quanted': ['genequant_id', 'gene_id',
                                   'channel_id', 'quantvalue', 'amount_psms'],
              }

    def create_pdata_map(self):
        """This runs only once, returns the data which is not dependent on sets,
        in a dict with accessions as keys"""
        sql = """
SELECT g.gene_acc, GROUP_CONCAT(p.protein_acc, ';'), IFNULL(aid.assoc_id, 'NA'), 
    IFNULL(pd.description, 'NA')
    FROM genes AS g
    LEFT OUTER JOIN ensg_proteins AS egp ON egp.gene_id=g.gene_id
    LEFT OUTER JOIN proteins AS p ON p.pacc_id=egp.pacc_id
    LEFT OUTER JOIN genename_proteins AS gnp ON gnp.pacc_id=egp.pacc_id
    LEFT OUTER JOIN associated_ids AS aid ON aid.gn_id=gnp.gn_id
    LEFT OUTER JOIN prot_desc AS pd ON pd.pacc_id=p.pacc_id
    GROUP BY g.gene_acc
    """
        cursor = self.get_cursor()
        pgdata = {}
        for gacc, pacc, aid, desc in cursor.execute(sql):
            pgdata[gacc] = {
                    ph.HEADER_PROTEINS: pacc,
                    ph.HEADER_GENEID: gacc,
                    ph.HEADER_GENENAME: aid,
                    ph.HEADER_DESCRIPTION: desc,
                    }
        return pgdata

    def merge_features(self):
  ### protein_acc on gene table is indexed??
        sql = """
    SELECT bs.set_name, g.gene_acc, COUNT(DISTINCT ppsm.psm_id), 
    COUNT (DISTINCT ps.pep_id), COUNT(DISTINCT uniq.pep_id), 
    IFNULL(gf.fdr, 'NA'), gpq.quant, fqpsm.amount_psms,
    GROUP_CONCAT(gqc.channel_name), GROUP_CONCAT(giq.quantvalue),
    GROUP_CONCAT(giq.amount_psms)
        FROM protein_psm AS ppsm
        INNER JOIN psms ON ppsm.psm_id=psms.psm_id 
        INNER JOIN mzml ON psms.spectra_id=mzml.spectra_id
        INNER JOIN mzmlfiles ON mzml.mzmlfile_id=mzmlfiles.mzmlfile_id
        INNER JOIN biosets AS bs ON mzmlfiles.set_id=bs.set_id
        INNER JOIN peptide_sequences AS ps ON psms.pep_id=ps.pep_id
        INNER JOIN gene_tables AS gt ON gt.set_id=bs.set_id
        INNER JOIN proteins AS p ON p.protein_acc=ppsm.protein_acc
        INNER JOIN ensg_proteins AS egp ON egp.pacc_id=p.pacc_id
        INNER JOIN genes AS g ON g.gene_id=egp.gene_id
        INNER JOIN gene_fdr AS gf ON gf.genetable_id=gt.genetable_id AND 
            gf.gene_id=g.gene_id
        LEFT OUTER JOIN gene_precur_quanted AS gpq ON gpq.genetable_id=gt.genetable_id AND 
            gpq.gene_id=g.gene_id
        LEFT OUTER JOIN gene_iso_fullpsms AS fqpsm ON fqpsm.genetable_id=gt.genetable_id AND 
            fqpsm.gene_id=g.gene_id
        LEFT OUTER JOIN genequant_channels AS gqc ON gqc.genetable_id=gt.genetable_id
        LEFT OUTER JOIN gene_iso_quanted AS giq ON giq.channel_id=gqc.channel_id AND
            giq.gene_id=g.gene_id

        LEFT OUTER JOIN (
                SELECT psmg.pep_id AS pep_id FROM (
                    SELECT psms.pep_id AS pep_id, COUNT (DISTINCT g.gene_acc) AS nrpep 
                        FROM protein_psm AS ppsm INNER JOIN psms USING(psm_id)
                        INNER JOIN proteins AS p ON p.protein_acc=ppsm.protein_acc
                        INNER JOIN ensg_proteins AS egp ON egp.pacc_id=p.pacc_id
                        INNER JOIN genes AS g ON g.gene_id=egp.gene_id
                        GROUP BY psms.pep_id
                    ) AS psmg WHERE psmg.nrpep==1
                ) AS uniq ON uniq.pep_id=ps.pep_id
        GROUP BY g.gene_id, bs.set_id
        """
        cursor = self.get_cursor()
        cursor.execute(sql)
        return cursor


class GeneTableAssocIDsDB(GeneTableDB):
    datatype = 'assoc'
    colmap = {'associated_ids': ['gn_id', 'assoc_id'],
              'assoc_precur_quanted': ['gn_id', 'genetable_id', 'quant'],
              'assoc_fdr': ['gn_id', 'genetable_id', 'fdr'],
              'assoc_iso_fullpsms': ['gn_id', 'genetable_id', 'amount_psms'],
              'genequant_channels': ['channel_id', 'genetable_id',
                                     'channel_name', 'amount_psms_name'],
              'assoc_iso_quanted': ['genequant_id', 'gn_id',
                                   'channel_id', 'quantvalue', 'amount_psms'],
              }

    def create_pdata_map(self):
        """This runs only once, returns the data which is not dependent on sets,
        in a dict with accessions as keys"""
        sql = """
SELECT gn.assoc_id, GROUP_CONCAT(p.protein_acc, ';'), IFNULL(g.gene_acc, 'NA'), 
    IFNULL(pd.description, 'NA')
    FROM associated_ids AS gn
    LEFT OUTER JOIN genename_proteins AS gnp ON gnp.gn_id=gn.gn_id
    LEFT OUTER JOIN proteins AS p ON p.pacc_id=gnp.pacc_id
    LEFT OUTER JOIN ensg_proteins AS egp ON egp.pacc_id=gnp.pacc_id
    LEFT OUTER JOIN genes AS g ON g.gene_id=egp.gene_id
    LEFT OUTER JOIN prot_desc AS pd ON pd.pacc_id=p.pacc_id
    GROUP BY gn.gn_id
    """
        cursor = self.get_cursor()
        pgdata = {}
        for aid, pacc, gacc, desc in cursor.execute(sql):
            pgdata[aid] = {
                    ph.HEADER_PROTEINS: pacc,
                    ph.HEADER_GENEID: gacc,
                    ph.HEADER_GENENAME: aid,
                    ph.HEADER_DESCRIPTION: desc,
                    }
        return pgdata

    def merge_features(self):
  ### protein_acc on gene table is indexed??
  ### check distinct stuff, pgroups etc need it, matched with content numbers (which cnanot habe it)
        sql = """
    SELECT bs.set_name, gn.assoc_id, COUNT(DISTINCT ppsm.psm_id), 
    COUNT (DISTINCT ps.pep_id), COUNT(DISTINCT uniq.pep_id), 
    IFNULL(gf.fdr, 'NA'), gpq.quant, fqpsm.amount_psms,
    GROUP_CONCAT(gqc.channel_name), GROUP_CONCAT(giq.quantvalue),
    GROUP_CONCAT(giq.amount_psms)
        FROM protein_psm AS ppsm
        INNER JOIN psms ON ppsm.psm_id=psms.psm_id 
        INNER JOIN mzml ON psms.spectra_id=mzml.spectra_id
        INNER JOIN mzmlfiles ON mzml.mzmlfile_id=mzmlfiles.mzmlfile_id
        INNER JOIN biosets AS bs ON mzmlfiles.set_id=bs.set_id
        INNER JOIN peptide_sequences AS ps ON psms.pep_id=ps.pep_id
        INNER JOIN gene_tables AS gt ON gt.set_id=bs.set_id
        INNER JOIN proteins AS p ON p.protein_acc=ppsm.protein_acc
        INNER JOIN genename_proteins AS gnp ON gnp.pacc_id=p.pacc_id
        INNER JOIN associated_ids AS gn ON gn.gn_id=gnp.gn_id
        INNER JOIN assoc_fdr AS gf ON gf.genetable_id=gt.genetable_id AND 
            gf.gn_id=gn.gn_id
        LEFT OUTER JOIN assoc_precur_quanted AS gpq ON gpq.genetable_id=gt.genetable_id AND 
            gpq.gn_id=gn.gn_id
        LEFT OUTER JOIN assoc_iso_fullpsms AS fqpsm ON fqpsm.genetable_id=gt.genetable_id AND 
            fqpsm.gn_id=gn.gn_id
        LEFT OUTER JOIN genequant_channels AS gqc ON gqc.genetable_id=gt.genetable_id
        LEFT OUTER JOIN assoc_iso_quanted AS giq ON giq.channel_id=gqc.channel_id AND
            giq.gn_id=gn.gn_id

        LEFT OUTER JOIN (
                SELECT psmg.pep_id AS pep_id FROM (
                    SELECT psms.pep_id AS pep_id, COUNT (DISTINCT gn.assoc_id) AS nrpep 
                        FROM protein_psm AS ppsm INNER JOIN psms USING(psm_id)
                        INNER JOIN proteins AS p ON p.protein_acc=ppsm.protein_acc
                        INNER JOIN genename_proteins AS gnp ON gnp.pacc_id=p.pacc_id
                        INNER JOIN associated_ids AS gn ON gn.gn_id=gnp.gn_id
                        GROUP BY psms.pep_id
                    ) AS psmg WHERE psmg.nrpep==1
                ) AS uniq ON uniq.pep_id=ps.pep_id
        GROUP BY gn.gn_id, bs.set_id
        """
        cursor = self.get_cursor()
        cursor.execute(sql)
        return cursor
