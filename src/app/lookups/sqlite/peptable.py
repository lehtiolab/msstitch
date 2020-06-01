from app.lookups.sqlite.protpeptable import ProtPepTable
from app.dataformats import peptable as ph


class PepTableProteinCentricDB(ProtPepTable):
    datatype = 'peptide'
    stdheaderfields = [
            ph.HEADER_NO_PSM,
            ph.HEADER_QVAL,
            ]
    singlefields = stdheaderfields + [
            ph.HEADER_AREA,
            ]
    colmap = {'peptide_sequences': ['pep_id', 'sequence'],
              'peptide_precur_quanted': ['pep_id', 'peptable_id', 'quant'],
              'peptide_fdr': ['pep_id', 'peptable_id', 'fdr'],
              'pepquant_channels': ['channel_id', 'peptable_id',
                                    'channel_name', 'amount_psms_name'],
              'peptide_iso_quanted': ['peptidequant_id', 'pep_id',
                                      'channel_id', 'quantvalue',
                                      'amount_psms'],
              }

    def add_tables(self, tabletypes=[]):
        self.create_tables(['peptide_tables', 'pepquant_channels',
                            'peptide_iso_quanted', 'peptide_precur_quanted',
                            'peptide_fdr'])

    def create_pdata_map(self):
        """This runs only once, returns the data which is not dependent on sets,
        in a dict with accessions as keys"""
        sql = """
    SELECT sub.pep_id, sub.seq, GROUP_CONCAT(sub.pacc, ';'),
    GROUP_CONCAT(sub.pgcnr, ';'),
    GROUP_CONCAT(IFNULL(pd.description, 'NA'), ';'),
    GROUP_CONCAT(pc.coverage, ';'), 
    GROUP_CONCAT(IFNULL(gsub.ensg, 'NA'), ';'), GROUP_CONCAT(IFNULL(gnsub.gn, 'NA'), ';')
    FROM (
        SELECT psms.pep_id AS pep_id, ps.sequence AS seq, pgm.master_id,
            pgm.pacc_id AS pacc_id, p.protein_acc AS pacc, COUNT(pgc.protein_acc) AS pgcnr
        FROM psm_protein_groups AS pp
        INNER JOIN psms ON psms.psm_id=pp.psm_id
        INNER JOIN protein_group_master AS pgm ON pgm.master_id=pp.master_id
        INNER JOIN proteins AS p ON p.pacc_id=pgm.pacc_id
        INNER JOIN protein_group_content AS pgc ON pgm.master_id=pgc.master_id
        INNER JOIN peptide_sequences AS ps ON psms.pep_id=ps.pep_id
        GROUP BY psms.pep_id, pgm.master_id
        ) AS sub
    LEFT OUTER JOIN (
        SELECT DISTINCT psms.pep_id AS pid, g.gene_acc AS ensg
        FROM genes as g
        INNER JOIN ensg_proteins AS egp ON egp.gene_id=g.gene_id
        INNER JOIN proteins AS p ON p.pacc_id=egp.pacc_id
        INNER JOIN protein_group_master AS pgm ON pgm.pacc_id=egp.pacc_id
        INNER JOIN psm_protein_groups AS ppg ON pgm.master_id=ppg.master_id
        INNER JOIN psms ON psms.psm_id=ppg.psm_id
        ) AS gsub ON gsub.pid=sub.pep_id
    LEFT OUTER JOIN (
        SELECT DISTINCT psms.pep_id AS pid, aid.assoc_id AS gn
        FROM associated_ids AS aid
        INNER JOIN genename_proteins AS gnp ON gnp.gn_id=aid.gn_id
        INNER JOIN proteins AS p ON p.pacc_id=gnp.pacc_id
        INNER JOIN protein_group_master AS pgm ON pgm.pacc_id=p.pacc_id
        INNER JOIN psm_protein_groups AS ppg ON pgm.master_id=ppg.master_id
        INNER JOIN psms ON psms.psm_id=ppg.psm_id
        ) AS gnsub ON gnsub.pid=sub.pep_id
    LEFT OUTER JOIN protein_coverage AS pc ON sub.pacc=pc.protein_acc
    LEFT OUTER JOIN prot_desc AS pd ON sub.pacc_id=pd.pacc_id
    GROUP BY sub.pep_id
    """
        cursor = self.get_cursor()
        pgdata = {}
        for pid, seq, paccs, pgroupnrs, descs, covs, gaccs, aids in cursor.execute(sql):
            pgdata[pid] = {
                    ph.HEADER_PEPTIDE: seq,
                    ph.HEADER_PROTEINS: paccs,
                    ph.HEADER_NO_CONTENTPROTEINS: pgroupnrs,
                    ph.HEADER_DESCRIPTIONS: descs,
                    ph.HEADER_COVERAGES: covs,
                    ph.HEADER_GENES: gaccs,
                    ph.HEADER_ASSOCIATED: aids,
                    }
        return pgdata

    def merge_features(self):
        sql = """
    SELECT bs.set_name, ps.pep_id, COUNT(DISTINCT psms.psm_id), pf.fdr, ppq.quant,
            GROUP_CONCAT(pqc.channel_name), GROUP_CONCAT(piq.quantvalue),
            GROUP_CONCAT(piq.amount_psms)
        FROM peptide_sequences AS ps
        JOIN biosets AS bs 
        INNER JOIN psms ON ps.pep_id=psms.pep_id 
        INNER JOIN peptide_tables AS pt ON pt.set_id=bs.set_id
        INNER JOIN peptide_fdr AS pf ON pf.peptable_id=pt.peptable_id AND 
            pf.pep_id=ps.pep_id
        LEFT OUTER JOIN peptide_precur_quanted AS ppq ON ppq.peptable_id=pt.peptable_id AND 
            ppq.pep_id=ps.pep_id
        LEFT OUTER JOIN pepquant_channels AS pqc ON pqc.peptable_id=pt.peptable_id
        LEFT OUTER JOIN peptide_iso_quanted AS piq ON piq.channel_id=pqc.channel_id AND
            piq.pep_id=ps.pep_id
        GROUP BY ps.pep_id, bs.set_id
        """
        cursor = self.get_cursor()
        cursor.execute(sql)
        return cursor



class PepTableGeneCentricDB(PepTableProteinCentricDB):
    datatype = 'peptide'

    def create_pdata_map(self):
        """This runs only once, returns the data which is not dependent on sets,
        in a dict with accessions as keys. Different from protein grouped
        data in that we dont pull proteins from the protein group master table.
        """

        sql = """
SELECT ps.pep_id, ps.sequence, GROUP_CONCAT(p.protein_acc, ';'),
        GROUP_CONCAT(IFNULL(pd.description, 'NA'), ';'), 
        IFNULL(gsub.ensg, 'NA'), IFNULL(gnsub.gn, 'NA')
    FROM protein_psm AS pp
    INNER JOIN psms ON psms.psm_id=pp.psm_id
    INNER JOIN peptide_sequences AS ps ON psms.pep_id=ps.pep_id
    INNER JOIN proteins AS p ON p.protein_acc=pp.protein_acc
    LEFT OUTER JOIN (
        SELECT gss.pid AS pid, gss.ensg AS ensg FROM (
            SELECT DISTINCT psms.pep_id AS pid, g.gene_acc AS ensg
            FROM genes as g
            INNER JOIN ensg_proteins AS egp ON egp.gene_id=g.gene_id
            INNER JOIN proteins AS p ON p.pacc_id=egp.pacc_id
            INNER JOIN protein_psm AS pp ON pp.protein_acc=p.protein_acc
            INNER JOIN psms ON psms.psm_id=pp.psm_id
            ) AS gss GROUP BY gss.pid
        ) AS gsub ON gsub.pid=ps.pep_id
    LEFT OUTER JOIN (
        SELECT gnss.pid AS pid, GROUP_CONCAT(gnss.gn) AS gn FROM (
            SELECT DISTINCT psms.pep_id AS pid, aid.assoc_id AS gn
            FROM associated_ids AS aid
            INNER JOIN genename_proteins AS gnp ON gnp.gn_id=aid.gn_id
            INNER JOIN proteins AS p ON p.pacc_id=gnp.pacc_id
            INNER JOIN protein_psm AS pp ON pp.protein_acc=p.protein_acc
            INNER JOIN psms ON psms.psm_id=pp.psm_id
            ) AS gnss
            GROUP BY gnss.pid
        ) AS gnsub ON gnsub.pid=ps.pep_id
    LEFT OUTER JOIN prot_desc AS pd ON p.pacc_id=pd.pacc_id
    LEFT OUTER JOIN protein_coverage AS pc ON p.protein_acc=pc.protein_acc
    GROUP BY ps.pep_id
    """
        cursor = self.get_cursor()
        pgdata = {}
        for pid, seq, paccs, descs, gaccs, aids in cursor.execute(sql):
            pgdata[pid] = {
                    ph.HEADER_PEPTIDE: seq,
                    ph.HEADER_GENES: gaccs,
                    ph.HEADER_ASSOCIATED: aids,
                    }
        return pgdata

    def merge_features(self):
        sql = """
    SELECT bs.set_name, ps.pep_id, pf.fdr, ppq.quant, COUNT(DISTINCT psms.psm_id),
            GROUP_CONCAT(pqc.channel_name), GROUP_CONCAT(piq.quantvalue),
            GROUP_CONCAT(piq.amount_psms)
        FROM peptide_sequences AS ps
        JOIN biosets AS bs 
        INNER JOIN psms ON ps.pep_id=psms.pep_id 
        INNER JOIN peptide_tables AS pt ON pt.set_id=bs.set_id
        INNER JOIN peptide_fdr AS pf ON pf.peptable_id=pt.peptable_id AND 
            pf.pep_id=ps.pep_id
        INNER JOIN peptide_precur_quanted AS ppq ON ppq.peptable_id=pt.peptable_id AND 
            ppq.pep_id=ps.pep_id
        INNER JOIN pepquant_channels AS pqc ON pqc.peptable_id=pt.peptable_id
        INNER JOIN peptide_iso_quanted AS piq ON piq.channel_id=pqc.channel_id AND
            piq.pep_id=ps.pep_id
        GROUP BY ps.pep_id, bs.set_id
        """
        cursor = self.get_cursor()
        cursor.execute(sql)
        return cursor



class PepTablePlainDB(PepTableProteinCentricDB):
    datatype = 'peptide'

    def create_pdata_map(self):
        """This runs only once, returns the data which is not dependent on sets,
        in a dict with accessions as keys. 
        In plain DB we only output peptides, not proteins etc
        """
        sql = 'SELECT ps.pep_id, ps.sequence FROM peptide_sequences AS ps'
        cursor = self.get_cursor()
        pgdata = {}
        for pid, seq in cursor.execute(sql):
            pgdata[pid] = {
                    ph.HEADER_PEPTIDE: seq,
                    }
        return pgdata

    def merge_features(self):
        sql = """
    SELECT bs.set_name, ps.pep_id, pf.fdr, ppq.quant, COUNT(DISTINCT psms.psm_id),
            GROUP_CONCAT(pqc.channel_name), GROUP_CONCAT(piq.quantvalue),
            GROUP_CONCAT(piq.amount_psms)
        FROM peptide_sequences AS ps
        JOIN biosets AS bs 
        INNER JOIN psms ON ps.pep_id=psms.pep_id 
        INNER JOIN peptide_tables AS pt ON pt.set_id=bs.set_id
        INNER JOIN peptide_fdr AS pf ON pf.peptable_id=pt.peptable_id AND 
            pf.pep_id=ps.pep_id
        INNER JOIN peptide_precur_quanted AS ppq ON ppq.peptable_id=pt.peptable_id AND 
            ppq.pep_id=ps.pep_id
        INNER JOIN pepquant_channels AS pqc ON pqc.peptable_id=pt.peptable_id
        INNER JOIN peptide_iso_quanted AS piq ON piq.channel_id=pqc.channel_id AND
            piq.pep_id=ps.pep_id
        GROUP BY ps.pep_id, bs.set_id
        """
        cursor = self.get_cursor()
        cursor.execute(sql)
        return cursor
