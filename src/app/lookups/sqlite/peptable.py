from app.lookups.sqlite.protpeptable import ProtPepTable
from app.dataformats import peptable as ph
from app.dataformats import prottable as proth


class PepTableProteinCentricDB(ProtPepTable):
    datatype = 'peptide'
    stdheaderfields = [
            ph.HEADER_NO_PSM,
            ph.HEADER_QVAL,
            ph.HEADER_PEP,
            ]
    singlefields = stdheaderfields + [
            ph.HEADER_FALSE_LOC_RATE,
            ph.HEADER_AREA,
            proth.HEADER_NO_FULLQ_PSMS,
            ]
    colmap = {'peptide_sequences': ['pep_id', 'sequence'],
              'peptide_precur_quanted': ['pep_id', 'peptable_id', 'quant'],
              'peptide_fdr': ['pep_id', 'peptable_id', 'fdr'],
              'peptide_pep': ['pep_id', 'peptable_id', 'pep'],
              'ptm_flr': ['pep_id', 'peptable_id', 'flr'],
              'peptide_iso_fullpsms': ['pep_id', 'peptable_id', 'amount_psms'],
              'pepquant_channels': ['channel_id', 'peptable_id',
                                    'channel_name', 'amount_psms_name'],
              'peptide_iso_quanted': ['peptidequant_id', 'pep_id',
                                      'channel_id', 'quantvalue',
                                      'amount_psms'],
              }

    def add_tables(self, tabletypes=[]):
        ttypes = ['fntable', 'isoqtable', 'isochtable', 'prectable',
                'fdrtable', 'peptable', 'fullqpsmtable', 'flrtable']
        self.create_tables([self.table_map[self.datatype][x] for x in ttypes
            if x in self.table_map[self.datatype]])

    def store_pep(self, pep):
        self.store_singlecol('peptable', pep)
        table = self.table_map[self.datatype]['peptable']
        self.singlecols_to_index.append(('pep_acc_ix', table, self.colmap[table][0]))
        self.singlecols_to_index.append(('pep_table_ix', table, self.colmap[table][1]))

    def create_pdata_map(self):
        """This runs only once, returns the data which is not dependent on sets,
        in a dict with accessions as keys"""

        # First protein SQL, some peptides have TOO MANY proteingroup matches 
        # for SQLite to output on a single line, therefore we loop rows 
        # instead of doing GROUP_CONCAT
        protsql = """
        SELECT sub.pep_id, sub.pacc, sub.pgcnr, IFNULL(pd.description, 'NA'), pc.coverage
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
            LEFT OUTER JOIN protein_coverage AS pc ON sub.pacc=pc.protein_acc
            LEFT OUTER JOIN prot_desc AS pd ON sub.pacc_id=pd.pacc_id
            """
        pgdata = {}
        cursor = self.get_cursor()
        for pid, pacc, pgroupnr, desc, cov in cursor.execute(protsql):
            if pid in pgdata:
                pgdata[pid][ph.HEADER_PROTEINS].append(pacc)
                pgdata[pid][ph.HEADER_NO_CONTENTPROTEINS].append(str(pgroupnr))
                pgdata[pid][ph.HEADER_DESCRIPTIONS].append(desc)
                pgdata[pid][ph.HEADER_COVERAGES].append(str(cov))
            else:
                pgdata[pid] = {
                    ph.HEADER_PROTEINS: [pacc],
                    ph.HEADER_NO_CONTENTPROTEINS: [str(pgroupnr)],
                    ph.HEADER_DESCRIPTIONS: [desc],
                    ph.HEADER_COVERAGES: [str(cov)],
                    ph.HEADER_GENES:  [],
                    ph.HEADER_ASSOCIATED: [],
                    }

        # Now get the genes
        ensg_sql = """
            SELECT DISTINCT psms.pep_id, IFNULL(g.gene_acc, 'NA')
            FROM psms
            INNER JOIN psm_protein_groups AS ppg ON psms.psm_id=ppg.psm_id
            INNER JOIN protein_group_master AS pgm ON pgm.master_id=ppg.master_id
            LEFT OUTER JOIN ensg_proteins AS egp ON egp.pacc_id=pgm.pacc_id
            LEFT OUTER JOIN genes AS g ON g.gene_id=egp.gene_id
            """

        gn_sql = """
            SELECT DISTINCT psms.pep_id, IFNULL(aid.assoc_id, 'NA')
            FROM psms
            INNER JOIN psm_protein_groups AS ppg ON psms.psm_id=ppg.psm_id
            INNER JOIN protein_group_master AS pgm ON pgm.master_id=ppg.master_id
            LEFT OUTER JOIN genename_proteins AS gnp ON gnp.pacc_id=pgm.pacc_id
            LEFT OUTER JOIN associated_ids AS aid ON aid.gn_id=gnp.gn_id
            """
        for pid, ensg in cursor.execute(ensg_sql):
            pgdata[pid][ph.HEADER_GENES].append(ensg)
        for pid, gene in cursor.execute(gn_sql):
            pgdata[pid][ph.HEADER_ASSOCIATED].append(gene)

        # Finish with creating strings of the lists
        for pid, headers in pgdata.items():
            for header, vals in headers.items():
                pgdata[pid][header] = ';'.join(vals)

        pepsql = 'SELECT pep_id, sequence FROM peptide_sequences'
        cursor = self.get_cursor()
        for pid, seq in cursor.execute(pepsql):
            pgdata[pid][ph.HEADER_PEPTIDE] = seq
        return pgdata

    def merge_features(self):
        sql = """
    SELECT bs.set_name, ps.pep_id, COUNT(DISTINCT psms.psm_id), pf.fdr, pp.pep,
            flr.flr, ppq.quant, fqpsm.amount_psms, GROUP_CONCAT(pqc.channel_name),
            GROUP_CONCAT(piq.quantvalue), GROUP_CONCAT(piq.amount_psms)
        FROM peptide_sequences AS ps
        INNER JOIN psms ON ps.pep_id=psms.pep_id 
        INNER JOIN mzml ON psms.spectra_id=mzml.spectra_id
        INNER JOIN mzmlfiles ON mzml.mzmlfile_id=mzmlfiles.mzmlfile_id
        INNER JOIN biosets AS bs ON mzmlfiles.set_id=bs.set_id
        INNER JOIN peptide_tables AS pt ON pt.set_id=bs.set_id
        INNER JOIN peptide_fdr AS pf ON pf.peptable_id=pt.peptable_id AND 
            pf.pep_id=ps.pep_id
        LEFT OUTER JOIN peptide_pep AS pp ON pp.peptable_id=pt.peptable_id AND 
            pp.pep_id=ps.pep_id
        LEFT OUTER JOIN ptm_flr AS flr ON flr.peptable_id=pt.peptable_id AND
            flr.pep_id=ps.pep_id
        LEFT OUTER JOIN peptide_precur_quanted AS ppq ON ppq.peptable_id=pt.peptable_id AND
            ppq.pep_id=ps.pep_id
        LEFT OUTER JOIN peptide_iso_fullpsms AS fqpsm ON fqpsm.peptable_id=pt.peptable_id AND
            fqpsm.pep_id=ps.pep_id
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
SELECT ps.pep_id, ps.sequence, GROUP_CONCAT(IFNULL(gsub.ensg, 'NA'), ';'),
    GROUP_CONCAT(IFNULL(gnsub.gn, 'NA'), ';')
    FROM peptide_sequences AS ps
    LEFT OUTER JOIN (
        SELECT DISTINCT psms.pep_id AS pid, g.gene_acc AS ensg
        FROM genes as g
        INNER JOIN ensg_proteins AS egp ON egp.gene_id=g.gene_id
        INNER JOIN proteins AS p ON p.pacc_id=egp.pacc_id
        INNER JOIN protein_psm AS pp ON pp.protein_acc=p.protein_acc
        INNER JOIN psms ON psms.psm_id=pp.psm_id
        ) AS gsub ON gsub.pid=ps.pep_id
    LEFT OUTER JOIN (
        SELECT DISTINCT psms.pep_id AS pid, aid.assoc_id AS gn
        FROM associated_ids AS aid
        INNER JOIN genename_proteins AS gnp ON gnp.gn_id=aid.gn_id
        INNER JOIN proteins AS p ON p.pacc_id=gnp.pacc_id
        INNER JOIN protein_psm AS pp ON pp.protein_acc=p.protein_acc
        INNER JOIN psms ON psms.psm_id=pp.psm_id
        ) AS gnsub ON gnsub.pid=ps.pep_id
    GROUP BY ps.pep_id
    """
        cursor = self.get_cursor()
        pgdata = {}
        for pid, seq, gaccs, aids in cursor.execute(sql):
            pgdata[pid] = {
                    ph.HEADER_PEPTIDE: seq,
                    ph.HEADER_GENES: gaccs,
                    ph.HEADER_ASSOCIATED: aids,
                    }
        return pgdata

    def merge_features(self):
        sql = """
    SELECT bs.set_name, ps.pep_id, COUNT(DISTINCT psms.psm_id), pf.fdr, pp.pep,
            flr.flr, ppq.quant, fqpsm.amount_psms, GROUP_CONCAT(pqc.channel_name),
            GROUP_CONCAT(piq.quantvalue), GROUP_CONCAT(piq.amount_psms)
        FROM peptide_sequences AS ps
        JOIN biosets AS bs 
        INNER JOIN psms ON ps.pep_id=psms.pep_id 
        INNER JOIN peptide_tables AS pt ON pt.set_id=bs.set_id
        INNER JOIN peptide_fdr AS pf ON pf.peptable_id=pt.peptable_id AND 
            pf.pep_id=ps.pep_id
        LEFT OUTER JOIN peptide_pep AS pp ON pp.peptable_id=pt.peptable_id AND 
            pp.pep_id=ps.pep_id
        LEFT OUTER JOIN ptm_flr AS flr ON flr.peptable_id=pt.peptable_id AND
            flr.pep_id=ps.pep_id
        LEFT OUTER JOIN peptide_precur_quanted AS ppq ON ppq.peptable_id=pt.peptable_id AND
            ppq.pep_id=ps.pep_id
        LEFT OUTER JOIN peptide_iso_fullpsms AS fqpsm ON fqpsm.peptable_id=pt.peptable_id AND
            fqpsm.pep_id=ps.pep_id
        LEFT OUTER JOIN pepquant_channels AS pqc ON pqc.peptable_id=pt.peptable_id
        LEFT OUTER JOIN peptide_iso_quanted AS piq ON piq.channel_id=pqc.channel_id AND
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
        sql = """
    SELECT ps.pep_id, ps.sequence, GROUP_CONCAT(ppp.protein_acc, ';')
    FROM peptide_sequences AS ps
    INNER JOIN (
            SELECT DISTINCT psms.pep_id, protein_acc
            FROM psms
            INNER JOIN protein_psm AS pp ON pp.psm_id=psms.psm_id
            ) AS ppp ON ppp.pep_id=ps.pep_id
    GROUP BY ps.pep_id
        """
        cursor = self.get_cursor()
        pgdata = {}
        for pid, seq, prots in cursor.execute(sql):
            pgdata[pid] = {
                    ph.HEADER_PEPTIDE: seq,
                    ph.HEADER_PROTEINS: prots,
                    }
        return pgdata

    def merge_features(self):
        sql = """
    SELECT bs.set_name, ps.pep_id, COUNT(DISTINCT psms.psm_id), pf.fdr, pp.pep,
            flr.flr, ppq.quant, fqpsm.amount_psms, GROUP_CONCAT(pqc.channel_name),
            GROUP_CONCAT(piq.quantvalue), GROUP_CONCAT(piq.amount_psms)
        FROM peptide_sequences AS ps
        JOIN biosets AS bs 
        INNER JOIN psms ON ps.pep_id=psms.pep_id 
        INNER JOIN peptide_tables AS pt ON pt.set_id=bs.set_id
        INNER JOIN peptide_fdr AS pf ON pf.peptable_id=pt.peptable_id AND 
            pf.pep_id=ps.pep_id
        LEFT OUTER JOIN peptide_pep AS pp ON pp.peptable_id=pt.peptable_id AND 
            pp.pep_id=ps.pep_id
        LEFT OUTER JOIN ptm_flr AS flr ON flr.peptable_id=pt.peptable_id AND flr.pep_id=ps.pep_id
        LEFT OUTER JOIN peptide_precur_quanted AS ppq ON ppq.peptable_id=pt.peptable_id AND
            ppq.pep_id=ps.pep_id
        LEFT OUTER JOIN peptide_iso_fullpsms AS fqpsm ON fqpsm.peptable_id=pt.peptable_id AND
            fqpsm.pep_id=ps.pep_id
        LEFT OUTER JOIN pepquant_channels AS pqc ON pqc.peptable_id=pt.peptable_id
        LEFT OUTER JOIN peptide_iso_quanted AS piq ON piq.channel_id=pqc.channel_id AND
            piq.pep_id=ps.pep_id
        GROUP BY ps.pep_id, bs.set_id
        """
        cursor = self.get_cursor()
        cursor.execute(sql)
        return cursor
