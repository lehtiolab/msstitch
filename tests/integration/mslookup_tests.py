from tests.integration import basetests

import shutil
import os
import sqlite3
import yaml


class TestTrypticLookup(basetests.BaseTest):
    command = 'seqspace'
    suffix = ''
    executable = 'mslookup.py'
    infilename = 'proteins.fasta'

    def all_seqs_in_db(self, dbfn, sequences, seqtype):
        db = sqlite3.connect(dbfn)
        print(dbfn)
        print(db)

        seqs_in_db = set()
        for seq in sequences:
            seqs_in_db.add(self.seq_in_db(db, seq, seqtype))
        db.close()
        return seqs_in_db == set([True])

    def query_db_assert(self, options=None, seqtype=None):
        if options is None:
            options = []
        with open(os.path.join(self.fixdir, 'peptides_trypsinized.yml')) as fp:
            tryp_sequences = yaml.load(fp)
        sequences = tryp_sequences['fully_tryptic']
        if seqtype is not None:
            sequences.extend(tryp_sequences[seqtype])
        self.run_command(options)
        print(os.listdir(self.workdir))
        self.assertTrue(self.all_seqs_in_db(self.resultfn,
                                            sequences, seqtype))

    def run_without_db(self, options=None, seqtype=None):
        self.resultfn = os.path.join(self.workdir,
                                     'msstitch_searchspace_db.sqlite')
        self.query_db_assert(options, seqtype)

    def run_with_existing_db(self, options=None, seqtype=None):
        if options is None:
            options = []
        self.resultfn = os.path.join(self.workdir, 'seqspace.db')
        options.extend(['--dbfile', self.resultfn])
        shutil.copy(os.path.join(self.fixdir, 'mzidtsv_db.sqlite'),
                    self.resultfn)
        self.query_db_assert(options, seqtype)

    def test_cutproline_nodb(self):
        self.run_without_db(['--cutproline'], 'proline_cuts')

    def test_cutproline_yesdb(self):
        self.run_with_existing_db(['--cutproline'], 'proline_cuts')

    def test_ntermwildcards_nodb(self):
        self.run_without_db(['--ntermwildcards'], 'ntermfalloff')

    def test_ntermwildcards_yes_db(self):
        self.run_with_existing_db(['--ntermwildcards'], 'ntermfalloff')

    def test_noflags_no_db(self):
        self.run_without_db()

    def test_noflags_yes_db(self):
        self.run_with_existing_db()
