from app.drivers import base
from app.readers import tsv as tsvreader
from app.actions import merge

from app.drivers.options import prottable_options, peptable_options, lookup_options
from app.dataformats import prottable as ph 
from app.dataformats import peptable as peph

class MergeDriver(base.PepProttableDriver):
    outsuffix = ''
    command = 'merge'
    commandhelp = ('Create merged gene/protein/peptide tables from multiple '
            'sample set tables')

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['setnames',
            'quantcolpattern', 'precursorquantcolpattern', 'fdrcolpattern',
            'flrcolpattern', 'multifiles', 'nogroup', 
            'featcol'], lookup_options))
        self.options.update(self.define_options(['lookupfn', 'mergecutoff'],
            prottable_options))
        self.options.update(self.define_options(['pepcolpattern',
            'genecentric'], peptable_options))

    def parse_input(self, **kwargs):
        super().parse_input(**kwargs)
        header = tsvreader.get_tsv_header(self.fn[0])
        self.header = [header[0]]
        if header[0] != peph.HEADER_PEPTIDE:
            self.is_peptidetable = False
            self.pepcolpattern = None # override input if any
        #else:
        if header[0] == peph.HEADER_PEPTIDE:
            self.header.append(peph.HEADER_BAREPEP)
            self.is_peptidetable = True
            if self.genecentric:
                self.lookuptype = 'peptidegenecentrictable'
                self.header.extend([peph.HEADER_GENES, peph.HEADER_ASSOCIATED])
            elif self.nogroup:
                self.lookuptype = 'peptidetableplain'
                self.header.extend([peph.HEADER_PROTEINS, peph.HEADER_STARTSTOP])
            else:
                self.header.extend([
                    peph.HEADER_PROTEINS, peph.HEADER_STARTSTOP, peph.HEADER_NO_CONTENTPROTEINS, 
                    peph.HEADER_DESCRIPTIONS, peph.HEADER_COVERAGES, peph.HEADER_GENES, 
                    peph.HEADER_ASSOCIATED])
                self.lookuptype = 'peptidetable'
        elif header[0] == ph.HEADER_PROTEIN and self.nogroup:
            self.lookuptype = 'plainprottable'
        elif header[0] == ph.HEADER_PROTEIN:
            self.lookuptype = 'prottable'
            self.header.extend([
                ph.HEADER_GENEID, ph.HEADER_GENENAME, ph.HEADER_DESCRIPTION,
                ph.HEADER_COVERAGE, ph.HEADER_CONTENTPROT, ph.HEADER_NO_PROTEIN])
        elif header[0] == ph.HEADER_GENENAME:
            self.lookuptype = 'associdtable'
            self.header.extend([
                ph.HEADER_GENEID, ph.HEADER_PROTEINS, ph.HEADER_DESCRIPTION])
        elif header[0] == ph.HEADER_GENEID:
            self.lookuptype = 'genetable'
            self.header.extend([
                ph.HEADER_GENENAME, ph.HEADER_PROTEINS, ph.HEADER_DESCRIPTION])

    def set_features(self):
        # FIXME need check to see same poolnames correlate with self.fn len
        self.lookup.add_tables(self.tabletypes)
        self.featcol = self.featcol - 1
        self.setnames = [x.replace('"', '') for x in self.setnames]
        merge.create_lookup(self.fn, self.lookup, self.setnames, self.featcol,
                self.precursorquantcolpattern, self.quantcolpattern,
                self.fdrcolpattern, self.pepcolpattern, self.flrcolpattern)
        for field in self.lookup.stdheaderfields:
            self.header.extend(['{}_{}'.format(x, field) for x in self.setnames])
        if self.fdrcolpattern and not self.is_peptidetable:
            self.header.extend(['{}_{}'.format(x, ph.HEADER_QVAL) for x in self.setnames])
        if self.is_peptidetable:
            if self.flrcolpattern:
                self.header.extend(['{}_{}'.format(x, peph.HEADER_FALSE_LOC_RATE) for x in self.setnames])
            if self.precursorquantcolpattern:
                self.header.extend(['{}_{}'.format(x, peph.HEADER_AREA) for x in self.setnames])
        elif self.precursorquantcolpattern:
            self.header.extend(['{}_{}'.format(x, ph.HEADER_AREA) for x in self.setnames])
        if self.quantcolpattern:
            stored_psmnrs = self.lookup.check_isoquant_psmnrs()
            channels = [x for x in self.lookup.get_isoquant_headernames(stored_psmnrs)]
            self.header.extend([x[0] for x in channels])
            if stored_psmnrs:
                self.header.extend([f'{setn}_{ph.HEADER_NO_FULLQ_PSMS}' for setn in self.setnames])
                self.header.extend([x[1] for x in channels])
        self.features = merge.build_proteintable(self.lookup, self.mergecutoff, self.quantcolpattern)
