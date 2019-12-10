from app.actions.mzidtsv import percolator as prep
from app.drivers.mzidtsv import MzidTSVDriver
from app.drivers.options import mzidtsv_options
from app.readers import tsv as tsvreader
from app.readers import mzidplus as mzidreader
from app.writers import tsv as writer


class MzidPercoTSVDriver(MzidTSVDriver):
    """
    Adds percolator data from mzid file to table.
    """
    outsuffix = '_fdr.tsv'
    command = 'percolator'
    commandhelp = (
            'Calculates FDR from percolator output and adds FDR and percolator '
            ' to the corresponding TSV PSM tables. FDR calculation method is TD-'
            'competition.'
            )

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['multifiles', 'mzidfns', 'percofn'], mzidtsv_options))

    def get_psms(self):
        self.percopsms = prep.calculate_target_decoy_competition(self.percofn)
                

    def write(self):
        # FIXME 'Strip', 'Fraction', 'missed_cleavage'
        for psmfn, mzidfn in zip(self.fn, self.mzidfns):
            oldheader = tsvreader.get_tsv_header(psmfn)
            header = prep.get_header_with_percolator(oldheader)
            outfn = self.create_outfilepath(psmfn, self.outsuffix)
            psms = tsvreader.generate_tsv_psms(psmfn, oldheader)
            mzns = mzidreader.get_mzid_namespace(mzidfn)
            mzidsr = mzidreader.mzid_spec_result_generator(mzidfn, mzns)
            psms_out = prep.add_fdr_to_mzidtsv(psms, mzidsr, mzns, self.percopsms)
            writer.write_tsv(header, psms_out, outfn)

