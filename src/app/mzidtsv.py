#!/usr/bin/python3

import sys

from app.drivers.mzidtsv import (spectra, percolator, proteingrouping,
                                 prot2gene, quant, splitmerge,
                                 filter_confidence, isonormalize)
from app.drivers import startup


def main():
    drivers = [spectra.TSVSpectraDriver(),
               percolator.MzidPercoTSVDriver(),
               splitmerge.MzidTSVConcatenateDriver(),
               splitmerge.MzidTSVSplitDriver(),
               filter_confidence.ConfidenceFilterDriver(),
               quant.TSVQuantDriver(),
               proteingrouping.ProteinGroupDriver(),
               prot2gene.TSVGeneFromProteinDriver(),
               isonormalize.MzidTSVIsoquantNormalizeDriver(),
               ]
    startup.start_msstitch(drivers, sys.argv)
