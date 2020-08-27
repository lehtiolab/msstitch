#!/usr/bin/python3

import sys
from app.drivers import (startup, lookups, sequence, psmtable, percolator,
        peptable, prottable, merge)


def main():
    drivers = [lookups.SpectraLookupDriver(),
               lookups.QuantLookupDriver(),
               lookups.SequenceLookupDriver(),

               sequence.DecoySeqDriver(),
               sequence.TrypsinizeDriver(),

               psmtable.TSVConcatenateDriver(),
               psmtable.TSVSplitDriver(),
               psmtable.Perco2PSMDriver(),
               psmtable.ConfidenceFilterDriver(),
               psmtable.PSMTableRefineDriver(),
               psmtable.IsoSummarizeDriver(),
               psmtable.DeleteSetDriver(),

               percolator.FilterSequences(),
               percolator.SplitProteinDriver(),

               peptable.CreatePeptableDriver(),
               prottable.ProteinsDriver(),
               prottable.GenesDriver(),
               prottable.ENSGDriver(),
               
               merge.MergeDriver(),
               ]
    startup.start_msstitch(drivers, sys.argv)
