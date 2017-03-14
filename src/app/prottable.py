#!/usr/bin/python3

import sys
from app.drivers import startup
from app.drivers.prottable import (probability, info, merge, isoquant, fdr,
                                   precursorarea, create_empty, bestpeptide)


def main():
    drivers = [info.AddProteinInfoDriver(),
               merge.BuildProteinTableDriver(),
               isoquant.AddIsobaricQuantDriver(),
               precursorarea.AddPrecursorAreaDriver(),
               probability.AddProteinProbability(),
               create_empty.CreateEmptyDriver(),
               bestpeptide.BestPeptidePerProtein(),
               fdr.ProttableFDRDriver(),
               fdr.PickedFDRDriver(),
               fdr.ProttableAddFDRDriver(),
               ]
    startup.start_msstitch(drivers, sys.argv)
