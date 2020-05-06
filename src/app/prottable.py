#!/usr/bin/python3

import sys
from app.drivers import startup
from app.drivers.prottable import (merge, fdr, precursorarea, 
        create_empty, bestpeptide)


def main():
    drivers = [merge.BuildProteinTableDriver(),
               precursorarea.AddPrecursorAreaDriver(),
               create_empty.CreateEmptyDriver(),
               bestpeptide.BestPeptidePerProtein(),
               fdr.ProttableFDRDriver(),
               fdr.PickedFDRDriver(),
               ]
    startup.start_msstitch(drivers, sys.argv)
