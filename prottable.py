#!/usr/bin/python3

import sys
from app.drivers import startup
from app.drivers.prottable import (probability, info, merge, isoquant,
                                   precursorarea, create_empty, bestpeptide,
                                   qvality, fdr)


drivers = [info.AddProteinInfoDriver(),
           merge.BuildProteinTableDriver(),
           isoquant.AddIsobaricQuantDriver(),
           precursorarea.AddPrecursorAreaDriver(),
           probability.AddProteinProbability(),
           create_empty.CreateEmptyDriver(),
           bestpeptide.BestPeptidePerProtein(),
           qvality.ProttableQvalityDriver(),
           qvality.PickedQvalityDriver(),
           fdr.ProttableFDRDriver(),
           ]
startup.start_msstitch(drivers, sys.argv)
