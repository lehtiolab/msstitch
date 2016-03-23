#!/usr/bin/python3

import sys
from app.drivers import startup
from app.drivers.peptable import (psmtopeptable, merge, model_qvals, isoquant)


def main():
    drivers = [psmtopeptable.MzidTSVPeptableDriver(),
               isoquant.AddIsobaricQuantDriver(),
               model_qvals.ModelQValuesDriver(),
               merge.BuildPeptideTableDriver(),
               ]
    startup.start_msstitch(drivers, sys.argv)
