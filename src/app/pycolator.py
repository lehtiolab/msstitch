#!/usr/bin/python3

import sys

from app.drivers.pycolator import (splitmerge, filters, stats, qvality)
from app.drivers import startup


def main():
    drivers = [splitmerge.SplitTDDriver(),
               splitmerge.SplitProteinDriver(),
               splitmerge.MergeDriver(),
               filters.FilterUniquePeptides(),
               filters.FilterPeptideLength(),
               filters.FilterWholeProteinSequence(),
               filters.FilterPeptideSequence(),
               qvality.QvalityDriver(),
               stats.ReassignmentDriver(),
               ]
    startup.start_msstitch(drivers, sys.argv)
