#!/usr/bin/python3

import sys

from app.drivers.pycolator import (splitmerge, filters, stats, qvality)
from app.drivers import startup


def main():
    drivers = [splitmerge.SplitDriver(),
               splitmerge.MergeDriver(),
               filters.FilterUniquePeptides(),
               filters.FilterPeptideLength(),
               filters.FilterKnownPeptides(),
               qvality.QvalityDriver(),
               stats.ReassignmentDriver(),
               ]
    startup.start_msstitch(drivers, sys.argv)
