#!/usr/bin/python3

import sys

from app.drivers.pycolator import (splitmerge, filters)
from app.drivers import startup


def main():
    drivers = [splitmerge.SplitProteinDriver(),
               filters.FilterWholeProteinSequence(),
               filters.FilterPeptideSequence(),
               ]
    startup.start_msstitch(drivers, sys.argv)
