#!/usr/bin/env python3
import os

scriptdir = os.path.dirname(os.path.abspath(__file__))

# read the solar elemental abundances
solarmetallicity = 0.015
elsolarlogepsilon = {}
with open(os.path.join(scriptdir, 'solar_abund_ARAAedit.dat'), 'r') as solarinfile:
    for line in solarinfile:
        if len(line.split(' ')) >= 4:
            row = [line[0:4], line[4:8], line[8:16], line[16:24]]
            elcode = row[1].strip(' ')
            if elcode == 'p' and int(row[0]) == 1:
                elcode = 'h'
            elsolarlogepsilon[elcode] = float(row[3])
