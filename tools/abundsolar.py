#!/usr/bin/env python3

# read the solar elemental abundances
solarmetallicity = 0.015
elsolarlogepsilon = {}
with open('solar_abund_Asplund09.dat','r') as solarinfile:
    for line in solarinfile:
        if len(line.split(' ')) >= 4:
            row = [line[0:4], line[4:8], line[8:16], line[16:24]]
            elcode = row[1].strip(' ')
            if elcode == 'p' and int(row[0]) == 1:
                elcode = 'h'
            elsolarlogepsilon[elcode] = float(row[3])
