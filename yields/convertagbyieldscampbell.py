#!/usr/bin/env python
import math
import sys
import numpy as np
from collections import OrderedDict

elements = ['']+[line.split()[1] for line in open('data/atomic_symbols.dat')]

# initial mass fraction

agbYieldFiles = ('m1.25_pmz2e-3.dat','m2.5_pmz2e-3.dat','m3.5_pmz1e-3.dat',
                 'm5.0_nopmz.dat','m6.5_nopmz.dat')
agbModelNames = map(lambda x:x[:-4],agbYieldFiles)

yields = [OrderedDict({}) for k in range(len(agbModelNames))]
ejectaMass = [0 for k in range(len(agbModelNames))]
remnantMass = [0 for k in range(len(agbModelNames))]

for m in range(len(agbYieldFiles)):
    f = open('data/AGByieldsGC/' + agbYieldFiles[m],'rb')
    for line in f.readlines():
        if not line.startswith("#") and not line.startswith("  species"):
            yields[m][line[0:8].strip()] = float(line[31:45]) #absolute yield
            ejectaMass[m] += float(line[31:45])
    remnantMass[m] = float(agbModelNames[m][1:5].rstrip('_')) - ejectaMass[m]

fyields = open('yields.txt', 'w')
fyields.write("#test case: AGB models\n")
fyields.write("#description".ljust(25))
fyields.write("mass".rjust(14))
fyields.write("Z".rjust(14))
fyields.write("remnantmass".rjust(14))
fyields.write("\n[evmodels]\n")
fyields.write(str(len(agbModelNames)) + "\n")

for i in range(len(agbModelNames)):
    fyields.write((agbModelNames[i]).ljust(25))
    fyields.write(("%6.2f" % float(agbModelNames[i][1:4])).rjust(14))   # mass
    fyields.write(("%6.2e" % 0.0016).rjust(14))   # metallicity
    fyields.write(("%6.2f" % remnantMass[i]).rjust(14))   # remnant mass
    fyields.write("\n")

fyields.write("\n" + "#species".ljust(8))
fyields.write("type".rjust(10))
for mn in agbModelNames:
    fyields.write(mn.rjust(14))
fyields.write("\n")

fyields.write("[yields]\n")

for species in yields[0]:
    fyields.write(species.ljust(8))

    elcode = species.rstrip('0123456789').title()
    #if elcode in elements and elements.index(elcode) <= 26:
    fyields.write("absolute".rjust(10))
    #else:
    #    fyields.write("relative".rjust(10))

    for i in range(len(yields)):
        fyields.write(("%14.6e" % yields[i][species]).rjust(14))

    fyields.write("\n")
fyields.close()
