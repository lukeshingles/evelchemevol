#!/usr/bin/env python3
from collections import OrderedDict

elements = ('','H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si',
        'P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu',
        'Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc',
        'Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La',
        'Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',
        'Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At',
        'Rn','Fr','Ra','Ac')

if False:
    agbYieldFiles = ('yields-m1.0z001.pmz2e-3.dat',
    'yields-m1.25z001.pmz2e-3.dat',
    'yields-m1.5z001.pmz2e-3.dat',
    'yields-m2.0z001.pmz2e-3.dat',
    'yields-m2.25z001.pmz2e-3.dat',
    'yields-m2.5z001.pmz2e-3.dat',
    'yields-m2.75z001.pmz2e-3.dat',
    'yields-m3.0z001.pmz2e-3.dat',
    'yields-m3.25z001.dat',
    'yields-m3.5z001.dat',
    'yields-m4.0z001.dat',
    'yields-m4.5z001.dat',
    'yields-m5.0z001.dat',
    'yields-m6.0z001.dat',
    'yields-m7.0z001.dat')
else:
    agbYieldFiles = ('yields-m1.0z001.pmz2e-3.dat',
    'yields-m1.25z001.pmz2e-3.dat',
    'yields-m1.5z001.pmz2e-3.dat',
    'yields-m2.0z001.pmz2e-3.dat',
    'yields-m2.25z001.pmz2e-3.dat',
    'yields-m2.5z001.pmz2e-3.dat',
    'yields-m2.75z001.pmz2e-3.dat',
    'yields-m3.0z001.pmz2e-3.dat',
    'yields-m3.25z001.pmz1e-3.dat',
    'yields-m3.5z001.pmz1e-3.dat',
    'yields-m4.0z001.dat',
    'yields-m4.5z001.dat',
    'yields-m5.0z001.dat',
    'yields-m6.0z001.dat',
    'yields-m7.0z001.dat')
agbModelNames = [x[7:-4] for x in agbYieldFiles]

yields = [OrderedDict({}) for k in range(len(agbModelNames))]
ejectaMass = [0.0 for k in range(len(agbModelNames))]
initMass = [float(agbModelNames[k][1:].split('z')[0]) for k in range(len(agbModelNames))]
lifetime = [0.0 for k in range(len(agbModelNames))]

with open('data/fishlock-z001-20140503/lifetimes.dat','r') as flifetimes:
    for line in flifetimes:
        linesplit = line.split()
        for m in range(len(initMass)):
            if "%1.2f" % initMass[m] == linesplit[0]:
                lifetime[m] = float(linesplit[1])

for m in range(len(agbYieldFiles)):
        with open('data/fishlock-z001-20140503/' + agbYieldFiles[m],'r') as fyields:
            print("reading " + agbYieldFiles[m])
            for line in fyields:
                if not line.startswith("#") and not line.startswith("  species"):
                    linesplit = line.split()
                    yields[m][linesplit[0]] = float(linesplit[3]) #absolute yield
                    ejectaMass[m] += float(linesplit[3])

with open('yields.txt', 'w') as fileout:
    print("writing yields.txt")
    fileout.write("#AGB models from Fishlock et al. (2014)\n")
    fileout.write("#modelname".ljust(25))
    fileout.write("mass".rjust(14))
    fileout.write("Z".rjust(14))
    fileout.write("remnantmass".rjust(14))
    fileout.write("lifetime".rjust(14))
    fileout.write("\n[stellarmodels]\n")
    fileout.write(str(len(agbModelNames)) + "\n")

    for i in range(len(agbModelNames)):
        fileout.write((agbModelNames[i]).ljust(25))
        fileout.write(("%6.2f" % float(agbModelNames[i][1:5].rstrip('z'))).rjust(14)) # mass
        fileout.write(("%6.2e" % 0.001).rjust(14)) # metallicity
        fileout.write(("%6.3f" % (initMass[i]  - ejectaMass[i])).rjust(14))   # remnant mass
        fileout.write(("%.6e" % lifetime[i]).rjust(14))
        fileout.write("\n")

    fileout.write("\n" + "#species".ljust(8))
    fileout.write("type".rjust(10))
    for mn in agbModelNames:
        fileout.write(mn[:13].rjust(14))
    fileout.write("\n")

    fileout.write("[yields]\n")

    for species in yields[0]:
        fileout.write(species.ljust(8))

        elcode = species.rstrip('0123456789').title()
        if elcode in elements and elements.index(elcode) <= 26:
            fileout.write("absolute".rjust(10))
        else:
            fileout.write("relative".rjust(10))

        for i in range(len(yields)):
            fileout.write(("%14.6e" % yields[i][species]).rjust(14))

        fileout.write("\n")
