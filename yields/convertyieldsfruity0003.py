#!/usr/bin/env python3
from collections import OrderedDict
#import glob
import os

elements = ['']+[line.split()[1] for line in open('data/atomic_symbols.dat')]

#(yval, k14alpha, s15alpha) = ('24', '4', '0')

yieldfilenames = ['yields_tot_m1p5z3m4_000_20140820_73609.txt',
                  'yields_tot_m2p0z3m4_000_20140820_73609.txt',
                  'yields_tot_m2p5z3m4_000_20140820_73609.txt',
                  'yields_tot_m3p0z3m4_000_20140820_73609.txt',
                  'yields_tot_m4p0z3m4_000_20140820_73609.txt',
                  'yields_tot_m5p0z3m4_000_20140820_73609.txt',
                  'yields_tot_m6p0z3m4_000_20140820_73609.txt']
modelNames = list(map(lambda x: x.split('_')[2],yieldfilenames))
initMass = [1.5, 2.0, 2.5, 3, 4, 5, 6]

yields = [OrderedDict({}) for k in range(len(modelNames))]
ejectaMass = [0.0 for k in range(len(modelNames))]
metallicity = [3e-4 for k in range(len(modelNames))]
lifetime = [
    1.85e9, #1.5 Msun
    8.46e8, #2.0 Msun
    4.92e8, #2.5 Msun
    3.02e8, #3.0 Msun
    1.48e8, #4.0 Msun
    8.93e7, #5.0 Msun
    6.06e7, #6.0 Msun
]


for m, modelname in enumerate(modelNames):
    evmodelname = modelNames[m].split('a')[0]

    with open(os.path.join('data/fruity-z0003',yieldfilenames[m]),'r') as fyields:
        print("reading " + yieldfilenames[m])
        fyields.readline()
        for line in fyields:
            row = line.split()
            elmsymbol = row[0].rstrip('0123456789')

            # decays of unstable nuclides
            if row[0] == 'Zr93': elmsymbol = 'Nb'
            if row[0] == 'Tc99': elmsymbol = 'Ru'
            if row[0] == 'Pd07': elmsymbol = 'Ag'
            if row[0] == 'Cs35': elmsymbol = 'Ba'
            if row[0] == 'Cs37': elmsymbol = 'Ba'
            if row[0] == 'Pb05': elmsymbol = 'Tl'
            if row[0] == 'Bi10': elmsymbol = 'Pb'
            if row[0] == 'Po10': elmsymbol = 'Pb'

            massnumber = int(row[1])    #(A) the number of nucleons
            massyield = float(row[3])   #absolute stellar yield in solar masses
            ejectaMass[m] += massyield

            speciesname = '{:}{:d}'.format(elmsymbol.lower(),massnumber)

            if elmsymbol == 'H':
                elmsymbol = 'p'
                speciesname = 'p'

            if speciesname in yields[m]:
                yields[m][speciesname] = yields[m][speciesname] + massyield
            else:
                yields[m][speciesname] = massyield

with open('yields.txt', 'w') as fileout:
    print("writing yields.txt")
    fileout.write("#test case: FRUITY database models from Straniero et al. 2014 \n")
    fileout.write("#modelname".ljust(25))
    fileout.write("mass".rjust(14))
    fileout.write("Z".rjust(14))
    fileout.write("remnantmass".rjust(14))
    fileout.write("lifetime".rjust(14))
    fileout.write("\n[stellarmodels]\n")
    fileout.write(str(len(modelNames)) + "\n")

    for i in range(len(modelNames)):
        fileout.write((chr(ord('A') + i) + ':' + modelNames[i]).ljust(25)[:25])
        fileout.write(("{0:6.2f}".format(initMass[i])).rjust(14))   # mass
        fileout.write(("{0:6.2e}".format(metallicity[i])).rjust(14))   # metallicity
        fileout.write(("{0:6.3f}".format((initMass[i]  - ejectaMass[i]))).rjust(14))   # remnant mass
        fileout.write(("{0:.6e}".format(lifetime[i])).rjust(14))
        fileout.write("\n")

    fileout.write("\n" + "#species".ljust(8))
    fileout.write("type".rjust(10))
    for m in range(len(modelNames)):
#        shortname = (modelNames[m].split('z')[0] + 'y' + modelNames[m].split('y')[1])
#        shortname = shortname.split('s')[0]
#        fileout.write(shortname.rjust(14)[:14])
        fileout.write(chr(ord('A') + m).rjust(14))
    fileout.write("\n")

    fileout.write("[yields]\n")

    for species in yields[0]:
        fileout.write(species.ljust(8))

        elcode = species.rstrip('0123456789').title()
        fileout.write("absolute".rjust(10))

        for i in range(len(yields)):
            if species in yields[i]:
                yieldout = yields[i][species]
            else:
                yieldout = 0.0
            fileout.write(("{0:14.6e}".format(yieldout)).rjust(14))

        fileout.write("\n")
