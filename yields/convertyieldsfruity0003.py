#!/usr/bin/env python3
from collections import OrderedDict
# import glob
import os

from enum import Enum
yield_types = Enum('yield_types', 'relative absolute')

yield_type = yield_types.relative

elements = ['']+[line.split()[1] for line in open('data/atomic_symbols.dat')]

# (yval, k14alpha, s15alpha) = ('24', '4', '0')

total_yield_filenames = [
    'yields_tot_m1p5z3m4_000_20140820_73609.txt',
    'yields_tot_m2p0z3m4_000_20140820_73609.txt',
    'yields_tot_m2p5z3m4_000_20140820_73609.txt',
    'yields_tot_m3p0z3m4_000_20140820_73609.txt',
    'yields_tot_m4p0z3m4_000_20140820_73609.txt',
    'yields_tot_m5p0z3m4_000_20140820_73609.txt',
    'yields_tot_m6p0z3m4_000_20140820_73609.txt']

net_yield_filenames = [f.replace('_tot_', '_net_').replace('20140820_73609', '20160414_155704') for f in total_yield_filenames]

if yield_type == yield_types.absolute:
    yield_filenames = total_yield_filenames
else:
    yield_filenames = net_yield_filenames

model_names = list(map(lambda x: x.split('_')[2], yield_filenames))
initial_masses = [1.5, 2.0, 2.5, 3, 4, 5, 6]

yields = [OrderedDict({}) for k in model_names]
ejecta_masses = [0.0 for k in model_names]
metallicities = [3e-4 for k in model_names]
lifetimes = [
    1.85e9,  # 1.5 Msun
    8.46e8,  # 2.0 Msun
    4.92e8,  # 2.5 Msun
    3.02e8,  # 3.0 Msun
    1.48e8,  # 4.0 Msun
    8.93e7,  # 5.0 Msun
    6.06e7,  # 6.0 Msun
]

for m, filename in enumerate(yield_filenames):
    with open(os.path.join('data/fruity-z0003', filename), 'r') as fyields:
        print("reading " + filename)
        fyields.readline()
        for line in fyields:
            row = line.split()
            elmsymbol = row[0].rstrip('0123456789')

            # decays of unstable nuclides
            if row[0] == 'Zr93':
                elmsymbol = 'Nb'
            if row[0] == 'Tc99':
                elmsymbol = 'Ru'
            if row[0] == 'Pd07':
                elmsymbol = 'Ag'
            if row[0] == 'Cs35':
                elmsymbol = 'Ba'
            if row[0] == 'Cs37':
                elmsymbol = 'Ba'
            if row[0] == 'Pb05':
                elmsymbol = 'Tl'
            if row[0] == 'Bi10':
                elmsymbol = 'Pb'
            if row[0] == 'Po10':
                elmsymbol = 'Pb'

            massnumber = int(row[1])   # (A) the number of nucleons
            massyield = float(row[3])  # absolute stellar yield in solar masses
            ejecta_masses[m] += massyield

            speciesname = '{:}{:d}'.format(elmsymbol.lower(), massnumber)

            if elmsymbol == 'H':
                elmsymbol = 'p'
                speciesname = 'p'

            if speciesname not in yields[m]:
                yields[m][speciesname] = 0.0
            yields[m][speciesname] += massyield

output_yield_file_name = 'yields.txt'
with open(output_yield_file_name, 'w') as fileout:
    print("writing {:} {:}".format(output_yield_file_name, yield_type))
    fileout.write("#test case: FRUITY database models from Straniero et al. 2014 \n")
    fileout.write("#modelname".ljust(25))
    fileout.write("mass".rjust(14))
    fileout.write("Z".rjust(14))
    fileout.write("remnantmass".rjust(14))
    fileout.write("lifetime".rjust(14))
    fileout.write("\n[stellarmodels]\n")
    fileout.write(str(len(model_names)) + "\n")

    for i, model_name in enumerate(model_names):
        fileout.write((chr(ord('A') + i) + ':' + model_name).ljust(25)[:25])
        fileout.write(("{0:6.2f}".format(initial_masses[i])).rjust(14))   # mass
        fileout.write(("{0:6.2e}".format(metallicities[i])).rjust(14))   # metallicity
        fileout.write(("{0:6.3f}".format((initial_masses[i] - ejecta_masses[i]))).rjust(14))   # remnant mass
        fileout.write(("{0:.6e}".format(lifetimes[i])).rjust(14))
        fileout.write("\n")

    fileout.write("\n" + "#species".ljust(8))
    fileout.write("type".rjust(10))
    for m, model_name in enumerate(model_names):
        fileout.write(chr(ord('A') + m).rjust(14))
    fileout.write("\n")

    fileout.write("[yields]\n")

    for species in yields[0]:
        fileout.write(species.ljust(8))

        if yield_type == yield_types.absolute:
            fileout.write("absolute".rjust(10))
        else:
            fileout.write("relative".rjust(10))

        for yields_one_model in yields:
            if species in yields_one_model:
                yieldout = yields_one_model[species]
            else:
                yieldout = 0.0
            fileout.write(("{0:14.6e}".format(yieldout)).rjust(14))

        fileout.write("\n")
