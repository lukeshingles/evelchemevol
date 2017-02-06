#!/usr/bin/env python3
from collections import OrderedDict
from enum import Enum

elements = ['']+[line.split()[1] for line in open('data/atomic_symbols.dat')]

lowest_initial_mass_with_c13_pocket = 3.00

yield_types = Enum('yield_types', 'relative absolute')

yield_type = yield_types.relative

yield_files = (
    'yields-m1.0z001.pmz2e-3.dat',
    'yields-m1.25z001.pmz2e-3.dat',
    'yields-m1.5z001.pmz2e-3.dat',
    'yields-m2.0z001.pmz2e-3.dat',
    'yields-m2.25z001.pmz2e-3.dat',
    'yields-m2.5z001.pmz2e-3.dat',
    'yields-m2.75z001.pmz2e-3.dat',
    'yields-m3.0z001.pmz2e-3.dat')

if lowest_initial_mass_with_c13_pocket >= 3.25:
    yield_files += ('yields-m3.25z001.pmz1e-3.dat',)
else:
    yield_files += ('yields-m3.25z001.dat',)

if lowest_initial_mass_with_c13_pocket >= 3.5:
    yield_files += ('yields-m3.5z001.pmz1e-3.dat',)
else:
    yield_files += ('yields-m3.5z001.dat',)

yield_files += (
    'yields-m4.0z001.dat',
    'yields-m4.5z001.dat',
    'yields-m5.0z001.dat',
    'yields-m6.0z001.dat',
    'yields-m7.0z001.dat')

model_names = [x[7:-4] for x in yield_files]

yields = [OrderedDict({}) for k in model_names]
netyields = [OrderedDict({}) for k in model_names]
ejecta_masses = [0.0 for k in model_names]
initial_masses = [float(model_name[1:].split('z')[0]) for model_name in model_names]
lifetimes = [0.0 for k in model_names]

with open('data/fishlock-z001-20140503/lifetimes.dat', 'r') as flifetimes:
    for line in flifetimes:
        linesplit = line.split()
        for m, init_mass_value in enumerate(initial_masses):
            if "{0:1.2f}".format(init_mass_value) == linesplit[0]:
                lifetimes[m] = float(linesplit[1])

for m, agb_yield_file in enumerate(yield_files):
    with open('data/fishlock-z001-20140503/' + agb_yield_file, 'r') as fyields:
        print("reading " + agb_yield_file)
        for line in fyields:
            if not line.startswith("#") and not line.startswith("  species"):
                linesplit = line.split()
                netyields[m][linesplit[0]] = float(linesplit[2])  # net/relative yield
                yields[m][linesplit[0]] = float(linesplit[3])  # absolute yield
                ejecta_masses[m] += float(linesplit[3])

output_yield_file_name = 'yields.txt'
with open(output_yield_file_name, 'w') as fileout:
    print("writing {:} {:}".format(output_yield_file_name, yield_type))
    fileout.write("#AGB yields from Fishlock et al. (2014)\n")
    fileout.write("#modelname".ljust(25))
    fileout.write("mass".rjust(14))
    fileout.write("Z".rjust(14))
    fileout.write("remnantmass".rjust(14))
    fileout.write("lifetime".rjust(14))
    fileout.write("\n[stellarmodels]\n")
    fileout.write(str(len(model_names)) + "\n")

    for i, model_name in enumerate(model_names):
        fileout.write((chr(ord('A') + i) + ':' + model_name).ljust(25)[:25])
        fileout.write(("{0:6.2f}".format(initial_masses[i])).rjust(14))  # mass
        fileout.write(("{0:6.2e}".format(0.001)).rjust(14))  # metallicity
        fileout.write(("{0:6.3f}".format((initial_masses[i] - ejecta_masses[i]))).rjust(14))  # remnant mass
        fileout.write(("{0:.6e}".format(lifetimes[i])).rjust(14))
        fileout.write("\n")

    fileout.write("\n" + "#species".ljust(8))
    fileout.write("type".rjust(10))
    for m, model_name in enumerate(model_names):
        # fileout.write(model_name[:13].rjust(14))
        fileout.write(chr(ord('A') + m).rjust(14))
    fileout.write("\n")

    fileout.write("[yields]\n")

    for species in yields[0]:
        fileout.write(species.ljust(8))

        # elcode = species.rstrip('0123456789').title()
        # if elcode in elements and elements.index(elcode) <= 26: #is this a good idea?
        if yield_type == yield_types.absolute:
            fileout.write("absolute".rjust(10))
            for yields_one_model in yields:
                fileout.write(("{0:14.6e}".format(yields_one_model[species])).rjust(14))
        else:
            fileout.write("relative".rjust(10))
            for yields_one_model in netyields:
                fileout.write(("{0:14.6e}".format(yields_one_model[species])).rjust(14))

        fileout.write("\n")
