#!/usr/bin/env python3
import collections
import math
from solarabundances import solarmetallicity, elsolarlogepsilon

while True:
    userinput = input('Enter helium mass fraction: Y=0.')
    outputheliumfrac = float('0.'+userinput)
    if userinput.isdigit() and 0.0 <= outputheliumfrac <= 100.0:
        break
    print("Invalid input")

#while True:
#    userinput = input('Enter metallicity: Z=0.')
#    outputmetallicity = float('0.'+userinput)
#    if userinput.isdigit() and 0.0 <= outputmetallicity <= 100.0:
#        break
#    print("Invalid input")

elfactor = {}

zfactor = 0.02
targetlogxtofe = {}
#targetlogxtofe is absolute, not relative to solar!

#from initialcompdata.abundngc1851 import zfactor, targetlogxtofe
#from initialcompdata.abundngc2808 import zfactor, targetlogxtofe
from initialcompdata.abundm2 import zfactor, targetlogxtofe
#from initialcompdata.abundngc5286 import zfactor, targetlogxtofe
#from initialcompdata.abundomegacen import zfactor, targetlogxtofe

#calculate zfactor based on outputmetallicity
#zfactor = outputmetallicity / solarmetallicity
#or calculate outputmetallicity based on zfactor
outputmetallicity = zfactor * solarmetallicity

outputhydrogenfrac = (1.0 - outputheliumfrac - outputmetallicity)

# read the reference isotopic composition
refelnumberfrac = collections.OrderedDict()
with open('initial_comp_solar_a4.dat','r') as infile:
    for line in infile:
        row = [line[0:7], line[7:21], line[21:28]]
        elcode = row[0].strip(' 0123456789')
        if row[0].strip(' ') in ['p','d']:
            elcode = 'h'

        if elcode not in refelnumberfrac:
            refelnumberfrac[elcode] = 0.0
        refelnumberfrac[elcode] += float(row[1])

print('reference: [Fe/H]={0:.4f}'.format(math.log10(refelnumberfrac['fe']/refelnumberfrac['h'])+12-elsolarlogepsilon['fe']))

print('output: [Fe/H]={0:.4f}'.format(math.log10(zfactor*refelnumberfrac['fe']/outputhydrogenfrac)+12-elsolarlogepsilon['fe']))

# iterate over the the elemental number fractions to get scaling factors
for elcode in refelnumberfrac:
    if elcode in elsolarlogepsilon:
        # these relate to the reference isotopic composition
        logxtofe = math.log10(refelnumberfrac[elcode]/refelnumberfrac['fe'])
        logxtofesolar = (elsolarlogepsilon[elcode]-elsolarlogepsilon['fe'])
        logxtoferelsolar = logxtofe - logxtofesolar

        consolelineout = '[{0:s}/Fe]={1:.4f}'.format(elcode.title(), logxtoferelsolar)

        if elcode in targetlogxtofe:
            elfactor[elcode] = 10 ** (targetlogxtofe[elcode] - logxtofe)
            consolelineout += "\tNEW VALUE: {0:.4f}".format(targetlogxtofe[elcode]-logxtofesolar)

        if elcode == 'h':
            print('X={0:.3f}'.format(outputhydrogenfrac))
#        elif elcode == 'he':
#            print('Y={0:.3f}'.format(outputheliumfrac))
        else:
            print(consolelineout)

# read the reference composition and produce the output by
# applying the scaling factors
with open('initial_comp_solar_a0.dat','r') as infile:
    with open('initial_comp.dat','w') as outfile:
        for line in infile:
            row = [line[0:7], line[7:21], line[21:28]]
            elcode = row[0].strip(' 0123456789')
            if row[0].strip(' ') == 'p':
                row[1] = '{0:14.4E}'.format(outputhydrogenfrac)
            elif row[0].strip(' ') == 'he4':
                row[1] = '{0:14.4E}'.format(outputheliumfrac/4.0)
            elif elcode in elfactor:
                row[1] = '{0:14.4E}'.format(float(row[1])*elfactor[elcode]*zfactor)
            else:
                row[1] = '{0:14.4E}'.format(float(row[1])*zfactor)

            outfile.write("".join(row) + "\n")
#            if elcode not in targetxtofe:
#                refelnumberfrac[elcode] = 0.0
#            refelnumberfrac[elcode] += float(row[1])
