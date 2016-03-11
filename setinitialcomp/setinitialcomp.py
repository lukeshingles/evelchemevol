#!/usr/bin/env python3
import collections
import math

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

zfactor = 0.02
targetlogxtofe = {}
#targetlogxtofe is absolute, not relative to solar!

if False:
    from ngc2808 import *

if True:
    zfactor = 10 ** -1.66
    # M2 r-only group mean abundances (Yong et al. 2014)
    # [Fe/H] is about ~-1.68
    # logxtofe = log epsilon(X) - log epsilon(Fe)
    targetlogxtofe = { 'o': 7.43 - 5.83,
                      'na': 4.63 - 5.83,
                      'mg': 6.30 - 5.83,
                      'al': 5.22 - 5.83,
                      'sr': 0.63 - 5.83,
                       'y': 0.36 - 5.83,
                      'zr': 1.17 - 5.83, # this is Zr II, also need to add Zr I?
                      'ba': 0.69 - 5.83,
                      'la':-0.47 - 5.83,
                      'ce':-0.10 - 5.83,
                      'pr':-0.88 - 5.83,
                      'nd':-0.10 - 5.83,
                      'pb': 0.15 - 5.83
                     }

if False:
    zfactor = 10 ** -1.92
    # mean of s-poor population in NGC5286
    # from Marino et al. (2015) 2015MNRAS.450..815M
    # [Fe/H] = -1.92
    # log X/Fe = [X/Fe] + log(X/Fe)_solar
    targetlogxtofe = { 'o': 0.58 + elsolarlogepsilon['o']  - elsolarlogepsilon['fe'],
                      'na': 0.18 + elsolarlogepsilon['na'] - elsolarlogepsilon['fe'],
                       'y': -0.04 +elsolarlogepsilon['y']  - elsolarlogepsilon['fe'],
                      'zr': 0.17 + elsolarlogepsilon['zr'] - elsolarlogepsilon['fe'],
                      'ba': 0.03 + elsolarlogepsilon['ba'] - elsolarlogepsilon['fe'],
                      'la': 0.29 + elsolarlogepsilon['la'] - elsolarlogepsilon['fe'],
                      'ce': 0.24 + elsolarlogepsilon['ce'] - elsolarlogepsilon['fe'],
                      'pr': 0.38 + elsolarlogepsilon['pr'] - elsolarlogepsilon['fe'],
                      'nd': 0.20 + elsolarlogepsilon['nd'] - elsolarlogepsilon['fe']
                      }

if False:
    zfactor = 10 ** -1.92
    # Smith et al. (2000) ROA 219 in Omega Centauri
    # [Fe/H] is about ~-1.7
    #logxtofe = log epsilon(X) - log epsilon(Fe)
    targetlogxtofe = {'rb': 1.34 - 6.25,
                       'y': 1.15 - 6.25,
                      'zr': 2.01 - 6.25,
                      'ba': 1.88 - 6.25,
                      'la': 0.75 - 6.25,
                      'ce': 0.42 - 6.25,
                      'pb': 0.40 + elsolarlogepsilon['pb'] - elsolarlogepsilon['fe'] #D'Orazi+2011 Leiden 60066
                      }


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
