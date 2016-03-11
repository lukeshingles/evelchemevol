#!/usr/bin/env python
import math
import sys
from collections import OrderedDict

# index in this list is the atomic number Z
elements = ['']+[line.split()[1] for line in open('data/atomic_symbols.dat')]

def absoluteYield(i,species):
    return yields[i][species] + (initialMassFrac[i][species] * ejectaMass[i])

# absolute elemental yield
def elYield(i,symbol):
    yieldsum = 0.0
    foundAtLeastOne = False
    for species in yields[i]:
        if species.rstrip('0123456789') == symbol.lower() or (symbol=='h' and species in ['p','d']):
            foundAtLeastOne = True
            yieldsum += yields[i][species] + (initialMassFrac[i][species] * ejectaMass[i])
#    if not foundAtLeastOne:
#        print "ERROR " + symbol + " yield not found in model " + modelname[i]
#        exit()
    return yieldsum

# initial mass elemental fraction for model i
def initElFrac(i,symbol):
    fracsum = 0.0
    for species in yields[i]:
        if species.rstrip('0123456789') == symbol:
            fracsum += initialMassFrac[i][species]
    return fracsum

def ejectaElMassFrac(i,symbol):
    return elYield(i,symbol) / ejectaMass[i]

# solar elemental mass fraction for model i
def solarElMassFrac(symbol):
    fracsum = 0.0
    for species in initcomp['1.4E-02']:
        if species.rstrip('0123456789') == symbol.lower() or (symbol=='h' and species in ['p','d']):
            fracsum += initcomp['1.4E-02'][species]
    return fracsum

def solarBracket(modelnum1,symbol1,modelnum2,symbol2):
    try:
        return math.log10(ejectaElMassFrac(modelnum1,symbol1) / ejectaElMassFrac(modelnum2,symbol2)) - math.log10(solarElMassFrac(symbol1) / solarElMassFrac(symbol2))
    except Exception:
        return -99.0

def limongiToHirschiSpeciesName(strIn):
    split = strIn.lstrip("^").split("^")
    if strIn == "^1^H":
        return "p"
    elif strIn == "^2^H":
        return "d"
    elif strIn == "^3^H":
        return "t"
    else:
        return split[1].lower() + split[0]

# lists begin at element 1 instead of 0
# yields are relative to initial comp and could be negative

hirschimodelcount = 29
limongimodelcount = 8

initz = ('1.4E-02','1.0E-03','1.0E-05','1.0E-07')

files = ['data/hirschi/iniab' + Z + 'As05.gva' for Z in initz]
elementNames = {}

initcomp = {}
for i in range(len(files)):
    initcomp[initz[i]] = {}
    f = open(files[i],'rb')
    for line in f.readlines():
        row = line.split()
        if len(row) != 4:
            continue
        else:
            massfrac = float(row[3])
            if row[1]=='1':
                if row[2]=='1':
                    initcomp[initz[i]]['p'] = massfrac
                    elementNames[1] = 'p'
                else:
                    initcomp[initz[i]]['d'] = massfrac
            else:
                initcomp[initz[i]][row[0] + row[2]] = massfrac
                elementNames[int(row[1])] = row[0]

modelname = ["" for k in range(hirschimodelcount+limongimodelcount+1)]
yields = [OrderedDict([]) for k in range(len(modelname))]

# initial mass fraction
initialMassFrac = [{} for k in range(len(modelname))]
ejectaMass = [0 for k in range(len(modelname))]
remnantMass = [0 for k in range(len(modelname))]

f = open('data/hirschi/yields.tab','rb')
for line in f.readlines():
    row = line.split()
    if len(row)==2:
        modelname[int(row[0])] = row[1][26:37]
    elif len(row)==19:
        ejectaMass[int(row[0])] = float(modelname[int(row[0])][1:4]) - float(row[6])
        remnantMass[int(row[0])] = float(row[6])
    elif len(row)==30:
        for i in range(1,len(row)):
            yields[i][row[0]] = float(row[i])
            if modelname[i][4:7] == 'z14' and row[0] in initcomp['1.4E-02']:
                initialMassFrac[i][row[0]] = initcomp['1.4E-02'][row[0]]
            elif modelname[i][4:7] == 'z01' and row[0] in initcomp['1.0E-03']:
                initialMassFrac[i][row[0]] = initcomp['1.0E-03'][row[0]]
            elif modelname[i][4:7] == 'zm5' and row[0] in initcomp['1.0E-05']:
                initialMassFrac[i][row[0]] = initcomp['1.0E-05'][row[0]]
            elif modelname[i][4:7] == 'zm7' and row[0] in initcomp['1.0E-07']:
                initialMassFrac[i][row[0]] = initcomp['1.0E-07'][row[0]]
            else:
                initialMassFrac[i][row[0]] = 0.0
                #print "Missing initial mass fraction of ", row[0]

f = open('data/limongi2012.txt','rb')
for line in f.readlines():
    row = line.split()
    if len(row) == 9 and row[0] == "#Param":
        for i in range(1,limongimodelcount+1):
            modelname[hirschimodelcount+i] = "Limongi-" + row[i]
    if len(row) == 9 and not line.startswith("#"):
        if row[0] == "M_cut_":
            for i in range(1,limongimodelcount+1):
                ejectaMass[hirschimodelcount+i] = float(modelname[hirschimodelcount+i][9:11]) - float(row[i])
                remnantMass[hirschimodelcount+i] = float(row[i])
        elif not row[0].endswith("_"):
            species = limongiToHirschiSpeciesName(row[0])
            if species == "ni56":
                species = "fe56"
            for i in range(1,limongimodelcount+1):
                if species in yields[i+hirschimodelcount].keys():
                    yields[i+hirschimodelcount][species] += float(row[i])
                else:
                    yields[i+hirschimodelcount][species] = float(row[i])
                initialMassFrac[i+hirschimodelcount][species] = 0.0

modeS0 = 0
modeS4 = 1
modeS5 = 2
modeS4CF88on10 = 3
modeS5CF88on10 = 4

modelabels = ["S0","S4","S5","S4 CF88/10","S5 CF88/10"]

if len(sys.argv) < 3 or sys.argv[1] not in ["z01","zm5"]:
    print "Use command-line options for metallicity (zm5 or z01) and mode from 0 to 4 (" + ", ".join(modelabels)
    sys.exit(0)

metallicity = sys.argv[1]

mode = int(sys.argv[2])
modelabel = modelabels[mode]

print "Metallicity: " + metallicity
print "Mode: " + modelabel

if metallicity == "z01" and modelabel not in ["S0","S4"]:
    print "Invalid mode for z01"
    sys.exit()

modelnamesuffix = ""
if mode == modeS0:
    if metallicity == "z01":
        selectedmodels = [3,9,15,28] #z01S0
        modelname[28] += "*S0/S4"
        print modelname[28]
    else:
        selectedmodels = [5,11,17,29] #zm5S0
        modelname[29] += "*S0/S4"
    modelnamesuffix = ""
elif mode == modeS5:
    selectedmodels = [6,12,18,29] #zm5S4
    modelnamesuffix = "*S5/S4"
elif mode == modeS4:
    if metallicity == "z01":
        selectedmodels = [4,10,16,28] #z01S4
    else:
        selectedmodels = [6,12,18,29] #zm5S4
    #modelnamesuffix = "(S4)"
elif mode == modeS4CF88on10:
    selectedmodels = [6,12,18,29] #zm5S4
    modelnamesuffix = "*CF88on10/CF88"
elif mode == modeS5CF88on10:
    selectedmodels = [6,12,18,29] #zm5S4
    modelnamesuffix = "*S5CF88on10/S4"
else:
    print "Unknown mode: " + str(mode)
    sys.exit(0)

fyields = open('yields.txt', 'w')
fyields.write("#test case: " + modelabel + "\n")
fyields.write("#description".ljust(25))
fyields.write("mass".rjust(14))
fyields.write("Z".rjust(14))
fyields.write("remnantmass".rjust(14))
fyields.write("\n[stellarmodels]\n")
fyields.write(str(len(selectedmodels)) + "\n")

for i in selectedmodels:
    remnantMassThisModel = remnantMass[i]

    if mode == modeS0:
        # the 40 Msun model with zero rotation has to be fudged
        if i == 29: #zm5
            remnantMassThisModel *= remnantMass[17] / remnantMass[18] #S0/S4 factor
        elif i == 28: #z01
            remnantMassThisModel *= remnantMass[15] / remnantMass[16] #S0/S4 factor
    elif mode == modeS4CF88on10:
        remnantMassThisModel *= remnantMass[19] / remnantMass[18] #CF88on10 factor
    elif mode == modeS5:
        remnantMassThisModel *= remnantMass[20] / remnantMass[18] #S5/S4 factor
    elif mode == modeS5CF88on10:
        remnantMassThisModel *= remnantMass[21] / remnantMass[18] #S5CF88on10/S4 factor

    fyields.write((modelname[i] + modelnamesuffix)[:25].ljust(25))
    fyields.write(("%6.2f" % float(modelname[i][1:4])).rjust(14))   # mass
    fyields.write(("%6.2e" % ([1e-7,1e-5,1e-3,1.4e-2][["zm7","zm5","z01","z14"].index(modelname[i][4:7])])).rjust(14))   # metallicity
    fyields.write(("%6.2f" % remnantMassThisModel).rjust(14))   # remnant mass
    fyields.write("\n")

fyields.write("\n" + "#species".ljust(8))
fyields.write("type".rjust(10))
for i in selectedmodels:
    fyields.write( ("model" + str(selectedmodels.index(i)+1)).rjust(14))
fyields.write("\n")

fyields.write("[yields]\n")

for species in yields[selectedmodels[0]]:
    fyields.write(species.ljust(8))

    elcode = species.rstrip('0123456789').title()
    #if elcode in elements and elements.index(elcode) <= 26:
    fyields.write("absolute".rjust(10))
    #else:
    #    fyields.write("relative".rjust(10))

    for i in selectedmodels:
        curyield = absoluteYield(i,species)
        if mode == modeS0:
            # the 40 Msun model with zero rotation has to be fudged
            if i == 29: #zm5
                curyield *= absoluteYield(17,species) / absoluteYield(18,species) #S0/S4 factor
            elif i == 28: #z01
                curyield *= absoluteYield(15,species) / absoluteYield(16,species) #S0/S4 factor
        elif mode == modeS4CF88on10:
            curyield *= absoluteYield(19,species) / absoluteYield(18,species) #CF88on10 factor
        elif mode == modeS5:
            curyield *= absoluteYield(20,species) / absoluteYield(18,species) #S5/S4 factor
        elif mode == modeS5CF88on10:
            curyield *= absoluteYield(21,species) / absoluteYield(18,species) #S5CF88on10/S4 factor

        if metallicity != "z01":
            #override yield with Limongi's yields for Z<=26 (Fe)
            if elcode in elements and elements.index(elcode) <= 26:
                if modelname[i].startswith('G015'):
                    curyield = yields[31][species]
                elif modelname[i].startswith('G020'):
                    curyield = yields[32][species]
                elif modelname[i].startswith('G025'):
                    curyield = yields[33][species]
                elif modelname[i].startswith('G040'):
                    #interpolate between 34 and 50 Msun models
                    curyield = 0.333*yields[35][species] + (1-0.333)*yields[36][species]

        fyields.write(("%14.6e" % curyield).rjust(14))

    fyields.write("\n")
fyields.close()
