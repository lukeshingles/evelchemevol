#!/usr/bin/env python
import math
import sys
import numpy as np

elements = ('','H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu',
        'Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La',
        'Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac')

# absolute elemental yield in solar masses  
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

# number ratio
def elNumbertoFe(i,symbol,femodelnum):
    yieldsum = 0.0
    foundAtLeastOne = False
    for species in yields[i]:
        if species.rstrip('0123456789') == symbol.lower() or (symbol=='h' and species in ['p','d']):
            foundAtLeastOne = True
            yieldsum += (yields[i][species] + (initialMassFrac[i][species] * ejectaMass[i])) / int(species.lstrip('abcdefghijklmnopqrstuvwxyz'))
    return yieldsum / (elYield(femodelnum,'fe') / 56.0)

def logEpsilon(i,symbol):
    yieldsum = 0.0
    foundAtLeastOne = False
    for species in yields[i]:
        if species.rstrip('0123456789') == symbol.lower() or (symbol=='h' and species in ['p','d']):
            foundAtLeastOne = True
            yieldsum += (yields[i][species] + (initialMassFrac[i][species] * ejectaMass[i])) / int(species.lstrip('abcdefghijklmnopqrstuvwxyz'))
    return math.log10(yieldsum / (elYield(femodelnum,'h'))) + 12.0

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

def productionFactor(i,symbol):
    return (elYield(i,symbol) / ejectaMass[i]) / initElFrac(i,symbol)

def massFracRelSolar(i,symbol):
    return (elYield(i,symbol) / ejectaMass[i]) / solarElMassFrac(symbol)

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
yields = [{} for k in range(len(modelname))]

# initial mass fraction
initialMassFrac = [{} for k in range(len(modelname))]
ejectaMass = [0 for k in range(len(modelname))]

f = open('data/hirschi/yields.tab','rb')
for line in f.readlines():
    row = line.split()
    if len(row)==2:
        modelname[int(row[0])] = row[1][26:37]
    elif len(row)==19:
        ejectaMass[int(row[0])] = float(modelname[int(row[0])][1:4]) - float(row[6])
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


print "model\t\t[Pb/Eu]\t[Y/Fe]\t[Ba/Fe]\t[Fe/H]\t[Y/Ba]"
for i in range(1,hirschimodelcount+1):
    #if not modelname[i].startswith('G025'):
    #    continue
    #if modelname[i][4:7] != 'zm5':
    #    continue
    
    if modelname[i].startswith('G015'):
        femodelnum = 31
    elif modelname[i].startswith('G020'):
        femodelnum = 32
    elif modelname[i].startswith('G025'):
        femodelnum = 33
    elif modelname[i].startswith('G040'):
        femodelnum = 36
    else:
        femodelnum = i
    #femodelnum = i

    print i,modelname[i] + '\t%0.2f' % solarBracket(i,'pb',i,'eu') + '\t%0.2f' % solarBracket(i,'y',femodelnum,'fe') + '\t%0.2f' % solarBracket(i,'ba',femodelnum,'fe') + '\t%0.2f' % solarBracket(femodelnum,'fe',femodelnum,'h') + '\t%0.2f' % solarBracket(i,'y',i,'ba')

relabel = {}
relabel = {'G025zm5S003': '$v_{rot}=0.0$', 'G025zm5S413': '$v_{rot}=0.4$', 'G025zm5S513': '$v_{rot}=0.5$', 'G025zm5S504': '$v_{rot}=0.5$ O17agCF88/10'}

colorList = ('black','blue','red','green')
linestyleList = ['-','--','-.',':']
markerList = ['s','D','o','^']
colorNumber = 0
for i in range(1,len(modelname)):
    if modelname[i].startswith('G025zm5') and not modelname[i] == 'G025zm5S404': #and not modelname[i].endswith('04'):
        femodelnum = 33
        listX = []
        listY = []
        for Z in range(37,83):
            XtoFe = solarBracket(i,elements[Z],femodelnum,'fe')
            #XtoFe = elNumbertoFe(i,elements[Z],femodelnum)
            
            #print "[" + elements[Z] + "/Fe]", XtoFe
            if XtoFe != -99 and XtoFe != 0:
                listX.append(Z)
                listY.append(XtoFe)
        if modelname[i] in relabel:
            linelabel = relabel[modelname[i]]
        else:
            linelabel = modelname[i]
        
        #print "#" + linelabel
        #print "#Z, [Fe/H]"
        #print "\n".join(map(lambda x,y:str(x)+",%0.2f"%y,listX, listY))

