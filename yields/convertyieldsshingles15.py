#!/usr/bin/env python3
from collections import OrderedDict
import struct

elements = ('','H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu',
        'Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La',
        'Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac')

#(yval, k14alpha, s15alpha) = ('24', '4', '0')
(yval, k14alpha, s15alpha) = ('40', '0', '0')

modelNames = [s.format(yval,k14alpha,s15alpha) for s in ['m1.7z0006y{0}a{1}pmz001s320','m2.36z0006y{0}a{1}pmz001s320','m3z0006y{0}a{2}pmz001s320','m4z0006y{0}a{2}s320','m5z0006y{0}a{2}s320','m6z0006y{0}a{2}s320']]
yields = [OrderedDict({}) for k in range(len(modelNames))]
ejectaMass = [0.0 for k in range(len(modelNames))]
initMass = [float(modelNames[k][1:].split('z')[0]) for k in range(len(modelNames))]
lifetime = [0.0 for k in range(len(modelNames))]
metallicity = [0.0 for k in range(len(modelNames))]

customlifetimes = {'m1.7z0006y24': 1.4355016E+09, 'm1.7z0006y40': 5.6321442E+08, 'm2.36z0006y24': 5.4127692E+08, 'm2.36z0006y40': 2.4111768E+08}

for m in range(len(modelNames)):
    evmodelname = modelNames[m].split('a')[0]
    metallicity[m] = 0.0006
    if evmodelname in customlifetimes:
        lifetime[m] = customlifetimes[evmodelname]
    else:
        with open("/Users/lukes/Dropbox/Papers (first author)/2015 He-enhanced IMAGB Stars/generateplots/evolution/" + evmodelname + "/evall.dat", mode='rb') as evfile: #or did you mean evall.dat?
            fileContent = evfile.read()

            bytepos = range(20,len(fileContent)-16,268)[-1]
            lifetime[m] = struct.unpack("<d", fileContent[bytepos+4:bytepos+12])[0]

    with open('/Users/lukes/Dropbox/Papers (first author)/2015 He-enhanced IMAGB Stars/generateplots/ppns320species/' + modelNames[m] + '/yields.txt','r') as fyields:
        print("reading " + modelNames[m])
        for line in fyields:
            if not line.startswith("#"):
                row = line.split()
                elcode = row[0].rstrip('0123456789').title()
                if elcode in elements and elements.index(elcode) <= 26:
                    yields[m][row[0]] = float(row[3]) #absolute yield
                else:
                    yields[m][row[0]] = float(row[2]) #relative yield
                ejectaMass[m] += float(row[3])

with open('data/kobayashi06snyields.txt','r') as k06yields:
    print("reading kobayashi06snyields.txt")

    k06startindex = len(modelNames)
    k06masslist = [13,15,18,20,25,30,40]
    for mass in k06masslist:
        modelNames.append('k06-m{0:d}'.format(mass))
        initMass.append(float(mass))
        lifetime.append(lifetime[k06startindex-1] * ((mass/initMass[k06startindex-1]) ** (-3.0)))
        metallicity.append(0.001)
        yields.append(OrderedDict({}))
    
    for line in k06yields:
        if line.startswith("0.001"):
            row = line.split()
            if row[1] == "M_cut_":
                for i in range(len(k06masslist)):
                    ejectaMass.append(k06masslist[i] - float(row[i+2]))
            elif row[1] != "M_final_":
                
                if row[1] in ['p','d']:
                    speciesname = row[1]
                else:
                    # turn ^26^Mg into mg26
                    speciesname = "".join(reversed(row[1][1:].lower().split('^')))
                for i in range(len(k06masslist)):
                    yields[k06startindex+i][speciesname] = float(row[2+i]) #absolute yield

with open('yields.txt', 'w') as fileout:
    print("writing yields.txt")
    fileout.write("#test case: AGB models\n")
    fileout.write("#modelname".ljust(25))
    fileout.write("mass".rjust(14))
    fileout.write("Z".rjust(14))
    fileout.write("remnantmass".rjust(14))
    fileout.write("lifetime".rjust(14))
    fileout.write("\n[stellarmodels]\n")
    fileout.write(str(len(modelNames)) + "\n")

    for i in range(len(modelNames)):
        fileout.write((chr(ord('A') + i) + ':' + modelNames[i]).ljust(25)[:25])
        fileout.write(("%6.2f" % initMass[i]).rjust(14))   # mass
        fileout.write(("%6.2e" % metallicity[i]).rjust(14))   # metallicity
        fileout.write(("%6.3f" % (initMass[i]  - ejectaMass[i])).rjust(14))   # remnant mass
        fileout.write(("%.6e" % lifetime[i]).rjust(14))
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
        if elcode in elements and elements.index(elcode) <= 26:
            fileout.write("absolute".rjust(10))
        else:
            fileout.write("relative".rjust(10))

        for i in range(len(yields)):
            if species in yields[i]:
                yieldout = yields[i][species]
            else:
                yieldout = 0.0
                if i == k06startindex:
                    print("Warning '" + species + "' undefined for Kobayashi+2006 models")
            fileout.write(("%14.6e" % yieldout).rjust(14))
 
        fileout.write("\n")