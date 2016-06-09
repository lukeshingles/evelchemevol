#!/usr/bin/env python3
from collections import OrderedDict
#import struct
import sys

elements = ['']+[line.split()[1] for line in open('data/atomic_symbols.dat')]

#(yval, k14alpha, s15alpha) = ('25', '4', '0') #Helium normal
(yval, k14alpha, s15alpha) = ('40', '0', '0') #Helium rich

model_names = [s.format(yval,k14alpha,s15alpha) for s in [
                        'm1.7z0006y{0}a{1}pmz001s320',
                        'm2.36z0006y{0}a{1}pmz001s320',
                        'm3z0006y{0}a{2}pmz001s320',
                        'm4z0006y{0}a{2}s320',
                        'm5z0006y{0}a{2}s320',
                        'm6z0006y{0}a{2}s320']]
yields = [OrderedDict({}) for k in model_names]
ejecta_masses = [0.0 for k in model_names]
initial_masses = [float(model_name[1:].split('z')[0]) for model_name in model_names]
lifetimes = [0.0 for k in model_names]
metallicities = [0.0 for k in model_names]

customlifetimes = {'m1.7z0006y24':  1.4355016E+09,
                   'm2.36z0006y24': 5.4127692E+08,
                   'm3z0006y24':    2.9003497E+08,
                   'm4z0006y24':    1.4625349E+08,
                   'm5z0006y24':    8.9311182E+07,
                   'm6z0006y24':    6.1418488E+07,
                   'm1.7z0006y40':  5.6321442E+08,
                   'm2.36z0006y40': 2.4111768E+08,
                   'm3z0006y40':    1.3478728E+08,
                   'm4z0006y40':    7.0936165E+07,
                   'm5z0006y40':    4.4938030E+07,
                   'm6z0006y40':    3.1594309E+07}

for m,modelname in enumerate(model_names):
    evmodelname = model_names[m].split('a')[0]
    metallicities[m] = 0.0006
    if evmodelname in customlifetimes:
        lifetimes[m] = customlifetimes[evmodelname]
    else:
        print("Fatal error: Lifetime needed for model '{:}'".format(model_names[m]))
        sys.exit()

        # with open("/Users/lukes/Dropbox/Papers (first author)/2015 He-enhanced IMAGB Stars/generateplots/evolution/"
        #           + evmodelname + "/evall.dat", mode='rb') as evfile: #or did you mean evall.dat?
        #     fileContent = evfile.read()
        #
        #     bytepos = range(20,len(fileContent)-16,268)[-1]
        #     lifetimes[m] = struct.unpack("<d", fileContent[bytepos+4:bytepos+12])[0]
        #     #print("                   '{:}':    {:13.7E},".format(evmodelname,lifetimes[m]))

    with open('data/shinglesetal2015/' + modelname + '/yields.txt','r') as fyields:
        print("reading " + modelname)
        for line in fyields:
            if not line.startswith("#"):
                row = line.split()
                elcode = row[0].rstrip('0123456789').title()
                if elcode in elements and elements.index(elcode) <= 26:
                    yields[m][row[0]] = float(row[3]) #absolute yield
                else:
                    yields[m][row[0]] = float(row[2]) #relative yield
                ejecta_masses[m] += float(row[3])

with open('data/kobayashi06snyields.txt','r') as k06yields:
    print("reading kobayashi06snyields.txt")

    k06startindex = len(model_names)
    k06masslist = [13,15,18,20,25,30,40]
    for mass in k06masslist:
        model_names.append('k06-m{0:d}'.format(mass))
        initial_masses.append(float(mass))
        lifetimes.append(lifetimes[k06startindex-1] * ((mass/initial_masses[k06startindex-1]) ** (-3.0)))
        metallicities.append(0.001)
        yields.append(OrderedDict({}))

    for line in k06yields:
        if line.startswith("0.001"):
            row = line.split()
            if row[1] == "M_cut_":
                for i,k06mass in enumerate(k06masslist):
                    ejecta_masses.append(k06mass - float(row[i+2]))
            elif row[1] != "M_final_":
                if row[1] in ['p','d']:
                    speciesname = row[1]
                else:
                    # turn ^26^Mg into mg26
                    speciesname = "".join(reversed(row[1][1:].lower().split('^')))
                for i,k06mass in enumerate(k06masslist):
                    yields[k06startindex+i][speciesname] = float(row[2+i]) #absolute yield

with open('yields.txt', 'w') as fileout:
    print("writing yields.txt")
    fileout.write("#test case: Shingles et al. 2015 AGB models and Kobayashi 2006 massive star yields\n")
    fileout.write("#modelname".ljust(25))
    fileout.write("mass".rjust(14))
    fileout.write("Z".rjust(14))
    fileout.write("remnantmass".rjust(14))
    fileout.write("lifetime".rjust(14))
    fileout.write("\n[stellarmodels]\n")
    fileout.write(str(len(model_names)) + "\n")

    for i,model_name in enumerate(model_names):
        fileout.write((chr(ord('A') + i) + ':' + model_name).ljust(25)[:25])
        fileout.write(("{0:6.2f}".format(initial_masses[i])).rjust(14))   # mass
        fileout.write(("{0:6.2e}".format(metallicities[i])).rjust(14))   # metallicity
        fileout.write(("{0:6.3f}".format((initial_masses[i]  - ejecta_masses[i]))).rjust(14))   # remnant mass
        fileout.write(("{0:.6e}".format(lifetimes[i])).rjust(14))
        fileout.write("\n")

    fileout.write("\n" + "#species".ljust(8))
    fileout.write("type".rjust(10))
    for m,model_name in enumerate(model_names):
#        shortname = (model_names[m].split('z')[0] + 'y' + model_names[m].split('y')[1])
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

        for i,yields_thismodel in enumerate(yields):
            if species in yields_thismodel:
                yieldout = yields_thismodel[species]
            else:
                yieldout = 0.0
                if i == k06startindex:
                    print("Warning: '" + species + "' undefined for Kobayashi+2006 models")
            fileout.write(("{0:14.6e}".format(yieldout)).rjust(14))

        fileout.write("\n")
