#!/usr/bin/env python
import matplotlib.pyplot as plt
#from matplotlib.lines import Line2D
#import math
import sys
#import numpy as np

for i in range(1):
    fig = plt.figure(figsize=(6, 6))

    ax = fig.add_axes([0.12, 0.12, 0.84, 0.84])
    if i == 0:
        fileSuffix = "Pb-Fe"
        axisy = "[Pb/Fe]"
        ymin = -3.0
        ymax = 3.0
        axisx = "[Fe/H]"
        xmin = -1.0
        xmax = 2.5
    elif i == 1:
        fileSuffix = "La-Fe"
        axisy = "[La/Fe]"
        ymin = -0.5
        ymax = 2.5
        axisx = "[Fe/H]"
        xmin = -0.85
        xmax = 0.8

    xvalues = []
    yvalues = []

    header_row = []

    axisxnum = -1
    axisynum = -1

    f = open('out-abundances.txt','rb')
    for line in f.readlines():
        row = line.split()
        if row[0].startswith("#"):
            header_row = row
            for i, header in enumerate(header_row):
                if header_row[i] == axisx:
                    axisxnum = i
                if header_row[i] == axisy:
                    axisynum = i
                    print i
            if axisxnum == -1 or axisynum == -1:
                print "Columns not found"
                sys.exit(1)
            print header_row[axisynum] + " vs. " + header_row[axisxnum]
        elif axisxnum != -1 and axisynum != -1 and len(row) >= len(header_row):
            xvalues.append(float(row[axisxnum]))
            yvalues.append(float(row[axisynum]))
            #print header_row[axisynum+1]
    #print yvalues[0]
    #for i in range(1,len(header_row)):
    #    colors = ['r','g','b','orange','black','purple']
    ax.plot(xvalues, yvalues, color='black', marker='None', lw=1.5)
    ax.plot(xvalues[0], yvalues[0], color='blue', marker='o', markersize=6,
            markeredgewidth=0,label="Initial")

    #ax.set_ylim(min(filter(lambda x: x>10**-20,y_arr[i]))/1.5,max(y_arr[i])*1.5)

    #ax.legend(loc=3,handlelength=2)

    fs = 10
    #plt.setp(plt.getp(plt.gca(), 'xticklabels'), fontsize=fs)
    #plt.setp(plt.getp(plt.gca(), 'yticklabels'), fontsize=fs)

    ax.set_xlabel(axisx, labelpad=8, fontsize=fs)
    ax.set_ylabel(axisy, labelpad=8, fontsize=fs)

    #ax.set_xlim(xmin,xmax)
    #ax.set_ylim(ymin,ymax)

    fig.savefig('cemodel-' + fileSuffix + '.pdf',format='pdf')
    plt.close()
