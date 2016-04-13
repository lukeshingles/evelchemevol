#!/usr/bin/env python3
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os

# Asplund et al. 2009 ARAA (CHECK AGAINST PAPER)
solarabund = {'pm': -9.0, 'tc': -9.0, 'er': 0.92, 'rh': 1.06, 'gd': 1.05,
    'pb': 2.04, 'te': 2.18, 'sm': 0.94, 'au': 0.8, 'po': -9.0, 'pt': 1.62,
    'sb': 1.01, 'lu': 0.09, 'ru': 1.76, 'ag': 1.2, 'ga': 3.08, 'ne': 7.97,
    'ce': 1.58, 's': 7.15, 'mo': 1.94, 'zr': 2.53, 'kr': 3.25, 'na': 6.24,
    'cr': 5.64, 'nd': 1.45, 'w': 0.65, 'si': 7.55, 'bi': 0.65, 'eu': 0.51,
    'co': 4.87, 'la': 1.17, 'al': 6.45, 'dy': 1.13, 'ca': 6.29, 'k': 5.08,
    'y': 2.21, 'b': 2.79, 'yb': 0.92, 'he': 10.93, 'hf': 0.71, 'ti': 4.91,
    'tb': 0.32, 'mn': 5.48, 're': 0.26, 'ar': 6.44, 'se': 3.34, 'br': 2.54,
    'os': 1.35, 'ta': -0.12, 'ir': 1.32, 'f': 4.42, 'ho': 0.47, 'c': 8.47,
    'sc': 3.05, 'cd': 1.71, 'li': 3.26, 'rb': 2.52, 'n': 7.87, 'cs': 1.08,
    'pd': 1.65, 'i': 1.55, 'as': 2.3, 'tl': 0.77, 'sn': 2.07, 'ni': 6.2,
    'ba': 2.18, 'o': 8.73, 'cu': 4.25, 'xe': 2.24, 'hg': 1.17, 'v': 3.96,
    'p': 5.41, 'sr': 2.88, 'ge': 3.58, 'tm': 0.12, 'fe': 7.54, 'mg': 7.64,
    'in': 0.76, 'be': 1.3, 'pr': 0.76, 'cl': 5.23, 'nb': 1.41, 'zn': 4.63}

fig, ax = plt.subplots(1, figsize=(5,3.5), sharex=True, tight_layout={"pad":0.3})

header_row = []
columnlist =['[Fe/H]','O','Na','Rb','Y','Zr','Ba','La','Ce','Nd']
with open('out-abundances.txt','r') as f:
    for line in f:
        row = line.split()
        if row[0].startswith("#"):
            header_row = row
            print(header_row)
            data = {x:[] for x in header_row}
        elif len(header_row) > 1 and len(row) >= len(header_row) and float(row[1]) > 0.0:
            for i, header in enumerate(header_row):
                if i == 1:
                    # don't plot every data point
                    #if len(data[0]) > 0 and math.log10(float(row[0]) * 10 ** -6) < data[0][-1] + 0.01:
                    #    break
                    data[header_row[i]].append(float(row[i]))
                else:
                    data[header_row[i]].append(float(row[i]))

colors = ['black','r','g','blue','purple','orange']
linestyle = ['-','--','-.',':']
for i, column in enumerate(columnlist):
    if column.startswith('['):
        yvalues = data[column]
        label = column
    else:
        yvalues = [X - Fe - (solarabund[column.lower()] - solarabund['fe']) for (X,Fe) in zip(data[column],data['Fe'])]
        label = column
#        label = '[{0}/Fe]'.format(column)

    ax.plot(data['time'], yvalues, color=colors[(i) % 6], marker='o',
            markersize=0, lw=1, linestyle=linestyle[int((i)/len(colors)) % 4],
            label=label)

#ax.set_ylim(min(filter(lambda x: x>10**-20,y_arr[i]))/1.5,max(y_arr[i])*1.5)

fs = 9

ax.set_title(os.getcwd().split('/')[-1])
ax.legend(loc='best',handlelength=2.5,frameon=False,ncol=2,prop={'size':fs})

ax.yaxis.set_minor_locator(ticker.MultipleLocator(base=0.1))
ax.yaxis.set_major_locator(ticker.MultipleLocator(base=0.5))

plt.setp(plt.getp(plt.gca(), 'xticklabels'), fontsize=fs)
plt.setp(plt.getp(plt.gca(), 'yticklabels'), fontsize=fs)
#ax.set_aspect(0.4)

ax.set_xlabel("time [years]", labelpad=8, fontsize=fs)
ax.set_ylabel('[X/Fe]', labelpad=8, fontsize=fs)

ax.set_xlim(xmin=1e7)
ax.set_xscale('log')

fig.savefig('cemodel-abund.pdf',format='pdf')
plt.close()
