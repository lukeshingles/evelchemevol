#!/usr/bin/env python
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import math
import sys

fig = plt.figure(figsize=(6, 6))

ax = fig.add_axes([0.08,  0.08, 0.88, 0.88])

header_row = []

f = open('../dataout/cemodel-abund.txt','rb')
for line in f.readlines():
    row = line.split()
    if row[0].startswith("#"):
        header_row = row
        print header_row
        data = [[] for x in range(len(header_row))]
    elif len(header_row) > 1 and len(row) >= len(header_row):
        for i in range(len(header_row)):
            if i == 0:
                if len(data[0]) > 0 and math.log10(float(row[i]) * 10 ** -9) < data[0][-1] + 0.01:
                    break
                data[i].append(math.log10(float(row[i]) * 10 ** -9))
            else:
                data[i].append(float(row[i]))

for i in range(1,len(header_row)):
    colors = ['black','r','g','blue','purple','orange']
    linestyle = ['-','--','-.',':']
    ax.plot(data[0], data[i], color=colors[(i-1) % 6], marker='o', markersize=0, lw=2, linestyle=linestyle[int((i-1)/len(colors)) % 4], label=header_row[i])

#ax.set_ylim(min(filter(lambda x: x>10**-20,y_arr[i]))/1.5,max(y_arr[i])*1.5)

ax.legend(loc=4,handlelength=2)

fs = 10
#plt.setp(plt.getp(plt.gca(), 'xticklabels'), fontsize=fs)
#plt.setp(plt.getp(plt.gca(), 'yticklabels'), fontsize=fs)
#ax.set_aspect(0.4)

ax.set_xlabel("log10(time [Gyr])", labelpad=8, fontsize=fs)
#ax.set_ylabel('Y', labelpad=8, fontsize=fs)

#ax.set_xlim(xmin=-1.8)
#ax.set_xlim(xmax=0.3)

#ax.set_xlim(-1.41,0.0)
#ax.set_ylim(-1,2.5)
#ax.set_ylim(ymax=2.e6)
#ax.set_ylim(ymax=1e12)
fig.savefig('cemodel-abund.pdf',format='pdf')
plt.close()
