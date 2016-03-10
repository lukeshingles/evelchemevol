#!/usr/bin/env python
import matplotlib.pyplot as plt
import math

fig = plt.figure()

ax = fig.add_axes([0.09,  0.08, 0.87, 0.88])

header_row = []

f = open('../dataout/cemodel.txt','rb')
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
                #data[i].append(float(row[i]) * 10 ** -9)
            else:
                data[i].append(max(float(row[i]),10**-30))

for i in range(1,len(data)):
    colors = ['r','g','b','orange','black','purple']
    ax.plot(data[0], data[i], color=colors[i % len(colors)-1], marker='o', markersize=0, lw=3, label=header_row[i])

ax.legend(loc=3,handlelength=1)

fs = 10
#plt.setp(plt.getp(plt.gca(), 'xticklabels'), fontsize=fs)
#plt.setp(plt.getp(plt.gca(), 'yticklabels'), fontsize=fs)

ax.set_xlabel("log10(time [Gyr])", labelpad=8, fontsize=fs)
#ax.set_ylabel('Y', labelpad=8, fontsize=fs)

ax.set_yscale('log')

#ax.set_xlim(xmin=2e3)
ax.set_xlim(xmax=0.4)

ax.set_ylim(ymin=10**-8)
#ax.set_ylim(ymax=2.e6)
#ax.set_ylim(ymax=1e12)
fig.savefig('cemodel.pdf',format='pdf')
plt.close()
