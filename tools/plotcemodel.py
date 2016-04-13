#!/usr/bin/env python3
import matplotlib.pyplot as plt
import os

fig, ax = plt.subplots(1, figsize=(5,3.5), sharex=True, tight_layout={"pad":0.3})

header_row = []

with open('out-cemodel.txt','r') as f:
    for line in f:
        row = line.split()
        if row[0].startswith("#"):
            header_row = row
            print(header_row)
            data = {x:[] for x in header_row}
        elif len(header_row) > 1 and len(row) >= len(header_row):
            if float(row[1]) != 0.0:
                for i, header in enumerate(header_row):
                    if i == 1:
                        #if len(data[0]) > 0 and float(row[0]) != 0.0:
                            #if math.log10(float(row[i]) * 10 ** -6) < data[0][-1] + 0.01:
                            #    break
                        data[header_row[i]].append(float(row[i]))
                    else:
                        if i > 0 and not "E" in row[i]:
                            row[header_row[i]] = 0.00
                        if header_row[i].startswith("["):
                            data[header_row[i]].append(10 ** float(row[i]))
                        else:
                            data[header_row[i]].append(max(float(row[i]),10**-30))

for i in range(2,len(data)):
    column = header_row[i]
    if header_row[i] == "[Fe/H]":
        header_row[i] = "10^[Fe/H]"
    colors = ['r','g','b','orange','black','purple']
    ax.plot(data['time'], data[column], color=colors[(i-1) % len(colors)-1],
            marker='o', markersize=0, lw=1.5, label=header_row[i])

fs = 9
ax.legend(loc='best',handlelength=2.5,frameon=False,ncol=1,prop={'size':fs-1})

plt.setp(plt.getp(plt.gca(), 'xticklabels'), fontsize=fs)
plt.setp(plt.getp(plt.gca(), 'yticklabels'), fontsize=fs)

ax.set_xlabel("Time [years]", labelpad=8, fontsize=fs)

ax.set_xscale('log')
ax.set_yscale('log')

ax.set_title(os.getcwd().split('/')[-1])

ax.set_ylim(ymin=10**-8)
fig.savefig('cemodel.pdf',format='pdf')
plt.close()
