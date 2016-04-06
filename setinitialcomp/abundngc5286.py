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

