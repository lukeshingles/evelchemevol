from abundsolar import elsolarlogepsilon

zfactor = 10 ** -1.92
# Smith et al. (2000) ROA 219 in Omega Centauri
# [Fe/H] is about ~-1.7
#logxtofe = log epsilon(X) - log epsilon(Fe)
targetlogxtofe = {'rb': 1.34 - 6.25,
                   'y': 1.15 - 6.25,
                  'zr': 2.01 - 6.25,
                  'ba': 1.88 - 6.25,
                  'la': 0.75 - 6.25,
                  'ce': 0.42 - 6.25,
                  'pb': 0.40 + elsolarlogepsilon['pb'] - elsolarlogepsilon['fe'] #D'Orazi+2011 Leiden 60066
                  }
