zfactor = 10 ** -1.66
# M2 r-only group mean abundances (Yong et al. 2014)
# [Fe/H] is about ~-1.68
# logxtofe = log epsilon(X) - log epsilon(Fe)
targetlogxtofe = { 'o': 7.43 - 5.83,
                  'na': 4.63 - 5.83,
                  'mg': 6.30 - 5.83,
                  'al': 5.22 - 5.83,
                  'sr': 0.63 - 5.83,
                   'y': 0.36 - 5.83,
                  'zr': 1.17 - 5.83, # this is Zr II, also need to add Zr I?
                  'ba': 0.69 - 5.83,
                  'la':-0.47 - 5.83,
                  'ce':-0.10 - 5.83,
                  'pr':-0.88 - 5.83,
                  'nd':-0.10 - 5.83,
                  'pb': 0.15 - 5.83
                 }

