import h5py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('mplstyle.txt')

# have to open the input file with h5py (g4 doesn't write pandas-ready hdf5)
g4sfile = h5py.File('../g4simple.hdf5', 'r')
g4sntuple = g4sfile['default_ntuples']['g4sntuple']

# build a pandas DataFrame from the hdf5 datasets we will use
g4sdf = pd.DataFrame(np.array(g4sntuple['event']['pages']), columns=['event'])
g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['step' ]['pages']), columns=['step' ]), lsuffix='_caller', rsuffix='_other')
g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['Edep' ]['pages']), columns=['Edep' ]), lsuffix='_caller', rsuffix='_other')
g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['volID']['pages']), columns=['volID']), lsuffix='_caller', rsuffix='_other')
g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['iRep' ]['pages']), columns=['iRep' ]), lsuffix='_caller', rsuffix='_other')

# apply E cut / detID cut and sum Edeps for each event using loc, groupby, and sum
# write directly into output dataframe
detector_hits = g4sdf.loc[(g4sdf.Edep > 0) & (g4sdf.volID == 1)]
procdf = pd.DataFrame(detector_hits.groupby(['event', 'volID', 'iRep'], as_index=False)['Edep'].sum())
procdf = procdf.rename(columns={'iRep': 'detID', 'Edep': 'energy', 'event': 'ievt'})

# write to output file
# procdf.to_hdf('processed.hdf5', key='procdf', mode='w')

bins = 100
plt.figure()
for det, detdf in procdf.groupby('detID'):
    detdf['energy'].hist(bins=bins, histtype='step', label='detector {}'.format(det))
plt.xlabel('energy [MeV]')
plt.gca().set_ylim(0.1, plt.gca().get_ylim()[1])
plt.gca().grid(False)
# plt.gca().axvline(1, lw=1, ls='--', c='k')
plt.yscale('log')
plt.legend(frameon=False, loc='upper right')
plt.show()

# bins = 100
# plt.figure()
# dat = []
# for det, detdf in procdf.groupby('ievt'):
#     # print(detdf)
#     tot = 0
#     for e in detdf['energy'].tolist():
#       tot += e
#     dat.append(tot)
# plt.hist(dat, bins=bins, histtype='step', label='all detectors')
# plt.xlabel('energy [MeV]')
# plt.gca().set_ylim(0.1, plt.gca().get_ylim()[1])
# plt.gca().grid(False)
# # plt.gca().axvline(1, lw=1, ls='--', c='k')
# plt.yscale('log')
# plt.legend(frameon=False, loc='upper right')
# plt.show()
