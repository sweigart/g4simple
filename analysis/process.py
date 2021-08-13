import h5py
import pandas as pd
import numpy as np

# have to open the input file with h5py (g4 doesn't write pandas-ready hdf5)
g4sfile = h5py.File('g4simpleout.hdf5', 'r')
g4sntuple = g4sfile['default_ntuples']['g4sntuple']

# build a pandas DataFrame from the hdf5 datasets we will use
g4sdf = pd.DataFrame(np.array(g4sntuple['event']['pages']), columns=['event'])
g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['step' ]['pages']), columns=['step' ]), lsuffix = '_caller', rsuffix = '_other')
g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['Edep' ]['pages']), columns=['Edep' ]), lsuffix = '_caller', rsuffix = '_other')
g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['volID']['pages']), columns=['volID']), lsuffix = '_caller', rsuffix = '_other')
g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['iRep' ]['pages']), columns=['iRep' ]), lsuffix = '_caller', rsuffix = '_other')

# apply E cut / detID cut and sum Edeps for each event using loc, groupby, and sum
# write directly into output dataframe
detector_hits = g4sdf.loc[(g4sdf.Edep>0)&(g4sdf.volID==1)]
procdf = pd.DataFrame(detector_hits.groupby(['event','volID','iRep'], as_index=False)['Edep'].sum())
procdf = procdf.rename(columns={'iRep':'detID', 'Edep':'energy'})

# write to output file
procdf.to_hdf('processed.hdf5', key='procdf', mode='w')

