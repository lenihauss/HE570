# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 16:16:06 2022

@author: hhauss
"""
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import glob

#path to script location
script_path = os.path.abspath('__file__') # i.e. /path/to/dir/script.py
script_dir = os.path.split(script_path)[0] #i.e. /path/to/dir/

######load CPICS data (from Saskia, M7 longtime deployment only, BOP15 and BOP16. No pressure data. Use UVP pressure file to match.
#super annoying format (comma-seperated, but also spaces everywhere)
rel_path_cpics = "CPICS" # relative path to data file
all_files = glob.glob(os.path.join(rel_path_cpics, "*.txt"))
cpics = pd.concat((pd.read_csv(f, ",", ).replace('\s', '', regex=True) for f in all_files), ignore_index=True)
#format datetime
cpics["Date"] = pd.to_datetime(cpics['Date'], format='%Y/%m/%d')
cpics["Time"]= cpics[" Time"].str.split(".", expand = True)[0]
cpics["Time"] = pd.to_datetime(cpics['Time'], format='%H:%M:%S').dt.time
cpics["datetime"] = cpics["Date"] + pd.to_timedelta(cpics['Time'].astype(str))

#load depth data from UVP pressure sensor
rel_path_pressure = "CPICS/UVP5_pressure_files" # relative path to data file
all_files = glob.glob(os.path.join(rel_path_pressure, "*.txt"))
pressure = pd.concat((pd.read_csv(f, "\t") for f in all_files), ignore_index=True)
pressure["datetime"] = pd.to_datetime(pressure['datetime'])
#merge pressure into cpics by secondly pressure data
cpics = pd.merge(cpics, pressure[['datetime','pressure']], on=['datetime'])
cpics = cpics[['datetime','pressure',' ESD_in_um']]

plt.title('CPICS Histogram')
plt.hist(cpics[" ESD_in_um"], bins=10)
plt.xlabel('ESD (in um???)')

cpics.columns = ['depth_min' ,100.0,126.0,158.0,200.0,251.0,316.0,398.0,501.0,631.0,794.0,1000.0,1259.0,1585.0,1995.0,2512.0,3162.0,3981.0,5012.0,6310.0,7943.0,10000.0,12589.0,15849.0,19953.0,25119.0,31623.0,39811.0,50119.0,63096.0,79433.0,100000.0,125893.0,158489.0,199526.0,251189.0,316228.0,398107.0,501187.0,630957.0,794328.0,1000000.0,1258925.0,1584893.0,1995262.0,2511886.0,3162278.0,3981072.0,5011872.0,6309573.0,7943282.0,10000000.0
]
##add 5m depth bins and label them mid_depth, calculate mean in depth bins
cpics['mid_depth'] = cpics['depth_min']+2.5
##SELECT ONE DEPTH BIN
cpics = cpics[(cpics.mid_depth == 102.5)]
cpics = cpics[[100.0,126.0,158.0,200.0,251.0,316.0,398.0,501.0,631.0,794.0,1000.0,1259.0]]

cpics= cpics.T.rename_axis('mean_ESD',axis=0).reset_index()
cpics.columns = ['mean_ESD', 'abundance']
cpics['mean_ESD'] = cpics['mean_ESD'].astype(float)/1000
cpics['mean_Area'] = (cpics['mean_ESD'])**2*3.14/4 ###in mm2
cpics['mean_biovolume'] = ((cpics['mean_ESD'])**3)*3.14/6   ###in mm3
#cpics['normalized_biovolume']=cpics['total_biovolume']/(cpics['mean_ESD']/1000)
#pisco['abundance']=pisco['total_biovolume']/pisco['mean_biovolume']
pisco['normalized_abundance']=pisco['abundance']/pisco['mean_Area'] 
pisco.replace(0,np.nan, inplace=True)

