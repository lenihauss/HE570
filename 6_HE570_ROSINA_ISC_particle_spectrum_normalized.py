import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

#path to script location
script_path = os.path.abspath('__file__') # i.e. /path/to/dir/script.py
script_dir = os.path.split(script_path)[0] #i.e. /path/to/dir/


######load ISC data (from Christian)
rel_path_isc = "ROSINA/ISC-P14-M7-61/BinnedData/ParticlesBinned.txt" # relative path to data file
abs_file_path_isc = os.path.join(script_dir, rel_path_isc)
isc = pd.read_csv(abs_file_path_isc, "\s+", header=None)
rel_path_isc_depth = "ROSINA/ISC-P14-M7-61/BinnedData//Depth.txt" # relative path to depth file
abs_file_path_isc_depth = os.path.join(script_dir, rel_path_isc_depth)
isc_depth = pd.read_csv(abs_file_path_isc_depth, header=None)
isc =pd.concat([isc_depth.reset_index(drop=True),isc.reset_index(drop=True)], axis=1)

isc.columns = ['depth', 0.15853,
0.18758,
0.22195,
0.26261,
0.31072,
0.36765,
0.43501,
0.51472,
0.60902,
0.72060,
0.85263,
1.00880,
1.19370,
1.41240,
1.67120,
1.97730,
2.33960,
2.76830,
3.27550,
3.87560,
'larger']
##add 5m depth bins and label them mid_depth, calculate mean in depth bins
bins = np.arange(0, 505, 5).tolist()
labels = np.arange(2.5, 500, 5).tolist()
isc['mid_depth'] = pd.cut(isc['depth'], bins=bins, labels=labels)
isc= isc.groupby("mid_depth").mean().reset_index()

##SELECT ONE DEPTH BIN
isc = isc[(isc.mid_depth == 102.5)]
isc = isc[[0.15853,
0.18758,
0.22195,
0.26261,
0.31072,
0.36765,
0.43501,
0.51472,
0.60902,
0.72060,
0.85263,
1.00880,
1.19370,
1.41240,
1.67120,
1.97730,
2.33960,
2.76830,
3.27550,
3.87560]]

isc= isc.T.rename_axis('mean_ESD',axis=0).reset_index()
isc.columns = ['mean_ESD', 'abundance']
isc['mean_ESD'] = isc['mean_ESD'].astype(float)
isc['mean_Area'] = (isc['mean_ESD'])**2*3.14/4 ###in mm2
isc['mean_biovolume'] = ((isc['mean_ESD'])**3)*3.14/6   ###in mm3
isc['normalized_abundance']=isc['abundance']/isc['mean_Area'] 
isc.replace(0,np.nan, inplace=True)


######load rosina data (from Christian)
rel_path_rosina = "ROSINA/ROSINA_P14-M7-61/Profile14_BinnedData_210521-1815.txt" # relative path to data file
abs_file_path_rosina = os.path.join(script_dir, rel_path_rosina)
rosina = pd.read_csv(abs_file_path_rosina, "\s+", skiprows=26, header=None, usecols=range(18))
rel_path_rosina_depth = "ROSINA/ROSINA_P14-M7-61/BinnedCTD/Depth.txt" # relative path to data file
abs_file_path_rosina_depth = os.path.join(script_dir, rel_path_rosina_depth)
rosina_depth = pd.read_csv(abs_file_path_rosina_depth, "\s+",header=None)
rosina =pd.concat([rosina_depth.reset_index(drop=True),rosina.reset_index(drop=True)], axis=1)

rosina.columns = ['depth', 'image', 0.04875,  0.061425,  0.0773955,  0.0975183,  0.122873,  0.15482,  0.195073,  0.245792,  0.309698,  0.39022,  0.491677,  0.619513,  0.780587,  0.983539,  1.23926,  1.56147,  1.96745]
##add 5m depth bins and label them mid_depth, calculate mean in depth bins
bins = np.arange(0, 505, 5).tolist()
labels = np.arange(2.5, 500, 5).tolist()
rosina['mid_depth'] = pd.cut(rosina['depth'], bins=bins, labels=labels)
rosina= rosina.groupby("mid_depth").mean().reset_index()

##SELECT ONE DEPTH BIN
rosina = rosina[(rosina.mid_depth == 102.5)]
rosina = rosina[[0.04875,  0.061425,  0.0773955,  0.0975183,  0.122873,  0.15482,  0.195073,  0.245792,  0.309698,  0.39022,  0.491677,  0.619513,  0.780587,  0.983539,  1.23926,  1.56147,  1.96745]]

rosina= rosina.T.rename_axis('mean_ESD',axis=0).reset_index()
rosina.columns = ['mean_ESD', 'abundance']
rosina['mean_ESD'] = rosina['mean_ESD'].astype(float)
rosina['mean_Area'] = (rosina['mean_ESD'])**2*3.14/4 ###in mm2
rosina['mean_biovolume'] = ((rosina['mean_ESD'])**3)*3.14/6   ###in mm3
rosina['normalized_abundance']=rosina['abundance']/rosina['mean_Area'] 
rosina.replace(0,np.nan, inplace=True)


##plot normalized abundance vs. area, log-log scale
fig = plt.figure(1, figsize=(7, 7))
ax1 = fig.add_subplot(1,1,1)
ax1.plot(rosina.mean_Area,rosina.normalized_abundance, c = 'green',  label = 'ROSINA')
ax1.plot(isc.mean_Area,isc.normalized_abundance, c = 'cyan',  label = 'ISC')
ax1.set_xlabel("Area ($mm^{2}$)")
ax1.set_ylabel("Normalized Abundance (# $mm^{-2} L^{-1}$)")
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.legend()
plt.show()
