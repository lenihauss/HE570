import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

#path to script location
script_path = os.path.abspath('__file__') # i.e. /path/to/dir/script.py
script_dir = os.path.split(script_path)[0] #i.e. /path/to/dir/

######load rosina data (from Christian)
rel_path_rosina = "ROSINA/ROSINA_P14-M7-61/Profile14_BinnedData_210521-1815.txt" # relative path to data file
abs_file_path_rosina = os.path.join(script_dir, rel_path_rosina)
rosina = pd.read_csv(abs_file_path_rosina, "\s+", skiprows=26, header=None, usecols=range(21))
rel_path_rosina_depth = "ROSINA/ROSINA_P14-M7-61/BinnedCTD/Depth.txt" # relative path to data file
abs_file_path_rosina_depth = os.path.join(script_dir, rel_path_rosina_depth)
rosina_depth = pd.read_csv(abs_file_path_rosina_depth, "\s+",header=None)
rosina =pd.concat([rosina_depth.reset_index(drop=True),rosina.reset_index(drop=True)], axis=1)

rosina.columns = ['depth', 'image', 0.04875,  0.061425,  0.0773955,  0.0975183,  0.122873,  0.15482,  0.195073,  0.245792,  0.309698,  0.39022,  0.491677,  0.619513,  0.780587,  0.983539,  1.23926,  1.56147,  1.96745,  2.47898,  3.12352,  3.93564]
##add 5m depth bins and label them mid_depth, calculate mean in depth bins
bins = np.arange(0, 505, 5).tolist()
labels = np.arange(2.5, 500, 5).tolist()
rosina['mid_depth'] = pd.cut(rosina['depth'], bins=bins, labels=labels)
rosina= rosina.groupby("mid_depth").mean().reset_index()

##SELECT ONE DEPTH BIN
rosina = rosina[(rosina.mid_depth == 102.5)]
rosina = rosina[[0.04875,  0.061425,  0.0773955,  0.0975183,  0.122873,  0.15482,  0.195073,  0.245792,  0.309698,  0.39022,  0.491677,  0.619513,  0.780587,  0.983539,  1.23926,  1.56147,  1.96745,  2.47898,  3.12352,  3.93564]]

rosina= rosina.T.rename_axis('mean_ESD',axis=0).reset_index()
rosina.columns = ['mean_ESD', 'total_biovolume']
rosina['mean_ESD'] = rosina['mean_ESD'].astype(float)
rosina['mean_Area'] = (rosina['mean_ESD'])**2*3.14/4 ###in mm2
rosina['mean_biovolume'] = ((rosina['mean_ESD'])**3)*3.14/6   ###in mm3
rosina['normalized_biovolume']=rosina['total_biovolume']/(rosina['mean_ESD']/1000)
rosina['abundance']=rosina['total_biovolume']/rosina['mean_biovolume']
rosina['normalized_abundance']=rosina['abundance']/rosina['mean_Area'] 


####
##plot normalized abundance vs. area, log-log scale

fig = plt.figure(1, figsize=(7, 7))
ax1 = fig.add_subplot(1,1,1)
ax1.plot(rosina.mean_Area,rosina.normalized_abundance, c = 'red', linestyle='dashed', label = 'rosina')
ax1.set_xlabel("Area (mm2)")
ax1.set_ylabel("Normalized Abundance (#/mm2/L)")
ax1.set_xscale('log')
ax1.set_yscale('log')
plt.show()
##plot normalized abundance vs. area, log-log scale for rosina and the two UVPs
fig = plt.figure(1, figsize=(7, 7))
ax1 = fig.add_subplot(1,1,1)
#ax1.plot(uvp5_95_100.mean_Area,uvp5_95_100.normalized_abundance, c = 'black', label = 'UVP5')
ax1.plot(uvp6_95_100.mean_ESD_um/1000,uvp6_95_100.normalized_biovolume, c = 'blue',  label = 'UVP6-HF')
ax1.plot(rosina.mean_ESD/1000,rosina.normalized_biovolume, c = 'red', linestyle='dashed', label = 'rosina')
ax1.set_xlabel("ESD (mm)")
ax1.set_ylabel("Normalized Biovolume (mm3/L/mm)")
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.legend()

plt.show()
####

def lingregress(df):
    x_name = df.columns[3] ##mean_Area
    y_name = df.columns[5]
    x = df[x_name]
    y = df[y_name]
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    r_sq = r_value ** 2

    plt.scatter(x,y)
    plt.xlabel(x_name)
    plt.ylabel(y_name)
    plt.title('slope = %.2f , R squared %.2f' % (slope, r_sq))
    plt.plot(x,intercept + slope * x,c='r')

    print('slope = %.2f , R squared %.2f' % (slope, r_sq))
    plt.show()

lingregress(df)

