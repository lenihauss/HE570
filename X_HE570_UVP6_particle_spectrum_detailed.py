import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy import stats


######load UVP particle data (downloaded as .tsv from EcoPart 2022/03/17
script_path = os.path.abspath('__file__') # i.e. /path/to/dir/script.py
script_dir = os.path.split(script_path)[0] #i.e. /path/to/dir/
rel_path_meta = "UVP6hf_detailed/Export_metadata_summary.tsv" # relative path to data file
rel_path_particles = "UVP6hf_detailed/PAR_Aggregated.tsv" # relative path to data file

abs_file_path_meta = os.path.join(script_dir, rel_path_meta)
abs_file_path_particles = os.path.join(script_dir, rel_path_particles)

meta = pd.read_csv(abs_file_path_meta, "\t")
part = pd.read_csv(abs_file_path_particles, "\t")
#merge lat/lon info from metadata file to the dataframe
part = pd.merge(part, meta[['profile','Latitude' ,'Longitude']], on=['profile'])
#rename some variables for easier column names
part.rename(columns={
'Depth [m]': 'depth', 
'LPM (102-128 �m) [# l-1]'	    :113,
'LPM (128-161 �m) [# l-1]'	    :140,
'LPM (161-203 �m) [# l-1]'	    :180,
'LPM (203-256 �m) [# l-1]'	    :225,
'LPM (256-323 �m) [# l-1]'	    :280,
'LPM (323-406 �m) [# l-1]'	    :370,
'LPM (406-512 �m) [# l-1]'	    :450,
'LPM (512-645 �m) [# l-1]'	    :590,
'LPM (645-813 �m) [# l-1]'	    :700,
'LPM (0.813-1.02 mm) [# l-1]'	:900,
'LPM (1.02-1.29 mm) [# l-1]'	:1100
}, inplace=True)
##too small for UVP

#LPM (1-1.26 �m) [# l-1]	
#LPM (1.26-1.59 �m) [# l-1]	
#LPM (1.59-2 �m) [# l-1]	
#LPM (2-2.52 �m) [# l-1]	
#LPM (2.52-3.17 �m) [# l-1]	
#LPM (3.17-4 �m) [# l-1]	
#LPM (4-5.04 �m) [# l-1]	
#LPM (5.04-6.35 �m) [# l-1]	
#LPM (6.35-8 �m) [# l-1]	
#LPM (8-10.1 �m) [# l-1]	
#LPM (10.1-12.7 �m) [# l-1]	
#LPM (12.7-16 �m) [# l-1]	
#LPM (16-20.2 �m) [# l-1]	
#LPM (20.2-25.4 �m) [# l-1]	
#LPM (25.4-32 �m) [# l-1]	
#LPM (32-40.3 �m) [# l-1]	
#'LPM (40.3-50.8 �m) [# l-1]'	:45,
#'LPM (50.8-64 �m) [# l-1]'	    :58,
#'LPM (64-80.6 �m) [# l-1]'	    :75,
#'LPM (80.6-102 �m) [# l-1]'	    :90,
##too large for UVP
#LPM (1.29-1.63 mm) [# l-1]	
#LPM (1.63-2.05 mm) [# l-1]	
#LPM (2.05-2.58 mm) [# l-1]	
#LPM (2.58-3.25 mm) [# l-1]	
#LPM (3.25-4.1 mm) [# l-1]	
#LPM (4.1-5.16 mm) [# l-1]	
#LPM (5.16-6.5 mm) [# l-1]	
#LPM (6.5-8.19 mm) [# l-1]

print(part.columns)
##SELECT ONE PROFILE AND DEPTH BIN
lpart = part[(part.profile == 'bop_001') & (part.depth == 7.5)]
lpart = lpart[[113,140,180,225,280,370,450,590,700,900,1100]]
df= lpart.T.rename_axis('size',axis=0).reset_index()
df.columns = ['size', 'abundance']
df['log2abundance'] = np.log2(df.abundance)
print(df)

def lingregress(df):
    x_name = df.columns[0]
    y_name = df.columns[2]
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

