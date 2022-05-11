import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy import stats


######load UVP particle data (downloaded as .tsv from EcoPart 2022/03/17
script_path = os.path.abspath('__file__') # i.e. /path/to/dir/script.py
script_dir = os.path.split(script_path)[0] #i.e. /path/to/dir/
rel_path_meta = "UVP5_reduced/Export_metadata_summary.tsv" # relative path to data file
rel_path_particles = "UVP5_reduced/PAR_Aggregated.tsv" # relative path to data file

abs_file_path_meta = os.path.join(script_dir, rel_path_meta)
abs_file_path_particles = os.path.join(script_dir, rel_path_particles)

meta = pd.read_csv(abs_file_path_meta, "\t")
part = pd.read_csv(abs_file_path_particles, "\t")


#merge lat/lon info from metadata file to the dataframe
part = pd.merge(part, meta[['profile','Latitude' ,'Longitude']], on=['profile'])
#rename some variables for easier column names
part.rename(columns={
'Depth [m]': 'depth', 
'LPM (64-128 �m) [# l-1]': 96, 
'LPM (128-256 �m) [# l-1]': 192, 
'LPM (256-512 �m) [# l-1]': 384, 
'LPM (0.512-1.02 mm) [# l-1]': 766, 
'LPM (1.02-2.05 mm) [# l-1]': 1535
}, inplace=True)
##mostly empty size bins
#'LPM (2.05-4.1 mm) [# l-1]': 3075,
#'LPM (4.1-8.19 mm) [# l-1]': 6145
#LPM (8.19-16.4 mm) [# l-1]
#LPM (>16.4 mm) [# l-1]

print(part.columns)
##SELECT ONE PROFILE AND DEPTH BIN
lpart = part[(part.profile == 'bop_001') & (part.depth == 7.5)]
lpart = lpart[[96, 192, 384, 766, 1535]]
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


