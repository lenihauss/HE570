import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy import stats


######load UVP particle data (downloaded as .tsv from EcoPart 2022/03/17
script_path = os.path.abspath('__file__') # i.e. /path/to/dir/script.py
script_dir = os.path.split(script_path)[0] #i.e. /path/to/dir/
rel_path_particles = "LISST/lisst200.txt" # relative path to data file

abs_file_path_particles = os.path.join(script_dir, rel_path_particles)

df = pd.read_csv(abs_file_path_particles, "\t")

##add 5m depth bins and label them mid_depth
bins = np.arange(0, 505, 5).tolist()
labels = np.arange(2.5, 500, 5).tolist()
df['mid_depth'] = pd.cut(df['depth'], bins=bins, labels=labels)
part= df.groupby("mid_depth").mean().reset_index()
print(part)

print(part.columns)
##SELECT ONE DEPTH BIN
lpart = part[(part.mid_depth == 102.5)]
lpart = lpart[['1.21', '1.60', '1.89', '2.23', '2.63', '3.11',
       '3.67', '4.33', '5.11', '6.03', '7.11', '8.39', '9.90', '11.70',
       '13.80', '16.30', '19.20', '22.70', '26.70', '31.60', '37.20', '43.90',
       '51.90', '61.20', '72.20', '85.20', '101.00', '119.00', '140.00',
       '165.00', '195.00', '230.00', '273.00', '324.00', '396.00', '459.00']]
print(lpart)
df= lpart.T.rename_axis('mean_ESD',axis=0).reset_index()
df.columns = ['mean_ESD', 'total_biovolume']
df['mean_ESD'] = df['mean_ESD'].astype(float)
df['mean_Area'] = (df['mean_ESD']/1000)**2*3.14/4 ###in mm2
df['mean_Biovolume'] = ((df['mean_ESD']/1000)**3)*3.14/6   ###in mm3
df['abundance']=df['total_biovolume']/df['mean_Biovolume']
print(df)

df['normalized_abundance']=df['abundance']/df['mean_Area'] 

print(df)

##plot normalized abundance vs. area, log-log scale
fig = plt.figure(1, figsize=(7, 7))
ax1 = fig.add_subplot(1,2,1)
ax1.plot(df.mean_Area,df.normalized_abundance)
#ax1.set_xlabel("Area")
#ax1.set_ylabel("abundance")
ax1.set_xscale('log')
ax1.set_yscale('log')
plt.show()
####

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

