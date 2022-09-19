import pandas as pd
import matplotlib.pyplot as plt
import os


######load UVP particle data (downloaded as .tsv from EcoPart 2022/05/09
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
part.rename(columns={'Depth [m]': 'depth', 'LPM (64-128 �m) [# l-1]': 'LPM_64_128', 'LPM (128-256 �m) [# l-1]': 'LPM_128_256', 'LPM (256-512 �m) [# l-1]': 'LPM_256_512'}, inplace=True)
print(part.columns)
##SELECT LUREJORD CENTRAL STATION ONLY
lpart= part.query('60.69 <= Latitude <= 60.698 & 5.14 <= Longitude <= 5.16')
#CALCULATE MEAN AND STDEV
lpart = lpart.groupby(['depth'], as_index=False).agg(
                      {'LPM_64_128':['mean','std'],'LPM_128_256':['mean','std'], 'LPM_256_512':['mean','std']})
lpart.columns = ['depth','LPM_64_128','sd_LPM_64_128','LPM_128_256','sd_LPM_128_256', 'LPM_256_512', 'sd_LPM_256_512']
lpart.reindex(columns=sorted(lpart.columns))
print(lpart)

##SELECT MASFORD CENTRAL STATION ONLY
mpart= part.query('60.87 <= Latitude <= 60.875 & 5.41 <= Longitude <= 5.42')
#CALCULATE MEAN AND STDEV
mpart = mpart.groupby(['depth'], as_index=False).agg(
                      {'LPM_64_128':['mean','std'],'LPM_128_256':['mean','std'], 'LPM_256_512':['mean','std']})
mpart.columns = ['depth','LPM_64_128','sd_LPM_64_128','LPM_128_256','sd_LPM_128_256', 'LPM_256_512', 'sd_LPM_256_512']
mpart.reindex(columns=sorted(mpart.columns))
print(mpart)

###PLOT
fig = plt.figure(1, figsize=(10, 6))

#######MASFJORDEN
ax1 = fig.add_subplot(1,2,1)
plt.title("A", x=0.05, y=0.92, color="black", fontweight='bold', fontsize = 14)
ax1.set_ylim(480, 0)
ax1.set_xlim(0, 270)
ax1.fill_betweenx(mpart.depth, mpart.LPM_64_128-mpart.sd_LPM_64_128, mpart.LPM_64_128+mpart.sd_LPM_64_128, facecolor='grey', alpha=0.3)
ax1.plot(mpart.LPM_64_128, mpart.depth, c='black', label = '64-128 µm')
ax1.fill_betweenx(mpart.depth, mpart.LPM_128_256-mpart.sd_LPM_128_256, mpart.LPM_128_256+mpart.sd_LPM_128_256, facecolor='red', alpha=0.3)
ax1.plot(mpart.LPM_128_256, mpart.depth, c='red', label = '128-256 µm')
ax1.fill_betweenx(mpart.depth, mpart.LPM_256_512-mpart.sd_LPM_256_512, mpart.LPM_256_512+mpart.sd_LPM_256_512, facecolor='blue', alpha=0.3)
ax1.plot(mpart.LPM_256_512, mpart.depth, c='blue', label = '256-512 µm')
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.set_xlabel('Particles ($\ L^{-1}$)',fontsize=12)
ax1.set_ylabel('Depth (m)',fontsize=12)



#######LUREFJORDEN
ax2 = fig.add_subplot(1,2,2)
plt.title("B", x=0.05, y=0.92, color="black", fontweight='bold', fontsize = 14)
ax2.set_ylim(480, 0)
ax2.set_xlim(0, 270)
ax2.fill_betweenx(lpart.depth, lpart.LPM_64_128-lpart.sd_LPM_64_128, lpart.LPM_64_128+lpart.sd_LPM_64_128, facecolor='grey', alpha=0.3)
ax2.plot(lpart.LPM_64_128, lpart.depth, c='black', label = '64-128 µm')
ax2.fill_betweenx(lpart.depth, lpart.LPM_128_256-lpart.sd_LPM_128_256, lpart.LPM_128_256+lpart.sd_LPM_128_256, facecolor='red', alpha=0.3)
ax2.plot(lpart.LPM_128_256, lpart.depth, c='red', label = '128-256 µm')
ax2.fill_betweenx(lpart.depth, lpart.LPM_256_512-lpart.sd_LPM_256_512, lpart.LPM_256_512+lpart.sd_LPM_256_512, facecolor='blue', alpha=0.3)
ax2.plot(lpart.LPM_256_512, lpart.depth, c='blue', label = '256-512 µm')
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.set_xlabel('Particles ($\ L^{-1}$)',fontsize=12)
ax2.set_ylabel('Depth (m)',fontsize=12)
plt.legend(loc="lower right", frameon = False)
###show and save

plt.tight_layout()
plt.savefig('HE570_centralstations_meanParticles.pdf')
plt.show()

