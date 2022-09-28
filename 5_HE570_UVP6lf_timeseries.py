import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import glob

#path to script location
script_path = os.path.abspath('__file__') # i.e. /path/to/dir/script.py
script_dir = os.path.split(script_path)[0] #i.e. /path/to/dir/


######load UVP6lf particle data (downloaded as .tsv from EcoPart 2022/03/17
rel_path_uvp6 = "UVP6lf/L7_1" # relative path to data file
all_files = glob.glob(os.path.join(rel_path_uvp6, "*.tsv"))

uvp6 = pd.concat((pd.read_csv(f, "\t").assign(depth=os.path.basename(f).split('_')[5]) for f in all_files), ignore_index=True)


aa=0.002342
exp=1.136
imgvolume=0.6
uvp6['area_mm']=aa*(uvp6['area']**exp)
uvp6['ESD_mm']=(4*(uvp6['area_mm']/3.14))**0.5
##exclude all counts larger than 1mm and smaller than 100um
uvp6 = uvp6[uvp6['ESD_mm'] < 1]
uvp6 = uvp6[uvp6['ESD_mm'] >0.1]
uvp6['sampled_volume']=uvp6['imgcount']*imgvolume
uvp6['abundance_L']=uvp6['nbr']/uvp6['sampled_volume']
uvp6['normalized_abundance']=uvp6['abundance_L']/uvp6['area_mm'] ### number /mm2/L
#uvp6['log2abundance'] = np.log2(uvp6.abundance)
uvp6.replace(0,np.nan, inplace=True)
uvp6["datetime"] = pd.to_datetime(uvp6['datetime'], format='%Y%m%d%H%M%S')
##Calculate flux 
uvp6['flux_per_particle_mgC_m_d'] =2.8649 * (uvp6['ESD_mm']/10)**2.24    ###Kiko et al. 2017, divide ESD by 10 to get to cm
uvp6['flux_mgC_m2_d']= uvp6['flux_per_particle_mgC_m_d'] * uvp6['abundance_L']*1000  ###*1000 to get to abundance per m3
##Integrate over size classes
uvp6_flux= uvp6.groupby(['depth','datetime'])['flux_mgC_m2_d'].sum().reset_index()
#uvp6_flux= uvp6_flux.groupby(['depth'])['flux_mgC_m2_d'].mean().reset_index()
#print(uvp6_flux)
####
##plot flux over time
#subset depths
flux_30 = uvp6_flux[(uvp6_flux.depth == '0030')]
flux_80 = uvp6_flux[(uvp6_flux.depth == '0080')]
flux_180 = uvp6_flux[(uvp6_flux.depth == '0180')]
flux_350 = uvp6_flux[(uvp6_flux.depth == '0350')]


fig = plt.figure(1, figsize=(15, 7))
ax1 = fig.add_subplot(1,1,1)
ax1.scatter(flux_30.datetime,flux_30.flux_mgC_m2_d, c='lightblue', label= "30m")
ax1.scatter(flux_80.datetime,flux_80.flux_mgC_m2_d, c='cyan', label= "80m")
ax1.scatter(flux_180.datetime,flux_180.flux_mgC_m2_d, c='blue', label= "180m")
ax1.scatter(flux_350.datetime,flux_350.flux_mgC_m2_d, c='black', label= "350m")
ax1.set_xlabel("Time")
ax1.set_ylabel("Flux (mgC_m_d)")
ax1.legend()
#plt.savefig('HE570_M7_1_flux.pdf')
plt.show()



