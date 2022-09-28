import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import glob

#path to script location
script_path = os.path.abspath('__file__') # i.e. /path/to/dir/script.py
script_dir = os.path.split(script_path)[0] #i.e. /path/to/dir/


######load UVP6lf particle data (downloaded as .tsv from EcoPart 2022/03/17
rel_path_uvp6 = "UVP6lf/M7_1" # relative path to data file
all_files = glob.glob(os.path.join(rel_path_uvp6, "*.tsv"))

uvp6 = pd.concat((pd.read_csv(f, "\t").assign(depth=os.path.basename(f).split('_')[5]) for f in all_files), ignore_index=True)

##calculate size (area and ESD) from pixels
aa=0.002342
exp=1.136
imgvolume=0.6
uvp6['area_mm']=aa*(uvp6['area']**exp)
uvp6['ESD_mm']=(4*(uvp6['area_mm']/3.14))**0.5
uvp6['sampled_volume']=uvp6['imgcount']*imgvolume
uvp6['abundance_L']=uvp6['nbr']/uvp6['sampled_volume']

uvp6.replace(0,np.nan, inplace=True)
uvp6["datetime"] = pd.to_datetime(uvp6['datetime'], format='%Y%m%d%H%M%S')
##add size bins bins and label them with mid ESD like on the other UVPs
bins = [0.0508, 0.064, 0.0806, 0.102, 0.128, 0.161, 0.203, 0.256, 0.323, 0.406
        , 0.512, 0.645, 0.813, 1.02, 1.29, 1.63]
labels = [0.0574, 0.0723, 0.0913, 0.115, 0.1445, 0.182, 0.2295, 0.2895, 0.3645, 0.459, 0.5785
          , 0.729, 0.9165, 1.115, 1.46]
uvp6['ESD_bin_mm'] = pd.cut(uvp6['ESD_mm'], bins=bins, labels=labels)
uvp6= uvp6.groupby(['depth','datetime','ESD_bin_mm']).sum().reset_index()
uvp6= uvp6.groupby(['depth','ESD_bin_mm']).mean().reset_index()
uvp6=uvp6[['depth','ESD_bin_mm', 'abundance_L']]
uvp6['ESD_bin_mm'] = uvp6['ESD_bin_mm'].astype(float)
uvp6['mean_Area'] = (uvp6['ESD_bin_mm'])**2*3.14/4 
uvp6['normalized_abundance']=uvp6['abundance_L']/uvp6['mean_Area'] ### number /mm2/L
uvp6.replace(0,np.nan, inplace=True)
####
#subset depths
m_uvp6_30 = uvp6[(uvp6.depth == '0030')]
m_uvp6_80 = uvp6[(uvp6.depth == '0080')]
m_uvp6_180 = uvp6[(uvp6.depth == '0180')]
m_uvp6_350 = uvp6[(uvp6.depth == '0350')]

##plot normalized abundance vs. area, log-log scale
fig = plt.figure(1, figsize=(7, 7))
ax1 = fig.add_subplot(1,2,1)
ax1.set_title('Masfjord')
ax1.set_ylim([0.001, 1000000])
ax1.set_xlim([0.002, 3])
ax1.plot(m_uvp6_30.mean_Area,m_uvp6_30.normalized_abundance, linestyle='dashed', c = 'lightblue', label = 'UVP6 30-35m')
ax1.plot(m_uvp6_80.mean_Area,m_uvp6_80.normalized_abundance, c = 'cyan', linestyle='dashed',label = 'UVP6 80-85m')
ax1.plot(m_uvp6_180.mean_Area,m_uvp6_180.normalized_abundance, c = 'blue', linestyle='dashed',label = 'UVP6 180-185m')
ax1.plot(m_uvp6_350.mean_Area,m_uvp6_350.normalized_abundance, c = 'black', linestyle='dashed',label = 'UVP6 395-400m')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.set_xlabel("Area ($mm^{2}$)")
ax1.set_ylabel("Normalized Abundance (# $mm^{-2} L^{-1}$)")
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.legend(frameon=False)

plt.savefig('plots/HE570_M7_normalizedabundance_UVP6hf.pdf')
plt.show()


