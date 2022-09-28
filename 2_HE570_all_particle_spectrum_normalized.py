import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

#path to script location
script_path = os.path.abspath('__file__') # i.e. /path/to/dir/script.py
script_dir = os.path.split(script_path)[0] #i.e. /path/to/dir/

######load PISCO data (from Veit)
rel_path_pisco = "PISCO/HE570_bop_014.csv" # relative path to data file
abs_file_path_pisco = os.path.join(script_dir, rel_path_pisco)
pisco = pd.read_csv(abs_file_path_pisco, ",", header=None, skiprows=1)

pisco.columns = ['depth_min' ,100.0,126.0,158.0,200.0,251.0,316.0,398.0,501.0,631.0,794.0,1000.0,1259.0,1585.0,1995.0,2512.0,3162.0,3981.0,5012.0,6310.0,7943.0,10000.0,12589.0,15849.0,19953.0,25119.0,31623.0,39811.0,50119.0,63096.0,79433.0,100000.0,125893.0,158489.0,199526.0,251189.0,316228.0,398107.0,501187.0,630957.0,794328.0,1000000.0,1258925.0,1584893.0,1995262.0,2511886.0,3162278.0,3981072.0,5011872.0,6309573.0,7943282.0,10000000.0
]
##add 5m depth bins and label them mid_depth, calculate mean in depth bins
pisco['mid_depth'] = pisco['depth_min']+2.5
##SELECT ONE DEPTH BIN
pisco = pisco[(pisco.mid_depth == 82.5)]
pisco = pisco[[100.0,126.0,158.0,200.0,251.0,316.0,398.0,501.0,631.0,794.0,1000.0,1259.0]]

pisco= pisco.T.rename_axis('mean_ESD',axis=0).reset_index()
pisco.columns = ['mean_ESD', 'abundance']
pisco['mean_ESD'] = pisco['mean_ESD'].astype(float)/1000
pisco['mean_Area'] = (pisco['mean_ESD'])**2*3.14/4 ###in mm2
pisco['mean_biovolume'] = ((pisco['mean_ESD'])**3)*3.14/6   ###in mm3
#pisco['normalized_biovolume']=pisco['total_biovolume']/(pisco['mean_ESD']/1000)
#pisco['abundance']=pisco['total_biovolume']/pisco['mean_biovolume']
pisco['normalized_abundance']=pisco['abundance']/pisco['mean_Area'] 
pisco.replace(0,np.nan, inplace=True)

######load ISC data (from Christian)
rel_path_isc = "ROSINA/ISC-P14-M7-61/BinnedData/ParticlesBinned.txt" # relative path to data file
abs_file_path_isc = os.path.join(script_dir, rel_path_isc)
isc = pd.read_csv(abs_file_path_isc, "\s+", header=None)
rel_path_isc_depth = "ROSINA/ISC-P14-M7-61/BinnedData//Depth.txt" # relative path to data file
abs_file_path_isc_depth = os.path.join(script_dir, rel_path_isc_depth)
isc_depth = pd.read_csv(abs_file_path_isc_depth, header=None)
isc =pd.concat([isc_depth.reset_index(drop=True),isc.reset_index(drop=True)], axis=1)

isc.columns = ['depth', 0.05046265
,0.059708213
,0.070647711
,0.083591499
,0.098906795
,0.117028098
,0.138469513
,0.163839337
,0.192539253
,0.229094142
,0.271400245
,0.3211251
,0.379960343
,0.449575141
,0.53194448
,0.629405197
,0.745064131
,0.881167277
,1.042611181
,1.233634187,
'larger']

##add 5m depth bins and label them mid_depth, calculate mean in depth bins
bins = np.arange(0, 505, 5).tolist()
labels = np.arange(2.5, 500, 5).tolist()
isc['mid_depth'] = pd.cut(isc['depth'], bins=bins, labels=labels)
isc= isc.groupby("mid_depth").mean().reset_index()

##SELECT ONE DEPTH BIN
isc = isc[(isc.mid_depth == 82.5)]
isc = isc[[0.05046265
,0.059708213
,0.070647711
,0.083591499
,0.098906795
,0.117028098
,0.138469513
,0.163839337
,0.192539253
,0.229094142
,0.271400245
,0.3211251
,0.379960343
,0.449575141
,0.53194448
,0.629405197
,0.745064131
,0.881167277
,1.042611181
,1.233634187]]

isc= isc.T.rename_axis('mean_ESD',axis=0).reset_index()
isc.columns = ['mean_ESD', 'abundance']
isc['mean_ESD'] = isc['mean_ESD'].astype(float)
isc['mean_Area'] = (isc['mean_ESD'])**2*3.14/4 ###in mm2
isc['mean_biovolume'] = ((isc['mean_ESD'])**3)*3.14/6   ###in mm3
isc['mean_biovolume'] = ((isc['mean_ESD'])**3)*3.14/6   ###in mm3
isc['total_biovolume'] = isc['mean_biovolume']*isc['abundance']  #in mm3/L
isc['normalized_biovolume']=isc['total_biovolume']/(isc['mean_ESD'])
isc['normalized_abundance']=isc['abundance']/isc['mean_Area'] 
isc.replace(0,np.nan, inplace=True)

######load LISST data (from Götz
rel_path_lisst = "LISST/lisst200.txt" # relative path to data file
abs_file_path_lisst = os.path.join(script_dir, rel_path_lisst)
lisst = pd.read_csv(abs_file_path_lisst, "\t")

##add 5m depth bins and label them mid_depth, calculate mean in depth bins
bins = np.arange(0, 505, 5).tolist()
labels = np.arange(2.5, 500, 5).tolist()
lisst['mid_depth'] = pd.cut(lisst['depth'], bins=bins, labels=labels)
lisst= lisst.groupby("mid_depth").mean().reset_index()

##SELECT ONE DEPTH BIN
lisst = lisst[(lisst.mid_depth == 82.5)]
lisst = lisst[['2.23', '2.63', '3.11',
       '3.67', '4.33', '5.11', '6.03', '7.11', '8.39', '9.90', '11.70',
       '13.80', '16.30', '19.20', '22.70', '26.70', '31.60', '37.20', '43.90',
       '51.90', '61.20', '72.20', '85.20', '101.00', '119.00', '140.00',
       '165.00', '195.00', '230.00', '273.00', '324.00', '396.00', '459.00']]

lisst= lisst.T.rename_axis('mean_ESD',axis=0).reset_index()
lisst.columns = ['mean_ESD', 'total_biovolume']
lisst['mean_ESD'] = lisst['mean_ESD'].astype(float)
lisst['mean_Area'] = (lisst['mean_ESD']/1000)**2*3.14/4 ###in mm2
lisst['mean_biovolume'] = ((lisst['mean_ESD']/1000)**3)*3.14/6   ###in mm3
lisst['normalized_biovolume']=lisst['total_biovolume']/(lisst['mean_ESD']/1000)
lisst['abundance']=lisst['total_biovolume']/lisst['mean_biovolume']
lisst['normalized_abundance']=lisst['abundance']/lisst['mean_Area'] 

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
rosina = rosina[(rosina.mid_depth == 82.5)]
rosina = rosina[[0.04875,  0.061425,  0.0773955,  0.0975183,  0.122873,  0.15482,  0.195073,  0.245792,  0.309698,  0.39022,  0.491677,  0.619513,  0.780587,  0.983539,  1.23926,  1.56147,  1.96745]]

rosina= rosina.T.rename_axis('mean_ESD',axis=0).reset_index()
rosina.columns = ['mean_ESD', 'abundance']
rosina['mean_ESD'] = rosina['mean_ESD'].astype(float)
rosina['mean_Area'] = (rosina['mean_ESD'])**2*3.14/4 ###in mm2
rosina['mean_biovolume'] = ((rosina['mean_ESD'])**3)*3.14/6   ###in mm3
rosina['total_biovolume'] = rosina['mean_biovolume']*rosina['abundance']  #in mm3/L
rosina['normalized_biovolume']=rosina['total_biovolume']/(rosina['mean_ESD'])
#rosina['abundance']=rosina['total_biovolume']/rosina['mean_biovolume']
rosina['normalized_abundance']=rosina['abundance']/rosina['mean_Area'] 
rosina.replace(0,np.nan, inplace=True)

######load UVP6 particle data (downloaded as .tsv from EcoPart 2022/03/17
rel_path_uvp6_meta = "UVP6hf_detailed/Export_metadata_summary.tsv" # relative path to data file
rel_path_uvp6 = "UVP6hf_detailed/PAR_Aggregated.tsv" # relative path to data file

abs_file_path_uvp6_meta = os.path.join(script_dir, rel_path_uvp6_meta)
abs_file_path_uvp6 = os.path.join(script_dir, rel_path_uvp6)

uvp6_meta = pd.read_csv(abs_file_path_uvp6_meta, "\t")
uvp6 = pd.read_csv(abs_file_path_uvp6, "\t")
#merge lat/lon info from metadata file to the dataframe
uvp6 = pd.merge(uvp6, uvp6_meta[['profile','Latitude' ,'Longitude']], on=['profile'])


uvp6.rename(columns={
'Depth [m]': 'depth', 
'yyyy-mm-dd hh:mm':'datetime',
'LPM (50.8-64 �m) [# l-1]'	    :57.4,
'LPM (80.6-102 �m) [# l-1]'	    :91.3,
'LPM (102-128 �m) [# l-1]'	    :115,
'LPM (128-161 �m) [# l-1]'	    :144.5,
'LPM (161-203 �m) [# l-1]'	    :182,
'LPM (203-256 �m) [# l-1]'	    :229.5,
'LPM (256-323 �m) [# l-1]'	    :289.5,
'LPM (323-406 �m) [# l-1]'	    :364.5,
'LPM (406-512 �m) [# l-1]'	    :459,
'LPM (512-645 �m) [# l-1]'	    :578.5,
'LPM (645-813 �m) [# l-1]'	    :729,
'LPM (0.813-1.02 mm) [# l-1]'	:916.5,
'LPM (1.02-1.29 mm) [# l-1]'	:1155,
'LPM (1.29-1.63 mm) [# l-1]' :1460
}, inplace=True)


print(uvp6.columns)
##SELECT ONE PROFILE AND DEPTH BIN
#subset columns to keep and pivot from wide to long format for calculations
uvp6 = uvp6.reset_index()
uvp6= pd.melt(uvp6, id_vars=['profile', 'depth', 'Latitude', 'Longitude', 'datetime'], value_vars=[57.4, 91.3, 115,144.5,182,229.5,289.5,364.5,459,578.5,729,916.5,1155, 1460])
uvp6.columns = ['profile', 'depth', 'Latitude', 'Longitude', 'datetime','mean_ESD_um', 'abundance_L']
print(uvp6.columns)
uvp6['mean_Area'] = (uvp6['mean_ESD_um']/1000)**2*3.14/4   ###Area in mm^2
  ###in mrosina['mean_biovolume'] = ((rosina['mean_ESD'])**3)*3.14/6   ###in m
uvp6['mean_biovolume'] = ((uvp6['mean_ESD_um']/1000)**3)*3.14/6   ###biovolume in mm^3
uvp6['total_biovolume'] = uvp6['mean_biovolume']*uvp6['abundance_L'] ###mm3 L-1
uvp6['normalized_biovolume']=uvp6['total_biovolume']/(uvp6['mean_ESD_um']/1000) ### mm3/L/mm 

uvp6['normalized_abundance']=uvp6['abundance_L']/uvp6['mean_Area'] ### number /mm2/L
#uvp6['log2abundance'] = np.log2(uvp6.abundance)
uvp6.replace(0,np.nan, inplace=True)

######load UVP5 particle data (downloaded as .tsv from Ecopart 2022/03/17
rel_path_uvp5_meta = "UVP5_detailed/Export_metadata_summary.tsv" # relative path to data file
rel_path_uvp5 = "UVP5_detailed/PAR_Aggregated.tsv" # relative path to data file

abs_file_path_uvp5_meta = os.path.join(script_dir, rel_path_uvp5_meta)
abs_file_path_uvp5 = os.path.join(script_dir, rel_path_uvp5)

uvp5meta = pd.read_csv(abs_file_path_uvp5_meta, "\t")
uvp5 = pd.read_csv(abs_file_path_uvp5, "\t")
#merge lat/lon info from metadata file to the dataframe
uvp5 = pd.merge(uvp5, uvp5meta[['profile','Latitude' ,'Longitude']], on=['profile'])
#rename column names for their mean ESD in um

uvp5.rename(columns={
'Depth [m]': 'depth', 
'yyyy-mm-dd hh:mm':'datetime',
'LPM (50.8-64 �m) [# l-1]'	    :57.4,
'LPM (80.6-102 �m) [# l-1]'	    :91.3,
'LPM (102-128 �m) [# l-1]'	    :115,
'LPM (128-161 �m) [# l-1]'	    :144.5,
'LPM (161-203 �m) [# l-1]'	    :182,
'LPM (203-256 �m) [# l-1]'	    :229.5,
'LPM (256-323 �m) [# l-1]'	    :289.5,
'LPM (323-406 �m) [# l-1]'	    :364.5,
'LPM (406-512 �m) [# l-1]'	    :459,
'LPM (512-645 �m) [# l-1]'	    :578.5,
'LPM (645-813 �m) [# l-1]'	    :729,
'LPM (0.813-1.02 mm) [# l-1]'	:916.5,
'LPM (1.02-1.29 mm) [# l-1]'	:1155,
'LPM (1.29-1.63 mm) [# l-1]' :1460
}, inplace=True)

##SELECT ONE PROFILE AND DEPTH BIN
#subset columns to keep and pivot from wide to long format for calculations
uvp5 = uvp5.reset_index()
uvp5= pd.melt(uvp5, id_vars=['profile', 'depth', 'Latitude', 'Longitude', 'datetime'], value_vars=[57.4, 91.3, 115,144.5,182,229.5,289.5,364.5,459,578.5,729,916.5,1155, 1460])
uvp5.columns = ['profile', 'depth', 'Latitude', 'Longitude', 'datetime','mean_ESD_um', 'abundance_L']
uvp5['mean_Area'] = (uvp5['mean_ESD_um']/1000)**2*3.14/4   ###Area in mm^2
uvp5['mean_biovolume'] = ((uvp5['mean_ESD_um']/1000)**3)*3.14/6   ###biovolume in mm^3
uvp5['total_biovolume'] = uvp5['mean_biovolume']*uvp5['abundance_L'] ###mm3 L-1
uvp5['normalized_biovolume']=uvp5['total_biovolume']/(uvp5['mean_ESD_um']/1000) ### mm3/L/mm 

uvp5['normalized_abundance']=uvp5['abundance_L']/uvp5['mean_Area'] 
#uvp5['log2abundance'] = np.log2(uvp5.abundance)
uvp5.replace(0,np.nan, inplace=True)


##UVP Masfjord only
m_uvp6 = uvp6[(uvp6.profile == 'bop_he570_014')]
m_uvp6_30_35 = m_uvp6[(m_uvp6.depth == 32.5)]
m_uvp6_80_85 = m_uvp6[(m_uvp6.depth == 82.5)]
m_uvp6_180_185 = m_uvp6[(m_uvp6.depth == 182.5)]
m_uvp6_350_455 = m_uvp6[(m_uvp6.depth == 352.5)]

m_uvp5 = uvp5[(uvp5.profile == 'bop_014')]
m_uvp5_30_35 = m_uvp5[(m_uvp5.depth == 32.5)]
m_uvp5_80_85 = m_uvp5[(m_uvp5.depth == 82.5)]
m_uvp5_180_185 = m_uvp5[(m_uvp5.depth == 182.5)]
m_uvp5_350_455 = m_uvp5[(m_uvp5.depth == 352.5)]

##UVP Lurefjord only
l_uvp6 = uvp6[(uvp6.profile == 'bop_he570_032')]
l_uvp6_30_35 = l_uvp6[(l_uvp6.depth == 32.5)]
l_uvp6_80_85 = l_uvp6[(l_uvp6.depth == 82.5)]
l_uvp6_180_185 = l_uvp6[(l_uvp6.depth == 182.5)]
l_uvp6_350_455 = l_uvp6[(l_uvp6.depth == 352.5)]

l_uvp5 = uvp5[(uvp5.profile == 'bop_032')]
l_uvp5_30_35 = l_uvp5[(l_uvp5.depth == 32.5)]
l_uvp5_80_85 = l_uvp5[(l_uvp5.depth == 82.5)]
l_uvp5_180_185 = l_uvp5[(l_uvp5.depth == 182.5)]
l_uvp5_350_455 = l_uvp5[(l_uvp5.depth == 352.5)]
####
####
##plot normalized abundance vs. area, log-log scale
fig = plt.figure(1, figsize=(7, 7))
ax1 = fig.add_subplot(1,2,1)
ax1.set_title('Masfjord')
ax1.set_ylim([0.001, 1000000])
ax1.set_xlim([0.002, 3])
ax1.plot(m_uvp5_30_35.mean_Area,m_uvp5_30_35.normalized_abundance, c = 'lightblue', label = 'UVP5 30-35m')
ax1.plot(m_uvp5_80_85.mean_Area,m_uvp5_80_85.normalized_abundance, c = 'cyan', label = 'UVP5 80-85m')
ax1.plot(m_uvp5_180_185.mean_Area,m_uvp5_180_185.normalized_abundance, c = 'blue', label = 'UVP5 180-185m')
ax1.plot(m_uvp5_350_455.mean_Area,m_uvp5_350_455.normalized_abundance, c = 'black', label = 'UVP5 395-400m')
ax1.plot(m_uvp6_30_35.mean_Area,m_uvp6_30_35.normalized_abundance, linestyle='dashed', c = 'lightblue', label = 'UVP6 30-35m')
ax1.plot(m_uvp6_80_85.mean_Area,m_uvp6_80_85.normalized_abundance, c = 'cyan', linestyle='dashed',label = 'UVP6 80-85m')
ax1.plot(m_uvp6_180_185.mean_Area,m_uvp6_180_185.normalized_abundance, c = 'blue', linestyle='dashed',label = 'UVP6 180-185m')
ax1.plot(m_uvp6_350_455.mean_Area,m_uvp6_350_455.normalized_abundance, c = 'black', linestyle='dashed',label = 'UVP6 395-400m')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.set_xlabel("Area ($mm^{2}$)")
ax1.set_ylabel("Normalized Abundance (# $mm^{-2} L^{-1}$)")
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.legend(frameon=False)
ax2 = fig.add_subplot(1,2,2)
ax2.set_title('Lurefjord')
ax2.set_ylim([0.001, 1000000])
ax2.set_xlim([0.002, 3])
ax2.plot(l_uvp5_30_35.mean_Area,l_uvp5_30_35.normalized_abundance, c = 'lightblue', label = '30-35m')
ax2.plot(l_uvp5_80_85.mean_Area,l_uvp5_80_85.normalized_abundance, c = 'cyan', label = '80-85m')
ax2.plot(l_uvp5_180_185.mean_Area,l_uvp5_180_185.normalized_abundance, c = 'blue', label = '180-185m')
ax2.plot(l_uvp5_350_455.mean_Area,l_uvp5_350_455.normalized_abundance, c = 'black', label = '395-400m')
ax2.plot(l_uvp6_30_35.mean_Area,l_uvp6_30_35.normalized_abundance, linestyle='dashed', c = 'lightblue', label = '30-35m')
ax2.plot(l_uvp6_80_85.mean_Area,l_uvp6_80_85.normalized_abundance, c = 'cyan', linestyle='dashed',label = 'UVP6 80-85m')
ax2.plot(l_uvp6_180_185.mean_Area,l_uvp6_180_185.normalized_abundance, c = 'blue', linestyle='dashed',label = 'UVP6 180-185m')
ax2.plot(l_uvp6_350_455.mean_Area,l_uvp6_350_455.normalized_abundance, c = 'black', linestyle='dashed',label = 'UVP6 395-400m')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.set_xlabel("Area ($mm^{2}$)")
ax2.set_xscale('log')
ax2.set_yscale('log')
plt.savefig('plots/HE570_M7_100m_normalizedabundance_UVP5and6.pdf')
plt.show()
##plot normalized abundance vs. Area, log-log scale
fig = plt.figure(1, figsize=(7, 7))
ax1 = fig.add_subplot(1,1,1)
#ax1.set_ylim([0.001, 1000000])
#ax1.set_xlim([0.001, 3])
ax1.plot(m_uvp5_80_85.mean_Area,m_uvp5_80_85.normalized_abundance, c = 'black', label = 'UVP5')
ax1.plot(m_uvp6_80_85.mean_Area,m_uvp6_80_85.normalized_abundance, c = 'blue',  label = 'UVP6-HF')
ax1.plot(rosina.mean_Area,rosina.normalized_abundance, c = 'green',  label = 'ROSINA')
ax1.plot(isc.mean_Area,isc.normalized_abundance, c = 'cyan',  label = 'ISC')
ax1.plot(pisco.mean_Area,pisco.normalized_abundance, c = 'orange',  label = 'PISCO')
ax1.plot(lisst.mean_Area,lisst.normalized_abundance, c = 'red', linestyle='dashed', label = 'LISST')

ax1.set_xlabel("Area ($mm^{2}$)")
ax1.set_ylabel("Normalized Abundance (# $mm^{-2} L^{-1}$)")
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.legend(frameon=False)
plt.savefig('plots/HE570_M7_80m_normalizedabundance_UVPs_ROSINA_ISC_PISCO_LISST.pdf')
plt.show()



##plot normalized biovolume vs. area, log-log scale for LISST and the two UVPs
fig = plt.figure(1, figsize=(7, 7))
ax1 = fig.add_subplot(1,1,1)
ax1.plot(m_uvp5_80_85.mean_ESD_um/1000,m_uvp5_80_85.normalized_biovolume, c = 'black', label = 'UVP5')
ax1.plot(m_uvp6_80_85.mean_ESD_um/1000,m_uvp6_80_85.normalized_biovolume, c = 'blue',  label = 'UVP6-HF')
ax1.plot(lisst.mean_ESD/1000,lisst.normalized_biovolume, c = 'red', linestyle='dashed', label = 'LISST')
ax1.plot(rosina.mean_ESD,rosina.normalized_biovolume, c = 'green', label = 'ROSINA')
ax1.plot(isc.mean_ESD,isc.normalized_biovolume, c = 'cyan', label = 'ISC')
ax1.set_xlabel("ESD (mm)")
ax1.set_ylabel("Normalized Biovolume ($mm^{3} L^{-1}mm^{-1}$)")
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.legend()
plt.savefig('plots/HE570_M7_80m_normalizedbiovolume_UVPs_ROSINA_ISC_PISCO_LISST.pdf')

plt.show()
####



