import pandas as pd
import matplotlib.pyplot as plt
import os


######load UVP particle data (downloaded as .tsv from Ecopart 2022/03/17
script_path = os.path.abspath('__file__') # i.e. /path/to/dir/script.py
script_dir = os.path.split(script_path)[0] #i.e. /path/to/dir/

trapuvp = pd.read_csv('UVP6_flux.txt', "\t")
trapuvp['depth'] =trapuvp['depth'].astype(float)

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

#subset columns to keep and pivot from wide to long format for calculations
uvp5 = uvp5.reset_index()
uvp5= pd.melt(uvp5, id_vars=['profile', 'depth', 'Latitude', 'Longitude', 'datetime'], value_vars=[115,144.5,182,229.5,289.5,364.5,459,578.5,729,916.5,1155, 1460])
uvp5.columns = ['profile', 'depth', 'Latitude', 'Longitude', 'datetime','mean_ESD_um', 'abundance_L']
##CHECK THESE CALCS; VALUES SEEM EXTREMELY LOW!!
uvp5['flux_per_particle_mgC_m_d'] =2.8649 * (uvp5['mean_ESD_um']/10000)**2.24    ###Kiko et al. 2017, divide ESD by 10000 to get to cm
uvp5['flux_mgC_m2_d']= uvp5['flux_per_particle_mgC_m_d'] * uvp5['abundance_L']*1000  ###*1000 to get to abundance per m3
##Integrate over size classes
uvp5_flux= uvp5.groupby(['profile', 'depth', 'Latitude', 'Longitude', 'datetime'])['flux_mgC_m2_d'].sum().reset_index()

print(uvp5_flux)

##SELECT MASFORD CENTRAL STATION ONLY
mflux= uvp5_flux.query('60.87 <= Latitude <= 60.875 & 5.405 <= Longitude <= 5.42')
#CALCULATE MEAN AND STDEV
mflux = mflux.groupby('depth').agg({'flux_mgC_m2_d': ['mean', 'std']}).reset_index()
mflux.columns = ['depth','flux_mgC_m2_d','SD']
mflux.reindex(columns=sorted(mflux.columns))

##SELECT LUREJORD CENTRAL STATION ONLY
lflux=uvp5_flux.query('60.69 <= Latitude <= 60.698 & 5.14 <= Longitude <= 5.16')
#CALCULATE MEAN AND STDEV
lflux = lflux.groupby('depth').agg({'flux_mgC_m2_d': ['mean', 'std']}).reset_index()
lflux.columns = ['depth','flux_mgC_m2_d','SD']
lflux.reindex(columns=sorted(lflux.columns))
print(lflux)


######load UVP6 particle data (downloaded as .tsv from EcoPart 2022/03/17
script_path = os.path.abspath('__file__') # i.e. /path/to/dir/script.py
script_dir = os.path.split(script_path)[0] #i.e. /path/to/dir/
rel_path_meta = "UVP6hf_detailed/Export_metadata_summary.tsv" # relative path to data file
rel_path_particles = "UVP6hf_detailed/PAR_Aggregated.tsv" # relative path to data file

abs_file_path_meta = os.path.join(script_dir, rel_path_meta)
abs_file_path = os.path.join(script_dir, rel_path_particles)

uvp6meta = pd.read_csv(abs_file_path_meta, "\t")
uvp6 = pd.read_csv(abs_file_path, "\t")
#merge lat/lon info from metadata file to the dataframe
uvp6 = pd.merge(uvp6, uvp6meta[['profile','Latitude' ,'Longitude']], on=['profile'])
#rename column names for their mean ESD in um
uvp6.rename(columns={
'Depth [m]': 'depth', 
'yyyy-mm-dd hh:mm':'datetime',
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

#subset columns to keep and pivot from wide to long format for calculations
uvp6 = uvp6.reset_index()
uvp6= pd.melt(uvp6, id_vars=['profile', 'depth', 'Latitude', 'Longitude', 'datetime'], value_vars=[115,144.5,182,229.5,289.5,364.5,459,578.5,729,916.5,1155, 1460])
uvp6.columns = ['profile', 'depth', 'Latitude', 'Longitude', 'datetime','mean_ESD_um', 'abundance_L']
##CHECK THESE CALCS; VALUES SEEM EXTREMELY LOW!!
uvp6['flux_per_particle_mgC_m_d'] =2.8649 * (uvp6['mean_ESD_um']/10000)**2.24    ###Kiko et al. 2017, divide ESD by 10000 to get to cm
uvp6['flux_mgC_m2_d']= uvp6['flux_per_particle_mgC_m_d'] * uvp6['abundance_L']*1000  ###*1000 to get to abundance per m3
##Integrate over size classes
uvp6_flux= uvp6.groupby(['profile', 'depth', 'Latitude', 'Longitude', 'datetime'])['flux_mgC_m2_d'].sum().reset_index()

print(uvp6_flux)

##SELECT MASFORD CENTRAL STATION ONLY
mflux6= uvp6_flux.query('60.87 <= Latitude <= 60.875 & 5.405 <= Longitude <= 5.42')
#CALCULATE MEAN AND STDEV
mflux6 = mflux6.groupby('depth').agg({'flux_mgC_m2_d': ['mean', 'std']}).reset_index()
mflux6.columns = ['depth','flux_mgC_m2_d','SD']
mflux6.reindex(columns=sorted(mflux6.columns))

##SELECT LUREJORD CENTRAL STATION ONLY
lflux6=uvp6_flux.query('60.69 <= Latitude <= 60.698 & 5.14 <= Longitude <= 5.16')
#CALCULATE MEAN AND STDEV
lflux6 = lflux6.groupby('depth').agg({'flux_mgC_m2_d': ['mean', 'std']}).reset_index()
lflux6.columns = ['depth','flux_mgC_m2_d','SD']
lflux6.reindex(columns=sorted(lflux6.columns))
print(lflux6)

###UVP6LF traps
trapL7 = trapuvp[(trapuvp.fjord == 'Lurefjord')]
trapM7 = trapuvp[(trapuvp.fjord == 'Masfjord')]

###PLOT
fig = plt.figure(1, figsize=(10, 6))

#######MASFJORDEN
ax1 = fig.add_subplot(1,2,1)
plt.title("A", x=0.05, y=0.92, color="black", fontweight='bold', fontsize = 14)
ax1.set_ylim(480, 0)
ax1.set_xlim(0, 250)
ax1.fill_betweenx(mflux.depth, mflux.flux_mgC_m2_d-mflux.SD, mflux.flux_mgC_m2_d+mflux.SD, facecolor='grey', alpha=0.3)
ax1.plot(mflux.flux_mgC_m2_d, mflux.depth, c='black')
ax1.fill_betweenx(mflux6.depth, mflux6.flux_mgC_m2_d-mflux6.SD, mflux6.flux_mgC_m2_d+mflux6.SD, facecolor='lightblue', alpha=0.3)
ax1.plot(mflux6.flux_mgC_m2_d, mflux6.depth, c='blue')
ax1.scatter(trapM7.flux_mgC_m2_d, trapM7.depth, label = 'UVP6-LF moored')
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.set_xlabel('Flux (mg C$\ m^{-2}\ d^{-1}$)',fontsize=12)
ax1.set_ylabel('Depth (m)',fontsize=12)



#######LUREFJORDEN
ax2 = fig.add_subplot(1,2,2)
plt.title("B", x=0.05, y=0.92, color="black", fontweight='bold', fontsize = 14)
ax2.set_ylim(480, 0)
ax2.set_xlim(0, 250)
ax2.fill_betweenx(lflux.depth, lflux.flux_mgC_m2_d-lflux.SD, lflux.flux_mgC_m2_d+lflux.SD, facecolor='grey', alpha=0.3)
ax2.plot(lflux.flux_mgC_m2_d, lflux.depth, c='black', label = 'UVP5')
ax2.fill_betweenx(lflux6.depth, lflux6.flux_mgC_m2_d-lflux.SD, lflux6.flux_mgC_m2_d+lflux6.SD, facecolor='lightblue', alpha=0.3)
ax2.plot(lflux6.flux_mgC_m2_d, lflux6.depth, c='blue', label = 'UVP6-HF')
ax2.scatter(trapL7.flux_mgC_m2_d, trapL7.depth, label = 'UVP6-LF moored')
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.set_xlabel('Flux (mg C$\ m^{-2}\ d^{-1}$)',fontsize=12)
ax2.set_ylabel('Depth (m)',fontsize=12)
ax2.legend()
###show and save

plt.tight_layout()
#os.chdir('V://Daten/Students/KeaWitting')
plt.savefig('HE570_centralstations_meanflux_UVP5UVP6.pdf')
plt.show()




