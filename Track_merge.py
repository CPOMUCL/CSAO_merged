import numpy as np
import datetime as dt
from dateutil.relativedelta import relativedelta
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from mpl_toolkits.basemap import Basemap
import glob
from netCDF4 import Dataset
import copy

import os
import sys
# sys.path.insert(0, '/Users/H/WAVES/geo_data_group/')
import grid_set as gs

import ant_plus


t_start = dt.datetime(2019,9,1) 
n_days = 30

### To make the merged files
MERGE_TRACKS = True

### To make the merged plots
PLOT_RANGE = False
PLOT_RANGE_ANOM = False
PLOT_EXTRA_ANOM = False
PLOT_INTERP = False

### Diagnostic to just give info on a single merged track
check_single_track = False
check_filter ='_' # change this string to limit printed info ie: 'range' 
check_filter ='ANOM' # Just ANOMALIES


#### Directories of the merged files and figures
merge_dir = '/home/hds/DATA/CSAO/MERGED/'
# merge_dir = './test_merge/'
# fig_dir = './ALL_FIGS/2019-09_merged_all/'
fig_dir = './ALL_FIGS/2019-09_Amundsen/'
# fig_dir = './ALL_FIGS/2019-09_West_Weddell/'

##### Area restriction (good for combining tracks)
lons = [-180,180] ### Default - whole Antarctic
lats = [-89,-50]

# lons = [-150,-80] ### Amundsen
# lats = [-73,-63]
# lons = [-61,-40] ### West Weddell
# lats = [-76,-63]

#### To combine certain data into a single file
COMBINE_FILES = True
# Combine_name = 'L2_ice_leads_Amundsen_'
# Combine_name = 'L2_ice_leads_West_Weddell_'
Combine_name = 'ST_and_range_ANOM_all_'
def Comb_ST(track):  
    ### sar_ocean  l2
#     return track.flag_surf_type_class_20_ku == 64 
#     #### Comb_ST = 128 #### sar_sea_ice 
#     return track.flag_surf_type_class_20_ku == 128
#     #### Comb_ST = 256 #### sar_lead
#     return track.flag_surf_type_class_20_ku == 256
    #### Comb_ST = 256 #### sar_ice and lead
#     return track.flag_surf_type_class_20_ku >= 128
    #### Comb_ST = 256 #### ocean sar_ice and lead
    return track.flag_surf_type_class_20_ku >= 32
#     #### Comb_ST = 256 #### CLS_lead
#     return track.flag_surf_type_20_ku_CLS == 3
#     #### Comb_ST = 256 #### LEGOS GPOD lead
#     return track.flag_surf_type_20_ku_LEGOS_GPOD == 3
#     #### Multiple options
#     #### L2 lead and CLS lead
#     mask = track.flag_surf_type_20_ku_CLS == 3
#     mask[track.flag_surf_type_class_20_ku != 256] = 0
#     return mask
Combine_vars= [
#  'ocean_tide_20_ku_ANOM_DTU_fes',
 'range_1_20_ku',
 'range_1_20_ku_ANOM_ISat_2step',
 'range_1_20_ku_ANOM_ISat_SWH_MSSfixed_Ice',
 'range_1_20_ku_ANOM_DTU_imp_th70',
 'range_1_20_ku_ANOM_DTU_th70',
 'range_1_20_ku_ANOM_DTU_gauss',
 'range_1_20_ku_ANOM_MSSL_D',
 'range_1_20_ku_ANOM_MSSL_S',
 'range_1_20_ku_ANOM_LEGOS_GPOD_All',
 'range_1_20_ku_ANOM_CLS',
#  'mean_sea_surf_20_ku_ANOM_DTU15',
#  'mean_sea_surf_20_ku_ANOM_LEGOS_GPOD_DTU15',
#  'mean_sea_surf_20_ku_ANOM_LEGOS_GPOD_DTU21',
#  'atm_geo_corrections_sum_ANOM_ISat',
#  'atm_geo_corrections_sum_ANOM_LEGOS_GPOD',
#  'mean_sea_surf_sea_ice_20_ku','geoid_20_ku',
#  'height_1_20_ku',
#  'ssha_interp_20_ku','ssha_20_ku','freeboard_20_ku',
#  'ssha_interp_20_ku_LEGOS_GPOD_Ice',
#  'ssha_interp_20_ku_LEGOS_GPOD_Lead',
#  'ssha_smooth_20_ku_LEGOS_GPOD_Ice',
#  'ssha_smooth_20_ku_LEGOS_GPOD_Lead',
#  'ssha_interp_20_ku_ANOM_LEGOS_GPOD',
 'flag_surf_type_class_20_ku',
 'flag_surf_type_20_ku_CLS',
 'flag_surf_type_20_ku_LEGOS_GPOD',
 'NSIDC_nasa',
                ]




#### etra output verbos flags
v_tload = False
v_tmask = False
v_save  = False


if PLOT_RANGE or PLOT_RANGE_ANOM or PLOT_EXTRA_ANOM or PLOT_INTERP:
    check_dir = os.path.dirname(fig_dir)
    # check it exists
    if not os.path.exists(check_dir):
        # make if it doesn't
        os.makedirs(check_dir)
        print('Creating directory: ',check_dir)
    else: print('Existing directory: ',check_dir)
    
if MERGE_TRACKS:
    check_dir = os.path.dirname(merge_dir)
    # check it exists
    save_root = merge_dir
    if not os.path.exists(check_dir):
        # make if it doesn't
        os.makedirs(check_dir)
        print('Creating directory: ',check_dir)
    else: print('Existing directory: ',check_dir)
        
read_comb_attr = True
if COMBINE_FILES:
    Comb_list = ant_plus.CS_track_list(Combine_name,Combine_vars)

### projection for plotting and perhaps gridded data too
m = Basemap(projection='spstere', lon_0=0.0, lat_0 = -90,boundinglat=-45)
# G = gs.grid_set(m)
# G.load_grid('/Users/H/WAVES/DOT_processing/grids/Polar_stereo_50km_SH.npz')
# G.load_mask('/Users/H/WAVES/DOT_processing/grids/Polar_stereo_50km_SH_mask.npz')
# G.reproject(m)

tw = relativedelta(days=1)
# tw = relativedelta(hours=1)

### NSIDC grid first
GIce = gs.grid_set(m)
# GIce.set_grid_lon_lat(np.flipud(Nlons),np.flipud(Nlats))
GIce.load_grid('NSIDC_gs_SH.npz')
GIce.reproject(m)
GIce.get_ptp()


for nd in range(n_days):
    ### START WITH A SINGLE DAY
    time = t_start + relativedelta(days=nd)


    #### here is some code for track to track comparisons for Ant+
    #### we start with L2I files as a base, all comparisons 
    #### will be 'virtually' added to this files to compare

    ### this is finding and sorting the L2 files
    ### SAR L2
#     l1b_root = '/Volumes/BU_extra/CryoSat/L2I/'
    l1b_root = '/cpnet/li1_cpdata/SATS/RA/CRY/L2I/'
    mode = 'SAR-A'
    # time = dt.datetime(2014,3,1)
    temp_files = glob.glob(l1b_root+mode+time.strftime('/%Y/%m/')+'*.nc')
    temp_files.sort()
    SARfiles = []
    SARtime = []
    SARtime_end = []
    tstring = time.strftime('%Y%m%d')
    for file in temp_files:
        if tstring in file:
            try:
                ftime = dt.datetime.strptime(
                    tstring+file.split(tstring)[1],
                    '%Y%m%dT%H%M%S_')
                SARtime.append(ant_plus.dt2CSt(ftime))
                SARfiles.append(file)
                try:
                    ftime = dt.datetime.strptime(
                    tstring+file.split(tstring)[2].split('_')[0],
                    '%Y%m%dT%H%M%S')
                    SARtime_end.append(ant_plus.dt2CSt(ftime))
                except ValueError:
                    SARtime_end.append(np.nan)
            except ValueError:
                pass
    SARtime = np.array(SARtime)
    SARtime_end = np.array(SARtime_end)
    Lmax = len(SARfiles)
    print(len(SARfiles))


    L2_vars = ['range_1_20_ku','alt_20_ku', 
           'flag_surf_type_class_20_ku','flag_height_20_ku',
           'height_1_20_ku','ssha_20_ku','ssha_interp_20_ku',
           'mean_sea_surf_sea_ice_20_ku','geoid_20_ku',
           'flag_freeboard_20_ku','freeboard_20_ku',
          ]
    L2_vars_1Hz= ['time_cor_01',
            'mod_wet_tropo_cor_01','mod_dry_tropo_cor_01', #### atmo corrections (also called geo cor)
            'pole_tide_01','solid_earth_tide_01','load_tide_01',
           'ocean_tide_01','ocean_tide_eq_01', #### geo corrections
           'iono_cor_gim_01','inv_bar_cor_01','hf_fluct_total_cor_01',
            ]

    ### CLS files - 
#     Cpath = '/Volumes/BU_extra/CryoSat/Antarctica+/CLS/'
    Cpath = '/home/hds/DATA2/Ant+/DATA/CLS/'
    Cfiles = ant_plus.add_file('CLS',Cpath,'%b_%Y/',multi_file=False)
    Cfiles.tvec = 'time_20_ku'
    # Cfiles.vars = ['range_rtk_cls_20_ku','surf_type_cls_20_ku']
    # Cfiles.names = ['CLS_range','CLS_surf_type']
    Cfiles.vars = ['range_rtk_cls_20_ku','surf_type_cls_20_ku','ssh_rtk_cls_20_ku']
    Cfiles.names = ['range_1_20_ku_CLS','flag_surf_type_20_ku_CLS','height_1_20_ku_CLS']


#     ### MSSL compare bas range - date from time gaps
    ### MSSL compare difuse- date from time gaps
    Mpath = '/raid6/userdata/hds/CSAO/MSSL/CSAO_RTRK_MSSL_2019/'
#     Mpath = '/Volumes/BU_extra/CryoSat/Antarctica+/MSSL/CSAO_RTRK_MSSL_2019/'
    Mfiles = ant_plus.add_file('MSSL',Mpath,'',multi_file=False,name_form = '%Y%m')
    Mfiles.tvec = 'time_20_ku'
    Mfiles.vars = ['range_20_ku']
    Mfiles.names = ['range_20_ku_MSSL']
#     Mfiles.vars = ['range_20_ku','flag_surf_type_disc_20_ku']
#     Mfiles.names = ['MSSL_range','MSSL_surf_type']


    ### MSSL compare difuse- date from time gaps
#     MDpath = '/Volumes/BU_extra/CryoSat/Antarctica+/MSSL/CSAO_DRSA_MSSL_2019/'
    MDpath = '/raid6/userdata/hds/CSAO/MSSL/CSAO_DRSA_MSSL_2019/'
    MDfiles = ant_plus.add_file('MSSL_D',MDpath,'',name_form = '%Y%m',multi_file=False)
    MDfiles.tvec = 'time_20_ku'
    MDfiles.vars = ['retracker_cor_20_ku']#,'retracker_cor_basd_20_ku']
    MDfiles.names = ['range_MSSL_Dcor']
#     MDfiles.vars = ['ssha_20_ku']
#     MDfiles.names = ['MSSL_D_ssha']

    ### MSSL compare specular - date from time gaps
#     MSpath = '/Volumes/BU_extra/CryoSat/Antarctica+/MSSL/CSAO_SRSA_MSSL_2019/'
    MSpath = '/raid6/userdata/hds/CSAO/MSSL/CSAO_SRSA_MSSL_2019/'
    MSfiles = ant_plus.add_file('MSSL_S',MSpath,'',name_form = '%Y%m',multi_file=False)
    MSfiles.tvec = 'time_20_ku'
    MSfiles.vars = ['retracker_cor_20_ku']
    MSfiles.names = ['range_MSSL_Scor']
#     MSfiles.vars = ['ssha_20_ku']
#     MSfiles.names = ['MSSL_S_ssha']
    # MSfiles.names = ['MSSL_S_cor']
    
    ### DTU files - date from filename
    # Dpath = '/Volumes/BU_extra/CryoSat/Antarctica+/2019_compressed/'
    Dpath = '/home/hds/DATA2/Ant+/DATA/DTU/'
    Dfiles = ant_plus.add_file('DTU',Dpath,'%Y/%m/')
    Dfiles.tvec = 'time_20_ku'
    Dfiles.vars = ['range_uncorr_imp_th70','range_uncorr_th70',
                   'range_uncorr_gauss','mss',
                  'otide_fes'] ### Tide to compare
    Dfiles.names = ['range_1_20_ku_DTU_imp_th70','range_1_20_ku_DTU_th70',
                    'range_1_20_ku_DTU_gauss','mean_sea_surf_20_ku_DTU15',
                   'ocean_tide_20_ku_DTU_fes']
    
    ### ISat files - date from filename
    # Ipath = '/Volumes/BU_extra/CryoSat/Antarctica+/IsardSat/2019_March_September_leads/'
    ILpath = '/raid6/userdata/hds/CSAO/IsardSat/Leads/'
    ILfiles = ant_plus.add_file('IsardSat_Leads',ILpath,'%m/',name_blurb='data/')
    ILfiles.tvec = 'time_20_ku'
    ILfiles.vars = ['retracked_range_analytical_2step_20_ku',
                  'total_geo_corr_20_ku']
    ILfiles.names = ['range_1_20_ku_ISat_2step',
                   'atm_geo_corrections_sum_ISat']

    ### ISat files - date from filename
#     ISpath = '/Volumes/BU_extra/CryoSat/Antarctica+/IsardSat/2019_March_September_Sea_Ice/'
    ISpath = '/raid6/userdata/hds/CSAO/IsardSat/Sea_Ice/'
    ISfiles = ant_plus.add_file('IsardSat_Ice',ISpath,'%m/',name_blurb='data/')
    ISfiles.tvec = 'time_20_ku'
    ISfiles.vars = [#'range_20_ku',
    #                 'retracked_range_threshold_20_ku',
    #                 'retracked_range_OCOG_20_ku',
                    'retracked_range_analytical_SWH_MSSfixed_20_ku',
                  'total_geo_corr_20_ku']
    ISfiles.names = [#'range_1_20_ku_ISat_Ice',
    #                  'range_1_20_ku_ISat_th_Ice',
    #                  'range_1_20_ku_ISat_OCOG_Ice',
                     'range_1_20_ku_ISat_SWH_MSSfixed_Ice',
                   'atm_geo_corrections_sum_ISat_Ice']
    ##### variable names that are loaded but we don't want to save
    var_remove = [
        'range_MSSL_Dcor',
        'range_MSSL_Scor',
        'range_20_ku_MSSL',
    'flag_floes_LEGOS_GPOD','flag_leads_LEGOS_GPOD','flag_open_ocean_LEGOS_GPOD',
    ]

    # All_files = [Dfiles,Cfiles,Ifiles,MDfiles,MSfiles]
#     All_files = [Dfiles,Cfiles,Mfiles,MDfiles,MSfiles]
    All_files = [Cfiles,Mfiles,MDfiles,MSfiles,Dfiles,ILfiles,ISfiles]
    # All_files = [Dfiles]
    for Af in All_files:
        Af.get_files_times(time,verbos=True)


    ### GPOD is special because it's in cycles sub folders
#     Gpath = '/Volumes/BU_extra/CryoSat/Antarctica+/GPOD_2020/'
    Gpath = '/home/hds/DATA/CSAO/LEGOS_GPOD/SAR/'
    Gfiles = ant_plus.add_file('GPOD',Gpath,'')
    mode = 'samosa_z'
    ### search through cycles first
    cycles = glob.glob(Gpath+'cycle*')
    cycles.sort()
    Gtime= []
    Gfiles.files=[]
    tstring = time.strftime('%Y%m%d')
    for cycle in cycles:
        print(cycle)
    #     temp_files = glob.glob(cycle+'/'+mode+'/fb_Z_SH/*.nc')
        temp_files = glob.glob(cycle+'/*.nc')
        temp_files.sort()
        ffound = 0
        for file in temp_files:
    #         temp_time=dt.datetime.strptime('201'+file.split(cycle)[1].split('_201')[1],'%Y%m%dT%H%M%S')
            if tstring in file:
                try:
                    ftime = dt.datetime.strptime(
                        tstring+file.split(cycle)[1].split(tstring)[1],
                        '%Y%m%dT%H%M%S_')
                    Gtime.append(ant_plus.dt2CSt(ftime))
                    Gfiles.files.append(file)
                    ffound+=1
                except ValueError:
                    pass
        if ffound >0: print('found '+str(ffound)+' files')
    # file = glob.glob(Gpath+mode+'/fb_Z_SH/fb_SRL_GPS_*_126_552_20191231T231346_20191231T231433.nc'
    # temp_files
    Gfiles.times = np.array(Gtime)
    Gfiles.tvec = 'time_tai_20hz'
    Gfiles.tshift = 1.0  -1e-7
    Gfiles.vars = [
        'all_corrs_ocean_20hz','mss_dtu_15_20hz','mss_dtu_21_20hz',
        'asa_raw_20hz',
        'ila_interp_20hz','sla_interp_20hz',
        'ila_smooth_20hz','sla_smooth_20hz',
        'flag_floes_meas_valid_20hz','flag_leads_meas_valid_20hz','flag_open_ocean_meas_valid_20hz',
                    ]
    Gfiles.names = ['atm_geo_corrections_sum_LEGOS_GPOD',
        'mean_sea_surf_20_ku_LEGOS_GPOD_DTU15','mean_sea_surf_20_ku_LEGOS_GPOD_DTU21',
        'ssha_20_ku_LEGOS_GPOD_All',
        'ssha_interp_20_ku_LEGOS_GPOD_Ice','ssha_interp_20_ku_LEGOS_GPOD_Lead',
        'ssha_smooth_20_ku_LEGOS_GPOD_Ice','ssha_smooth_20_ku_LEGOS_GPOD_Lead',
        'flag_floes_LEGOS_GPOD','flag_leads_LEGOS_GPOD','flag_open_ocean_LEGOS_GPOD',
                    ]

    All_files.append(Gfiles)

    ## ice conc too

#     ICE = ant_plus.NSIDC_nt('/Users/H/PREMELT/Ant+/DATA/NSIDC')
    ICE = ant_plus.NSIDC_nt('/home/hds/DATA2/NSIDC/nasa_daily/anto/')
#     time = dt.datetime(2019,9,25)
    gridIce = ICE.get_aice([time])

    Ilons = np.arange(-50,-35,0.1)
    Ilats = np.arange(-75,-60,0.1)
    testICE,iceinterp = GIce.GS2track(gridIce,Ilons,Ilats,save_array=True)


    ### find some data

    ### this variable name is a L2 variable that we check to see if the 
    ### data exists before loading
    ### if the variable is empty in a file (say no leads or sea ice) then we
    ### skip that file and keep searching
    f_check = 'range_1_20_ku'
    ### point in the L2 file list we start at
    Ltp = 0
    # Ltp =18
    # Ltp +=1
    ### no. of tracks to merge (make larger than 1 if you want to merge)
#     looklim = 1
#     fnum = 0
    looking =True
    tolerance = 100### path start in seconds

    while looking:
        ftime = SARtime[Ltp]
        ### here make the virtual L2 file
        track = ant_plus.CS2_track(SARfiles[Ltp],add_attr=True,
                                   lon_lims=lons,lat_lims=lats,count_lim=10)
        if track.use:
            new_files = np.array([Af.get_tp(ftime,tolerance,verbos=True)
                              for Af in All_files])
        else:
            Ltp += 1
            if (time+tw < ant_plus.CSt2dt(ftime)): 
                looking = False
                print('Next day: time limit')
            if Ltp == Lmax: 
                looking = False                    
                print('Next day: file limit')
            continue
        ### set L2 file
        if new_files.any():
            ### this is the variable from the l2 file we are interested in
            track.add_vars(L2_vars)
            track.add_vars_1Hz(L2_vars_1Hz)
            track.print_time() 
            print('L2 found: ',Ltp)
            ### add the gridded
            track.add_data_lat(track.lat_20_ku,
                               GIce.GS2track(iceinterp,track.lon_20_ku,track.lat_20_ku),
                               'NSIDC_nasa',verbos=True)
            track.NSIDC_nasa_attr = {
                'Ice_concentration':'NSIDC Nasa Team Ice Concentration',
                'Units':'Ice fraction',
            }
            ### loop the new files
            for load,Af in zip(new_files,All_files):
                if load:
                    Af.load_nc()
                    print(Af.name+' loading')
                    #### [list of data arrays] netCDF variables
                    vars_to_load = [Af.f_nc[vname][Af.F_pt] for vname in Af.vars]
                    attr_to_load = [ant_plus.add_ncattr(Af.f_nc[vname]) for vname in Af.vars]
                    ### shove them in the L2 file (virtual)
                    track.add_data_time(Af.f_nc[Af.tvec][Af.F_pt],
                           vars_to_load, Af.names, verbos=v_tload,
                            filt_less = -1e7, filt_great = 1e10,
                            nc_attr = attr_to_load)
                    Af.load_clean()
                else:
                    attr_to_load = [{'Warning':'Data empty'}]*len(Af.names)
                    track.add_data_time(np.array([0]),
                                [np.array([np.nan]) for i in range(len(Af.names))],
                                Af.names,nc_attr = attr_to_load)
            Ltp += 1
            ### masking options SAR
#             mask = track.flag_surf_type_class_20_ku.data[:] == 256
#             mask[track.flag_surf_type_class_20_ku.data[:] == 64] = True 
#             # mask = track.MSSL_surf_type.data[:] == 3
#             track.make_mask([],man_list = mask ,verbos=v_tmask)

            #### sum L2 corrections 
            track.geo_corrections_sum_L2 = ant_plus.add_geo_cor(track)
            track.geo_corrections_sum_L2_attr = {
                            'Description':'Summed L2 geo corrections',
                            'Original_variables':' , '.join([
                                'pole_tide_to_20_ku',
                                 'solid_earth_tide_to_20_ku',
                                 'load_tide_to_20_ku',
                                 'ocean_tide_to_20_ku',
                                 'ocean_tide_eq_to_20_ku']),
                            'coordinates'  :  'lon_poca_20_ku lat_poca_20_ku',
                            'Units':'m',
                        }
            track.vars.append('geo_corrections_sum_L2')
            track.atm_corrections_sum_L2 = ant_plus.add_atm_cor(track)
            track.atm_corrections_sum_L2_attr = {
                        'Description':'Summed L2 atmospheric corrections',
                        'Original_variables':' , '.join([
                            'mod_wet_tropo_cor_to_20_ku',
                             'mod_dry_tropo_cor_to_20_ku',
                             'iono_cor_gim_to_20_ku',
                             'inv_bar_cor_to_20_ku',
                             'hf_fluct_total_cor_to_20_ku']),
                        'coordinates'  :  'lon_poca_20_ku lat_poca_20_ku',
                        'Units':'m',
                    }
            track.vars.append('atm_corrections_sum_L2')
        

            #### HERE ARE THE OPTIONAL ANOMALIES
            #### 1 CLS - calculate anomaly
            #### 2 MSSL - Add Mfiles range_20_ku, then the anom is the difference
            track.range_1_20_ku_MSSL_D = track.range_20_ku_MSSL*1e-3+track.range_MSSL_Dcor
            track.range_1_20_ku_MSSL_S = track.range_20_ku_MSSL*1e-3+track.range_MSSL_Scor
            track.vars.append('range_1_20_ku_MSSL_D')
            track.vars.append('range_1_20_ku_MSSL_S')
            track.range_1_20_ku_MSSL_D_attr = track.range_MSSL_Dcor_attr
            track.range_1_20_ku_MSSL_S_attr = track.range_MSSL_Scor_attr
    
            #### 3 DTU - convert heights back to ranges
            track.range_1_20_ku_DTU_gauss = track.alt_20_ku - track.range_1_20_ku_DTU_gauss
            track.range_1_20_ku_DTU_th70 = track.alt_20_ku - track.range_1_20_ku_DTU_th70
            track.range_1_20_ku_DTU_imp_th70 = track.alt_20_ku - track.range_1_20_ku_DTU_imp_th70
            #### Tide anomaly
            track.ocean_tide_20_ku_ANOM_DTU_fes = (track.ocean_tide_20_ku_DTU_fes 
                                    - track.ocean_tide_eq_to_20_ku - track.ocean_tide_to_20_ku)
            track.ocean_tide_20_ku_ANOM_DTU_fes_attr = {
                'Description':'Anomaly from L2 ocean_tide_eq+ocean_tide to DTU_fes ocean tide',
                'coordinates'  :  'lon_poca_20_ku lat_poca_20_ku',
                'Units':'m',
                    }
            track.vars.append('ocean_tide_20_ku_ANOM_DTU_fes')
            #### 4 LEGOS GPOD - asa (ssha) back to range
            track.range_1_20_ku_LEGOS_GPOD_All = (track.alt_20_ku - 
                                           (track.ssha_20_ku_LEGOS_GPOD_All 
                                          + track.mean_sea_surf_20_ku_LEGOS_GPOD_DTU15
                                          + track.atm_geo_corrections_sum_LEGOS_GPOD))
            track.range_1_20_ku_LEGOS_GPOD_All_attr = track.ssha_20_ku_LEGOS_GPOD_All_attr
            track.range_1_20_ku_LEGOS_GPOD_All_attr['MERGE_NOTE'] = 'Converted back to range_1_20_ku equivalent \
                    using eqn: range = alt_20_ku - (ssha_20_ku_LEGOS_GPOD_All + mean_sea_surf_20_ku_LEGOS_GPOD_DTU15 \
                    + atm_geo_corrections_sum_LEGOS_GPOD)'
            track.vars.append('range_1_20_ku_LEGOS_GPOD_All')

            #### Range Anomalies
            p_vars = [p for p in track.vars if 'ISat' in p and 'range' in p]
            p_vars.extend([p for p in track.vars if 'DTU' in p and 'range' in p])
            p_vars.extend(['range_1_20_ku_MSSL_D','range_1_20_ku_MSSL_S',
                           'range_1_20_ku_LEGOS_GPOD_All','range_1_20_ku_CLS'])
            for v in p_vars:
                v_ANOM = 'range_1_20_ku_ANOM'+v.split('20_ku')[1]
    #             print(v_ANOM)
                v_ANOM_attr = 'range_1_20_ku_ANOM'+v.split('20_ku')[1]+'_attr'
                x = getattr(track,v)
                y = x - track.range_1_20_ku
                attr_dict = {
            'Description':'Anomaly from L2 range to '+v.split('20_ku')[1]+' retracker range',
            'coordinates'  :  'lon_poca_20_ku lat_poca_20_ku',
            'Units':'m',
                        }
                setattr(track,v_ANOM,y)
                setattr(track,v_ANOM_attr,attr_dict)
                track.vars.append(v_ANOM)
            #### MSS anomalies
            p_vars = [p for p in track.vars if 'mean_sea' in p and 'sea_ice' not in p]
            for v in p_vars:
                v_ANOM = 'mean_sea_surf_20_ku_ANOM'+v.split('20_ku')[1]
    #             print(v_ANOM)
                v_ANOM_attr = 'mean_sea_surf_20_ku_ANOM'+v.split('20_ku')[1]+'_attr'
                x = getattr(track,v)
                y = x - track.mean_sea_surf_sea_ice_20_ku
                attr_dict = {
            'Description':'Anomaly from L2 MSS to '+v.split('20_ku')[1]+' MSS',
            'coordinates'  :  'lon_poca_20_ku lat_poca_20_ku',
            'Units':'m',
                        }
                setattr(track,v_ANOM,y)
                setattr(track,v_ANOM_attr,attr_dict)
                track.vars.append(v_ANOM)
            #### 4 IsardSat - simple range comparison
            #### Other Anomalies
            track.atm_geo_corrections_sum_ANOM_ISat = (track.atm_geo_corrections_sum_ISat 
                                    - track.geo_corrections_sum_L2 - track.atm_corrections_sum_L2)
            track.atm_geo_corrections_sum_ANOM_ISat_attr = {
                'Description':'Anomaly from L2 atm/geo corrections to IsardSat atm/geo corrections',
                'coordinates'  :  'lon_poca_20_ku lat_poca_20_ku',
                'Units':'m',
                        }
            track.vars.append('atm_geo_corrections_sum_ANOM_ISat')
            #### GPOD corrections
            track.ssha_interp_20_ku_ANOM_LEGOS_GPOD = track.ssha_interp_20_ku_LEGOS_GPOD_Lead - track.ssha_interp_20_ku
            track.ssha_interp_20_ku_ANOM_LEGOS_GPOD_attr = {
                'Description':'Anomaly from L2 ssha_interp to LEGOS GPOD sla_smooth',
                'coordinates'  :  'lon_poca_20_ku lat_poca_20_ku',
                'Units':'m',
                        }
            track.vars.append('ssha_interp_20_ku_ANOM_LEGOS_GPOD')
            track.atm_geo_corrections_sum_ANOM_LEGOS_GPOD = (track.atm_geo_corrections_sum_LEGOS_GPOD
                                    - track.geo_corrections_sum_L2 - track.atm_corrections_sum_L2)
            track.atm_geo_corrections_sum_ANOM_LEGOS_GPOD_attr = {
                'Description':'Anomaly from L2 atm/geo corrections to LEGOS_GPOD atm/geo corrections',
                'coordinates'  :  'lon_poca_20_ku lat_poca_20_ku',
                'Units':'m',
                    }
            track.vars.append('atm_geo_corrections_sum_ANOM_LEGOS_GPOD')
        
            ### fixing the LEGOS GPOD surface flag
            temp_mask = np.zeros([track.n_u],dtype=int)
            temp_mask[track.flag_open_ocean_LEGOS_GPOD == 1.0] = 1
            temp_mask[track.flag_floes_LEGOS_GPOD == 1.0] = 2
            temp_mask[track.flag_leads_LEGOS_GPOD == 1.0] = 3
            track.flag_surf_type_20_ku_LEGOS_GPOD = temp_mask
            track.flag_surf_type_20_ku_LEGOS_GPOD_attr = {
    'coordinates': 'lon_20_ku lat_20_ku',
    'flag_values': '0b, 1b, 2b , 3b',
    'flag_meanings': 'unknown/mixed ocean sea_ice lead',
    'long_name': 'surface type flag',
    'comment': 'A 4-state surface type mask for Cryosat2 data from the LEGOS GPOD along track data.Flag file converted from the three original supplied flag files: flag_floes_meas_valid_20hz,flag_leads_meas_valid_20hz,flag_open_ocean_meas_valid_20h'
            }
            track.vars.append('flag_surf_type_20_ku_LEGOS_GPOD')

            #### remove unwanted vars from the list
            track.vars = [v for v in track.vars if v not in var_remove]

            if COMBINE_FILES and read_comb_attr:
                ### set attr
                for v in Combine_vars:
                    variable = getattr(Comb_list,v)
                    variable.attr = getattr(track,v+'_attr')
                read_comb_attr = False

        #### save the track
            if MERGE_TRACKS:
                track.check_attr()
                track.save_nc(save_root,l1b_root,merge_list=[Af.name for Af in All_files],verbos = v_save )
            if COMBINE_FILES:
#                 track.make_mask([],man_list=(track.flag_surf_type_class_20_ku == Comb_ST))
                track.make_mask([],man_list=Comb_ST(track))
                ### fill up the track list
                data_list = [getattr(track,v) for v in Combine_vars]
                Comb_list.addtrack(track,data_list)
                
            if PLOT_RANGE_ANOM:
                track.make_mask([],man_list=track.flag_surf_type_class_20_ku >= 64)
                f = plt.figure(figsize=[20,5])
                aplt = 0.6
                mksz = 0.8
                pstlye = '-'

                ax = plt.subplot2grid([1,4],[0,0])
                x,y = m(track.lon_20_ku,track.lat_20_ku)
    #             m.scatter(x,y,c=track.range_1_20_ku)
                m.scatter(x,y,c=track.alt_20_ku - track.range_1_20_ku)
                m.colorbar(location='bottom')
                m.drawcoastlines()
                # m.contour([GIce.xpts,GIce.ypts],gridIce[0],[15,85],cmap = 'k')
                m.contour(GIce.xpts,GIce.ypts,gridIce[0],[0.15,0.80],cmap='Blues_r',vmax=4)
                # plt.xlim(map_range(x,1.5))
                # plt.ylim(map_range(y,1.5))
                Icedesc = np.sum(np.diff(track.NSIDC_nasa))<0
                conc_lim1 = 0.15
                conc_lim2 = 0.85
                if Icedesc:
                    icelim=[
                            track.time_20_ku[np.argmax(track.NSIDC_nasa<conc_lim1)]-ftime,
                            track.time_20_ku[np.argmax(track.NSIDC_nasa<conc_lim2)]-ftime,
                            track.time_20_ku[0]-ftime]
                else:
                    icelim=[track.time_20_ku[np.argmax(track.NSIDC_nasa>conc_lim1)]-ftime,
                            track.time_20_ku[np.argmax(track.NSIDC_nasa>conc_lim2)]-ftime,
                            track.time_20_ku[track.n_u-1]-ftime]
                if (track.NSIDC_nasa<conc_lim2).all() and Icedesc:
                    icelim[1] = track.time_20_ku[0]-ftime
                elif (track.NSIDC_nasa<conc_lim2).all():
                    icelim[1] = track.time_20_ku[track.n_u-1]-ftime


                ax = plt.subplot2grid([1,4],[0,1],colspan=3)
                plt.axvspan(icelim[0],icelim[1],facecolor='k',alpha=0.1)
                plt.axvspan(icelim[1],icelim[2],facecolor='k',alpha=0.2)
                p_vars = [p for p in track.vars if 'ANOM' in p]
                p_vars = [p for p in     p_vars if 'range' in p]
                p_no = len(p_vars)
                colors = pl.cm.jet(np.linspace(0,0.9,p_no))
                ax.hlines([0.0],0,track.time_20_ku[-1]-ftime)
                for pn,v in enumerate(p_vars):
                    x = copy.copy(getattr(track,v))
                    x[track.mask==0] = np.nan
                    x[x.mask] = np.nan
                    plt.plot(track.time_20_ku-ftime,x,pstlye,alpha=aplt,markersize=mksz,
                        color=colors[pn])
#                 plt.ylim([-1.5,1.5])
                plt.ylim([-4,4])
                plt.ylabel('Range anomalies')
                plt.xlabel('Track time')
                xplot_r = [0,track.time_20_ku[track.n_u-1]-ftime]

                plt.legend(p_vars)
                # plt.ylim([10,pmax])
    #             Tstring = ant_plus.CSt2dt(ftime,string=True)
                Tstring = ant_plus.CSt2dt(ftime).strftime('%Y-%m-%dT%H:%M:%S')
                plt.title(Tstring+'_Anomalies_to_L2')

                plot_name = 'CSAO_'+Tstring.replace(':','_',100)+'_range_ANOM.png'
                print('Saving: '+plot_name)
                f.savefig(fig_dir+plot_name,bbox_inches='tight')
                plt.close()
                
            if PLOT_EXTRA_ANOM:
                track.make_mask([],man_list=track.flag_surf_type_class_20_ku >= 64)
                f = plt.figure(figsize=[20,5])
                aplt = 0.6
                mksz = 0.8
                pstlye = '-'

                ax = plt.subplot2grid([1,4],[0,0])
                x,y = m(track.lon_20_ku,track.lat_20_ku)
    #             m.scatter(x,y,c=track.range_1_20_ku)
                m.scatter(x,y,c=track.alt_20_ku - track.range_1_20_ku)
                m.colorbar(location='bottom')
                m.drawcoastlines()
                # m.contour([GIce.xpts,GIce.ypts],gridIce[0],[15,85],cmap = 'k')
                m.contour(GIce.xpts,GIce.ypts,gridIce[0],[0.15,0.80],cmap='Blues_r',vmax=4)
                # plt.xlim(map_range(x,1.5))
                # plt.ylim(map_range(y,1.5))
                Icedesc = np.sum(np.diff(track.NSIDC_nasa))<0
                conc_lim1 = 0.15
                conc_lim2 = 0.85
                if Icedesc:
                    icelim=[
                            track.time_20_ku[np.argmax(track.NSIDC_nasa<conc_lim1)]-ftime,
                            track.time_20_ku[np.argmax(track.NSIDC_nasa<conc_lim2)]-ftime,
                            track.time_20_ku[0]-ftime]
                else:
                    icelim=[track.time_20_ku[np.argmax(track.NSIDC_nasa>conc_lim1)]-ftime,
                            track.time_20_ku[np.argmax(track.NSIDC_nasa>conc_lim2)]-ftime,
                            track.time_20_ku[track.n_u-1]-ftime]
                if (track.NSIDC_nasa<conc_lim2).all() and Icedesc:
                    icelim[1] = track.time_20_ku[0]-ftime
                elif (track.NSIDC_nasa<conc_lim2).all():
                    icelim[1] = track.time_20_ku[track.n_u-1]-ftime


                ax = plt.subplot2grid([1,4],[0,1],colspan=3)
                plt.axvspan(icelim[0],icelim[1],facecolor='k',alpha=0.1)
                plt.axvspan(icelim[1],icelim[2],facecolor='k',alpha=0.2)
                p_vars = [p for p in track.vars if 'ANOM' in p]
                p_vars = [p for p in     p_vars if 'range' not in p]
                p_no = len(p_vars)
                colors = pl.cm.jet(np.linspace(0,0.9,p_no))
                ax.hlines([0.0],0,track.time_20_ku[-1]-ftime)
                for pn,v in enumerate(p_vars):
                    x = copy.copy(getattr(track,v))
                    x[track.mask==0] = np.nan
                    x[x.mask] = np.nan
                    plt.plot(track.time_20_ku-ftime,x,pstlye,alpha=aplt,markersize=mksz,
                        color=colors[pn])
                plt.ylim([-1.5,1.5])
                plt.title(Tstring+'_Anomalies_to_L2')
                plt.xlabel('Track time')
                xplot_r = [0,track.time_20_ku[track.n_u-1]-ftime]

                plt.legend(p_vars)
                # plt.ylim([10,pmax])
    #             Tstring = ant_plus.CSt2dt(ftime,string=True)
                Tstring = ant_plus.CSt2dt(ftime).strftime('%Y-%m-%dT%H:%M:%S')
                plt.title(Tstring+'_Range')

                plot_name = 'CSAO_'+Tstring.replace(':','_',100)+'_extra_ANOM.png'
                print('Saving: '+plot_name)
                f.savefig(fig_dir+plot_name,bbox_inches='tight')
                plt.close()


            if PLOT_RANGE:
                f = plt.figure(figsize=[20,5])
                aplt = 0.6
                mksz = 0.8
                pstlye = '-'

                ax = plt.subplot2grid([1,4],[0,0])
                x,y = m(track.lon_20_ku,track.lat_20_ku)
                m.scatter(x,y,c=track.alt_20_ku - track.range_1_20_ku)
                m.colorbar(location='bottom')
                m.drawcoastlines()
                # m.contour([GIce.xpts,GIce.ypts],gridIce[0],[15,85],cmap = 'k')
                m.contour(GIce.xpts,GIce.ypts,gridIce[0],[0.15,0.80],cmap='Blues_r',vmax=4)
                # plt.xlim(map_range(x,1.5))
                # plt.ylim(map_range(y,1.5))
                Icedesc = np.sum(np.diff(track.NSIDC_nasa))<0
                conc_lim1 = 0.15
                conc_lim2 = 0.85
                if Icedesc:
                    icelim=[
                            track.time_20_ku[np.argmax(track.NSIDC_nasa<conc_lim1)]-ftime,
                            track.time_20_ku[np.argmax(track.NSIDC_nasa<conc_lim2)]-ftime,
                            track.time_20_ku[0]-ftime]
                else:
                    icelim=[track.time_20_ku[np.argmax(track.NSIDC_nasa>conc_lim1)]-ftime,
                            track.time_20_ku[np.argmax(track.NSIDC_nasa>conc_lim2)]-ftime,
                            track.time_20_ku[track.n_u-1]-ftime]
                if (track.NSIDC_nasa<conc_lim2).all() and Icedesc:
                    icelim[1] = track.time_20_ku[0]-ftime
                elif (track.NSIDC_nasa<conc_lim2).all():
                    icelim[1] = track.time_20_ku[track.n_u-1]-ftime


                ax = plt.subplot2grid([1,4],[0,1],colspan=3)
                plt.axvspan(icelim[0],icelim[1],facecolor='k',alpha=0.1)
                plt.axvspan(icelim[1],icelim[2],facecolor='k',alpha=0.2)
                p_vars = [p for p in track.vars if 'range' in p]
                p_vars = [p for p in p_vars if 'ANOM' not in p]
                p_no = len(p_vars)
                colors = pl.cm.jet(np.linspace(0,0.9,p_no))
                pmin = [-10.0]
                pmax = [ 10.0]
                for pn,v in enumerate(p_vars):
                    x = copy.copy(getattr(track,v))
                    x[track.mask==0] = np.nan
                    x[x.mask] = np.nan
                    x = track.alt_20_ku - x
                    pmin.append(np.nanpercentile(x.compressed(),[ 3.1]))
                    pmax.append(np.nanpercentile(x.compressed(),[ 97]))
                    plt.plot(track.time_20_ku-ftime,x,pstlye,alpha=aplt,markersize=mksz,
                        color=colors[pn])
                plt.ylim([np.nanmin(pmin),np.nanmax(pmax)])
                # plt.ylim([10,pmax])
                plt.ylabel('Alt - Range')
                plt.xlabel('Track time')
                xplot_r = [0,track.time_20_ku[track.n_u-1]-ftime]

                plt.legend(p_vars)
                # plt.ylim([10,pmax])
                Tstring = ant_plus.CSt2dt(ftime).strftime('%Y-%m-%dT%H:%M:%S')
                plt.title(Tstring+'_Range')

                plot_name = 'CSAO_'+Tstring.replace(':','_',100)+'_range.png'
                print('Saving: '+plot_name)
                f.savefig(fig_dir+plot_name,bbox_inches='tight')
                plt.close()
            if PLOT_INTERP:
                track.make_mask([],man_list=track.flag_surf_type_class_20_ku >= 32)
                f = plt.figure(figsize=[20,5])
                aplt = 0.6
                mksz = 0.8
                pstlye = '-'

                ax = plt.subplot2grid([1,4],[0,0])
                x,y = m(track.lon_20_ku,track.lat_20_ku)
                m.scatter(x,y,c=track.ssha_interp_20_ku)
                m.colorbar(location='bottom')
                m.drawcoastlines()
                # m.contour([GIce.xpts,GIce.ypts],gridIce[0],[15,85],cmap = 'k')
                m.contour(GIce.xpts,GIce.ypts,gridIce[0],[0.15,0.80],cmap='Blues_r',vmax=4)
                # plt.xlim(map_range(x,1.5))
                # plt.ylim(map_range(y,1.5))
                Icedesc = np.sum(np.diff(track.NSIDC_nasa))<0
                conc_lim1 = 0.15
                conc_lim2 = 0.85
                if Icedesc:
                    icelim=[
                            track.time_20_ku[np.argmax(track.NSIDC_nasa<conc_lim1)]-ftime,
                            track.time_20_ku[np.argmax(track.NSIDC_nasa<conc_lim2)]-ftime,
                            track.time_20_ku[0]-ftime]
                else:
                    icelim=[track.time_20_ku[np.argmax(track.NSIDC_nasa>conc_lim1)]-ftime,
                            track.time_20_ku[np.argmax(track.NSIDC_nasa>conc_lim2)]-ftime,
                            track.time_20_ku[track.n_u-1]-ftime]
                if (track.NSIDC_nasa<conc_lim2).all() and Icedesc:
                    icelim[1] = track.time_20_ku[0]-ftime
                elif (track.NSIDC_nasa<conc_lim2).all():
                    icelim[1] = track.time_20_ku[track.n_u-1]-ftime


                ax = plt.subplot2grid([1,4],[0,1],colspan=3)
                plt.axvspan(icelim[0],icelim[1],facecolor='k',alpha=0.1)
                plt.axvspan(icelim[1],icelim[2],facecolor='k',alpha=0.2)
                p_vars = ['ssha_20_ku']+[p for p in track.vars if 'interp' in p]
                p_vars = p_vars+[p for p in track.vars if 'smooth' in p]
                p_no = len(p_vars)
                colors = pl.cm.jet(np.linspace(0,0.9,p_no))
                pmin = [np.nan]
                pmax = [np.nan]
                for pn,v in enumerate(p_vars):
                    x = copy.copy(getattr(track,v))
                    x[track.mask==0] = np.nan
                    x[x.mask] = np.nan
                    if v =='ssha_20_ku':
                        x[track.flag_surf_type_class_20_ku!=128] = np.nan
                        x[x<0.0] = np.nan
                    pmin.append(np.nanpercentile(x.compressed(),[ 3.1]))
                    pmax.append(np.nanpercentile(x.compressed(),[ 97]))
                    plt.plot(track.time_20_ku-ftime,x,pstlye,alpha=aplt,markersize=mksz,
                        color=colors[pn])
                
                plt.ylim([np.maximum(np.nanmin(pmin),-2.0)
                         ,np.minimum(np.nanmax(pmax),2.0)])
                # plt.ylim([10,pmax])
                plt.ylabel('Interpolated sea level anomaly')
                plt.xlabel('Track time')
                xplot_r = [0,track.time_20_ku[track.n_u-1]-ftime]

                plt.legend(p_vars)
                # plt.ylim([10,pmax])
                Tstring = ant_plus.CSt2dt(ftime).strftime('%Y-%m-%dT%H:%M:%S')
                plt.title(Tstring+'_Range')

                plot_name = 'CSAO_'+Tstring.replace(':','_',100)+'_interp_ssha_ice.png'
                print('Saving: '+plot_name)
                f.savefig(fig_dir+plot_name,bbox_inches='tight')
                plt.close()
            if check_single_track:
                track.make_mask([],man_list=track.flag_surf_type_class_20_ku >= 64)
                ### loop the variables
                for v in track.vars:
                    if check_filter in v:
                        print('-----')
                        print(v)
                        print('-----')
                        if hasattr(track,v+'_attr'):
                            attr = getattr(track,v+'_attr')
                            for key, value in attr.items():
                                print(key, ' : ', value)
                        x = copy.copy(getattr(track,v))
                        if type(x.mask) == np.ndarray:
                            x[x.mask] = np.nan
                        print('Distribution of '+v+' = '+' , '.join('{:.3}'.format(pc) 
                                     for pc in np.nanpercentile(x[track.mask],[25,50,75])))
                looking = False
        else:
            Ltp +=1
    #         f_nc.close()
        if (time+tw < ant_plus.CSt2dt(ftime)): 
            looking = False
            print('Next day: time limit')
        if Ltp == Lmax: 
            looking = False                    
            print('Next day: file limit')

if COMBINE_FILES:
    #### Save the combined file
    Comb_list.check_attr()
    Comb_list.save_nc(merge_dir,Combine_name,merge_list=[Af.name for Af in All_files],add_attr=True,verbos = True )