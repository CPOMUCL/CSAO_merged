
import numpy as np
import datetime as dt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from dateutil.relativedelta import relativedelta
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
# from mpl_toolkits.basemap import Basemap
import glob
from netCDF4 import Dataset
import copy
import gc

import os
import sys
# sys.path.insert(0, '/Users/h/Github/geo_data_group/')
import grid_set as gs

import ant_plus
# from imp import reload

#$### map projection
m = ccrs.SouthPolarStereo()
G= gs.grid_set(m)

G.load_grid('100km_for_CSAO.npz')

# ax = f.add_subplot(1,1,1,projection=m)
# ax.set_extent([-180, 180, -90, -50], ccrs.PlateCarree())
# G.set_grid_dxy(100e3,100e3,ax)
# plt.close()


dsearch = dt.datetime(2019,1,1)

Combdir = '/home/hds/DATA/CSAO/MERGED/'

###  '/Users/h/PREMELT/Ant+/Comb_files_for_Stine/DOT_SLA_whole_anto_2019-09-01T00-03_2019-09-30T23-42.nc'
SARfile = glob.glob(Combdir+'DOT_SLA_whole_anto_'+dsearch.strftime('%Y-%m')+'*.nc')[0]
print(SARfile)
trackSAR=ant_plus.CS2_track(SARfile,add_attr=True)


SINfile = glob.glob(Combdir+'ST_and_range_ANOM_all_SARIN_'+dsearch.strftime('%Y-%m')+'*.nc')[0]
print(SINfile)
trackSIN=ant_plus.CS2_track(SINfile,add_attr=True)

trackSIN.list_vars()#['ANOM','flag',])

LRMfile = glob.glob(Combdir+'ST_and_range_ANOM_all_LRM_'+dsearch.strftime('%Y-%m')+'*.nc')[0]
print(LRMfile)
trackLRM=ant_plus.CS2_track(LRMfile,add_attr=True)


trackLRM.list_vars()#['ANOM','flag',])

load_list = trackLRM.list_vars(['ANOM','flag','height','geoid','mean_sea_surf_sea_ice_20_ku'])
# load_list = track.list_vars(['mean_sea_surf_sea_ice_20_ku'])
trackLRM.add_vars(load_list)

load_list = trackSAR.list_vars(['ANOM','flag','height','geoid','mean_sea_surf_sea_ice_20_ku'])
# load_list = track.list_vars(['mean_sea_surf_sea_ice_20_ku'])
trackSAR.add_vars(load_list)

load_list = trackSIN.list_vars(['ANOM','flag','height','geoid','mean_sea_surf_sea_ice_20_ku'])
# load_list = track.list_vars(['mean_sea_surf_sea_ice_20_ku'])
trackSIN.add_vars(load_list)


trackLRM.range = '-'.join([
    ant_plus.CSt2dt(trackSAR.time_20_ku[0]).strftime('%Y%m%d'),
    ant_plus.CSt2dt(trackSAR.time_20_ku[-1]).strftime('%Y%m%d')])
trackSAR.range = '-'.join([
    ant_plus.CSt2dt(trackSAR.time_20_ku[0]).strftime('%Y%m%d'),
    ant_plus.CSt2dt(trackSAR.time_20_ku[-1]).strftime('%Y%m%d')])
trackSIN.range = '-'.join([
    ant_plus.CSt2dt(trackSIN.time_20_ku[0]).strftime('%Y%m%d'),
    ant_plus.CSt2dt(trackSIN.time_20_ku[-1]).strftime('%Y%m%d')])

### now a big loop to bin all the cases we want
### this generates the L2 plus range anomalies
trackinfo = [
    {"Name":'LRM',"track":trackLRM,"mask_no":2,
     "Cases":['Default_L2','height_3'],
     "Anomalies":[['None']]},
    {"Name":'SAR_O',"track":trackSAR,"mask_no":64,
     "Cases":['Default_L2',
 'range_1_20_ku_ANOM_ISat_2step',
 'range_1_20_ku_ANOM_ISat_SWH_MSSfixed_Ice',
 'range_1_20_ku_ANOM_LEGOS_GPOD_All' , 
             ],
     "Anomalies":[['None'],
 ['atm_geo_corrections_sum_ANOM_ISat',],
 ['mean_sea_surf_20_ku_ANOM_LEGOS_GPOD_DTU21' ], 
 ['atm_geo_corrections_sum_ANOM_ISat',
     'mean_sea_surf_20_ku_ANOM_LEGOS_GPOD_DTU21'], 
                 ]},
    {"Name":'SAR_L',"track":trackSAR,"mask_no":256,
     "Cases":['Default_L2',
 'range_1_20_ku_ANOM_ISat_2step',
 'range_1_20_ku_ANOM_ISat_SWH_MSSfixed_Ice',
 'range_1_20_ku_ANOM_LEGOS_GPOD_All' , 
 'range_1_20_ku_ANOM_DTU_imp_th70',
 'range_1_20_ku_ANOM_DTU_gauss',
 'range_1_20_ku_ANOM_CLS',
             ],
     "Anomalies":[['None'],
 ['atm_geo_corrections_sum_ANOM_ISat',],
 ['mean_sea_surf_20_ku_ANOM_LEGOS_GPOD_DTU21'], 
 ['atm_geo_corrections_sum_ANOM_ISat',
     'mean_sea_surf_20_ku_ANOM_LEGOS_GPOD_DTU21'], 
                 ]},
    {"Name":'SARIN',"track":trackSIN,"mask_no":272,
     "Cases":['Default_L2',
     'range_1_20_ku_ANOM_LEGOS_GPOD_All' , 
             ],
     "Anomalies":[['None'],
 ['atm_geo_corrections_sum_ANOM_LEGOS_GPOD'],
 ['mean_sea_surf_20_ku_ANOM_LEGOS_GPOD_DTU21' ], 
 ['atm_geo_corrections_sum_ANOM_LEGOS_GPOD',
     'mean_sea_surf_20_ku_ANOM_LEGOS_GPOD_DTU21'], 
                 ]},
]

### check the variables (exist some aren't in certain months)
for TI in trackinfo:
    cases_correct = []
    for v in TI["Cases"]:
        if v in TI["track"].vars+['Default_L2','height_3']:
            cases_correct.append(v)
    TI["Cases"] = cases_correct


### saving everthing SLA
### map for each - put it on the class - Grid
### for each anomaly, make SLA map, save as SLA_map_(var_name)
G.var_list = []
for TI in trackinfo:
    print('---- '+TI['Name'])
    msk = TI['track'].flag_surf_type_class_20_ku == TI['mask_no']
#     pvars = [v for v in TI['track'].vars if 'ANOM' in v] + TI["Cases"]
    pvars =  TI["Cases"]
    for v in pvars:
        vprint = ''.join([vs for vs in v])
        if '20_ku_' in v: vprint = v.split('20_ku_')[1]
        ### make attribute with all the info
        attr = {"Options_added":v}
        print('      '+v)
#         if 'GPOD' not in v: continue
        x = TI['track'].height_1_20_ku - TI['track'].mean_sea_surf_sea_ice_20_ku
        if 'Default' in v:
            pass
        elif 'height_3' in v:
            x = TI['track'].height_3_20_ku - TI['track'].mean_sea_surf_sea_ice_20_ku
        else:
            x = x - getattr(TI['track'],v).copy()
        ### here we add the extra correction/mss anomalies
        for an in TI['Anomalies']:
            xA = x.copy()
            attr = {"Options_added":v}
            ### loop subtracting the anomalies
#             print(' '.join([a.split('ANOM')[1] for a in an]))
            for a in an:
                if a != 'None':
                    xA = xA - getattr(TI['track'],a).copy()
                    attr["Options_added"] = attr["Options_added"]+' '+a
            print(attr["Options_added"])
            if 'CLS' in v:
                #### set CLS mask
                if TI['mask_no']==64: CLSmsk = 1 ## ocean
                if TI['mask_no']==256: CLSmsk = 3 ## lead
                msk = TI["track"].flag_surf_type_20_ku_CLS == CLSmsk
            if 'GPOD' in v and 'SAR_' in TI['Name']: ### only for SAR mode
                #### set GPOD mask
                if TI['mask_no']==64: CLSmsk = 1 ## ocean
                if TI['mask_no']==256: CLSmsk = 3 ## lead
                msk = TI["track"].flag_surf_type_20_ku_LEGOS_GPOD == CLSmsk
            hard_min = -10.0 ### extras to remove crazy outliers
            hard_max =  10.0
            xA[xA.mask] = np.nan
            msk[xA>hard_max] = 0
            msk[xA<hard_min] = 0
            xA = xA[msk]
            #### only if we have any data left
#             if np.sum(msk)>10:
            try:
                plot_bin_temp  = G.bin_list(xA,#trackSIN.mean_sea_surf_20_ku_ANOM_LEGOS_GPOD_DTU15,
                                            TI['track'].lon_20_ku[msk],TI['track'].lat_20_ku[msk],xy_order=0)
                bin_var = '_'.join([TI['Name'],'SLA_bin_range',vprint]+
                        [a.split('ANOM')[1] for a in an if 'ANOM' in a])
                G.var_list.append(bin_var)
                print(bin_var)
                setattr(G,bin_var,plot_bin_temp.copy())
                setattr(G,bin_var+'_attr',attr)
            except ValueError:
                pass


### save netcdf
# nc_name = '/Users/h/PREMELT/Ant+/LRM-SAR-SIN_comb/All_modes_SLA_bin_'+trackinfo[0]['track'].range+'.nc'
nc_name = Combdir+'All_modes_SLA_bin_'+trackinfo[0]['track'].range+'.nc'
f_nc = Dataset(nc_name, mode='w')

# Create the dimensions of the file
# this is time_20_ku
f_nc.createDimension('x', G.m)
f_nc.createDimension('y', G.n)

# Copy the global attributes
Fattr = {
    'Description':'CSAO binned SLA',
    'CSAO_code_by':'HDBS Heorton',
}
f_nc.setncatts(Fattr)

# Create the variables in the file

## lon/lat and all
for v in ['lons','lats']+G.var_list:
    fill_var = getattr(G,v)
    f_nc.createVariable(v, fill_var.dtype, ('x','y'))
    # Fill it
    f_nc.variables[v][:] = fill_var
    if hasattr(G,v+'_attr'):
        fill_var_attr = getattr(G,v+'_attr')
        f_nc.variables[v].setncatts(fill_var_attr)

# Save the file
f_nc.close()
print(nc_name)
