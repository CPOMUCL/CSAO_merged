### imports
import numpy as np
import datetime as dt
from netCDF4 import Dataset
from dateutil.relativedelta import relativedelta
import glob
import os


        
# class _NO_DEFAULT:
#     def __repr__(self):return "<no default>"
# _NO_DEFAULT = _NO_DEFAULT()

def CSt2dt(t,string=False):
    t0 = dt.datetime(2000,1,1)
    if string:
        return (t0+relativedelta(seconds=t-36.0)).strftime('%Y-%m-%dT%H:%M:%S.%f')
    else:
        return t0+relativedelta(seconds=t-36.0)
def dt2CSt(t):
    t0 = dt.datetime(2000,1,1)
    return (t-t0).total_seconds() + 36.0

### new code for attributes from nc_files
def add_ncattr(var):
    """
    for the given var, this returns a dict of the attributes
    """
    attr_list = var.ncattrs()
    return {a:var.getncattr(a) for a in attr_list}
    

# give track details - either time or cycle no. etc
# load Bl-D L2 file
class CS2_track():
    def __init__(self,file,lon_lims=[-180,180],lat_lims=[-90,90],add_attr = False,count_lim = 2,surface_type=None):
        self.name = file
        self.swh_calc = False
        self.add_attr = add_attr
        try:
            data_nc = Dataset(file)
        except OSError:
            self.use = False
            return
        lats = data_nc.variables["lat_20_ku"]
        lons = data_nc.variables["lon_20_ku"]
        lon_correct = (lons[:] > lon_lims[0])&(lons[:] < lon_lims[1] )
        lat_correct = (lats[:] > lat_lims[0])&(lats[:] < lat_lims[1] )
        d_u = lon_correct&lat_correct
        if surface_type is not None:
            ##### surface typ is list of number to accept
            surf = data_nc.variables["flag_surf_type_class_20_ku"]
            s_u =  np.array([surf[:] == s for s in surface_type]).any(axis=0)
            d_u = d_u&s_u
            
        if np.sum(d_u)>count_lim:
            self.lon_20_ku = lons[d_u]
            self.lat_20_ku = lats[d_u]
            self.time_20_ku = data_nc.variables["time_20_ku"][d_u]
            if self.add_attr:
                self.time_20_ku_attr = add_ncattr(data_nc.variables["time_20_ku"])
                self.lat_20_ku_attr = add_ncattr(data_nc.variables["lat_20_ku"])
                self.lon_20_ku_attr = add_ncattr(data_nc.variables["lon_20_ku"])
            self.use = True
            self.vars = ['lat_20_ku','lon_20_ku']
            self.d_u = d_u
            self.n_u = np.sum(d_u)
            data_nc.close()
        else:
            self.use = False
            data_nc.close()
    def list_vars(self,v_filter=['']):
        data_nc = Dataset(self.name)
        out_list =  [v for v in data_nc.variables.keys()
                  for st in v_filter if st in v ]
        data_nc.close()
        return out_list
        
        

    def add_vars(self,var,verbos=False,change_name=None):
        if change_name is not None: self.name = change_name
        if type(var) == str:
            var = [var]
        data_nc = Dataset(self.name)
        for v in var:
            if v in self.vars:
                if verbos: print(v+' exists already, ignoring')
                continue
            ##### read var from nc_file
            vnew = data_nc[v] 
            ##### check var_shape
            try:
                dnew = vnew[self.d_u]
            except IndexError:
                if verbos: print(v+' not consistent with preset size, ignoring')
                continue
            ##### set_attribute
            setattr(self, v, dnew)
            self.vars.append(v)
            if self.add_attr:
                setattr(self, v+'_attr', add_ncattr(vnew))
        data_nc.close()
        
        

    def add_vars_1Hz(self,var,method = 'const'):
        if type(var) == str:
            var = [var]
        data_nc = Dataset(self.name)
        #### get the locating array
        loc_array = data_nc['ind_first_meas_20hz_01']
        for v in var:
            ### check new data shape for L1B
            vshape = data_nc[v].shape
            if np.shape(vshape)[0]==2:
                ### intialise 2d 20 hz array
                dnew = np.ma.masked_all([self.n_u,vshape[1]])
            else:
                ### intialise 20 hz array
                dnew = np.ma.masked_all([self.n_u])
            ##### read var from nc_file
            vnew = data_nc[v] 
            ##### shift from 1hz to 20hz
            if method == 'const':
                ne = loc_array[-1]
                for n in range(loc_array.shape[0]-1):
                    ns = loc_array[n]
                    ne = loc_array[n+1]
                    if np.shape(vshape)[0]==2:
                        dnew[ns:ne,:] = vnew[n,:]
                    else:
                        dnew[ns:ne] = vnew[n]
                n = loc_array.shape[0]-1
                ns = loc_array[-1]
                if np.shape(vshape)[0]==2:
                    try:
                        dnew[ns:ne,:] = vnew[n,:]
                    except IndexError:
                        pass
                else:
                    try:
                        dnew[ns:ne] = vnew[n]
                    except IndexError:
                        pass
            ##### set_attribute
            new_name = v.replace('_01','_to_20_ku')
            if new_name in self.vars:
                print(new_name+' exists already, ignoring')
                continue
            setattr(self, new_name, dnew)
            self.vars.append(new_name)
            if self.add_attr:
                vattr = add_ncattr(vnew)
                vattr['Merge_note'] = 'Value moved from 1 Hz to 20 Hz using method = '+method
                setattr(self, new_name+'_attr', vattr)
        data_nc.close()

    def print_time(self):
        print(
          CSt2dt(self.time_20_ku[0],string=True)+'--'+
          CSt2dt(self.time_20_ku[-1],string=True))


    ### need a def to add another variable to the time axis
    ### we need a new variable added via string
    ### variable identical length to time/lon/lat
    ### match axis to axis making sure it all fits
    def add_data_time(self,time0,var,varnames,verbos=False,
                     filt_less = None,filt_great = None,tshift = -1e-7,
                     nc_attr = None):
        ### extra reshape
        reshape = False
        if np.shape(time0.shape)[0] > 1:
            time0 = time0.ravel()
            reshape = True
            if verbos: print('Reshaping inputs')
        vuse0 = (time0>=self.time_20_ku[0]-tshift) & (time0<self.time_20_ku[-1]-tshift) 
        if type(varnames) == str:
            var = [var]
            varnames = [varnames]
        for n,v in enumerate(var):
            if varnames[n] in self.vars:
                print(varnames[n]+' exists already, ignoring')
                continue
            ### need empty array
            dnew = np.ma.masked_all([self.n_u])

            ### need to match time to time
            ### remove times under self.time_20_ku[0] and over likewise
            if reshape:
                v = v.ravel()
            if type(v) == np.ma.core.MaskedArray:
                if type(v.mask) == np.ndarray:
                    vuse = vuse0 & np.isfinite(v) & ~v.mask
                else:
                    vuse = vuse0 & np.isfinite(v)
            else:
                vuse = vuse0 & np.isfinite(v)
            if filt_less:
                vuse = vuse & (v>filt_less)
            if filt_great:
                vuse = vuse & (v<filt_great)
            if verbos and n == 0: print('No.data in range: '+str(np.sum(vuse)))
            time = time0[vuse]
            locs = np.searchsorted(self.time_20_ku,time+tshift)
            #### delete repeated entries
            if n==0 and verbos:
                ### lat error
                print("add_data_time av position error")
                print(np.nanmedian(self.time_20_ku[locs]-time-tshift))
            
            dnew[locs] = v[vuse]
            dnew.mask[:] = True
#             if filt_less:
#                 ### remake locs only keeping those where v>filt_less
#                 locs_use = v[vuse]>filt_less
#                 locs_filt = locs[locs_use]
#                 dnew.mask[locs_filt] = False
#             else:
#                 dnew.mask[locs] = False
#             if filt_great:
#                 ### remake locs only keeping those where v<filt_great
#                 locs_use = v[vuse]<filt_great
#                 locs_filt = locs[locs_use]
#                 dnew.mask[locs_filt] = False
#             else:
            dnew.mask[locs] = False
            
            
            if verbos and n == 0: print('No.data found: '+str(locs.shape[0])) 

            ##### set_attribute
            setattr(self,varnames[n], dnew)
            self.vars.append(varnames[n])
            if self.add_attr and nc_attr is not None:
                setattr(self, varnames[n]+'_attr', nc_attr[n])
            elif self.add_attr:
                setattr(self, varnames[n]+'_attr',
                       {'Note':'No nc_attr info available'})

    def add_data_lat(self,lats0,var,varnames,verbos=False,
                     filt_less = None,filt_great = None,
                     nc_attr = None):
        reshape = False
        if np.shape(lats0.shape)[0] > 1:
            lats0 = lats0.ravel()
            reshape = True
        ## need to differentiate between increasign decreasing latitude
        ## increasing
        descend = np.sum(np.diff(self.lat_20_ku))<0.0
        if descend:
            vuse0 = (lats0<self.lat_20_ku[0]) & (lats0>=self.lat_20_ku[-1])
        else:
            vuse0 = (lats0>=self.lat_20_ku[0]) & (lats0<self.lat_20_ku[-1])
        if type(varnames) == str:
            var = [var]
            varnames = [varnames]
        for n,v in enumerate(var):
            if varnames[n] in self.vars:
                print(varnames[n]+' exists already, ignoring')
                continue
            ### need empty array
            dnew = np.ma.masked_all([self.n_u])

            if reshape:
                v = v.ravel()
            ### need to match lat to lat
            ### remove lats under self.lat_20_ku[0] and over likewise
            vuse = vuse0 & np.isfinite(v)
            if verbos and n == 0: print('No.data in range: '+str(np.sum(vuse)))
            lats = lats0[vuse]
            if descend:
                locs = np.searchsorted(-self.lat_20_ku,-lats)
            else:
                locs = np.searchsorted(self.lat_20_ku,lats)
            if n==0 and verbos:
                ### lat error
                print("add_data_lat av position error")
                print(np.nanmedian(self.lat_20_ku[locs]-lats))
            dnew[locs] = v[vuse]
            dnew.mask[:] = True
            if filt_less:
                ### remake locs only keeping those where v>filt_less
                locs_use = v[vuse]>filt_less
                locs = locs[locs_use]
            if filt_great:
                ### remake locs only keeping those where v<filt_great
                locs_use = v[vuse]<filt_great
                locs = locs[locs_use]
            dnew.mask[locs] = False
            if verbos and n == 0: print('No.data found: '+str(locs.shape[0]))

            ##### set_attribute
            setattr(self,varnames[n], dnew)
            self.vars.append(varnames[n])
            if self.add_attr and nc_attr is not None:
                setattr(self, varnames[n]+'_attr', nc_attr[n])
            elif self.add_attr:
                setattr(self, varnames[n]+'_attr',
                       {'Note':'No nc_attr info available'})



    def make_mask(self,vlist,man_list = None,verbos=False):
        ### for example vlist = ['cpom_ssh','height_sea_ice_lead_20_ku','MSSL_height_sea_ice_lead_20_ku']
        self.mask = np.ones_like(self.time_20_ku,dtype=bool)
        for v in vlist:
            vtemp = getattr(self,v) 
            self.mask[vtemp.mask] = 0
        if man_list is not None:
            self.mask[man_list==False] = 0
        if verbos: print('Keeping '+str(np.sum(self.mask))+' data')
    
    def check_attr(self):
        """
        filters all the attr before attempting to load them in to a nc file
        """
        for v in self.vars:
            if hasattr(self,v+'_attr'):
                attr_check = getattr(self,v+'_attr')
                try:
                    del attr_check['_FillValue']
                except KeyError:
                    pass
                try:
                    del attr_check['scale_factor']
                except KeyError:
                    pass
                try:
                    del attr_check['add_offset']
                except KeyError:
                    pass
                    
    

    def save_nc(self,save_dir,orig_root,merge_list = [''],verbos = False,mask=False):
        """
        saves the virtual track as a new track
        """
        ### file name
        ### use L2 file name,
        date_dir = CSt2dt(self.time_20_ku[0]).strftime('%Y/%m/')
        ### check if new dir exists, if not, make it 
        save_long_dir = os.path.dirname(save_dir+self.name.split(orig_root)[1].split(date_dir)[0]+date_dir)
        # check it exists
        if not os.path.exists(save_long_dir):
            # make if it doesn't
            os.makedirs(save_long_dir)
            print('Creating directory: ',save_long_dir)
        ### use exact same L2 names but in the save_dir
        file = save_dir+self.name.split(orig_root)[1]
        if verbos: print(file)
        
        ### make a netcdf file
        if verbos: print('Saving NetCDF: '+file)
        f_nc = Dataset(file, mode='w')

        # Create the dimensions of the file
        # this is time_20_ku
        if mask and hasattr(self,'mask'):
            f_length = self.mask.sum()
            write_mask = self.mask.data
        else: 
            f_length = self.time_20_ku.shape[0]
            write_mask = np.ones([f_length],dtype=bool)
        f_nc.createDimension('time_20_k', f_length)

        # Copy the global attributes
        Fattr = {
            'Description':'CSAO merged data',
            'CSAO_code_by':'HDBS Heorton',
            'Data_Partners':' '.join(merge_list),
            'Original_L2_file':'Same as this file name'
        }
        f_nc.setncatts(Fattr)

        # Create the variables in the file
        ## time
        v = 'time_20_ku'
        fill_var = getattr(self,v)
        f_nc.createVariable(v, fill_var.dtype, ('time_20_k'))

        ### fill it
        fill_var = getattr(self,v)
        f_nc.variables[v][:] = fill_var[write_mask]
        ## lon/lat
#         for v in self.vars[:2]:
#             if verbos: print(v)
#             fill_var = getattr(self,v)
#             f_nc.createVariable(v, fill_var.dtype, ('time_20_k'))
            
#             # Copy the variables values (as 'f4' eventually)
#             fill_var = getattr(self,v)
#             f_nc.variables[v][:] = fill_var[write_mask]
        ## all the others   
        for v in self.vars:
            if verbos: print(v)
            fill_var = getattr(self,v)
            f_nc.createVariable(v, fill_var.dtype, ('time_20_k'))

            # Copy the variable attributes
            fill_var_attr = getattr(self,v+'_attr')
            f_nc.variables[v].setncatts(fill_var_attr)

            # Fill it
            fill_var = getattr(self,v)
            f_nc.variables[v][:] = fill_var[write_mask]

        # Save the file
        f_nc.close()

        
        
### function that takes a CStrack, then depending on the mask
### and variables to copy, makes an abbrevated version
### keeping lon lat time, but only where it's mask is true
### actually we'll have an generalised object first
class CS_var_list():
    def __init__(self):
        self.vlist=[]
        self.vnumb = 0
    def append(self,varnew,vmask):
        self.vlist.append(varnew[vmask])
        self.vnumb += 1
    def ravel(self,ruse=[0,False]):
        ruse1 = [ruse[0],ruse[1]]
        if ruse[1] == False: ruse1[1] = self.vnumb
        return np.array([v for vl in self.vlist[ruse1[0]:ruse1[1]] for v in vl])
        

class CS_track_list():
    def __init__(self,name,varkeep):
        self.name = name
        for v in varkeep:
            setattr(self,v,CS_var_list())
        self.varadd = varkeep
        setattr(self,'time_20_ku',CS_var_list())
        setattr(self,'lon_20_ku',CS_var_list())
        setattr(self,'lat_20_ku',CS_var_list())
    def addtrack(self,CS2_track,varstokeep):
        if CS2_track.mask.sum()>0:
            self.time_20_ku.append(CS2_track.time_20_ku,CS2_track.mask)
            self.lon_20_ku.append(CS2_track.lon_20_ku,CS2_track.mask)
            self.lat_20_ku.append(CS2_track.lat_20_ku,CS2_track.mask)
            for n,v in enumerate(self.varadd):
                if n==0: print('Adding data to '+self.name)
                x = getattr(self,v)
                x.append(varstokeep[n],CS2_track.mask)  
    def shape(self):
        shape = 0
        for v in self.time_20_ku.vlist:
            shape += np.shape(v)[0]
        return shape
    
    def check_attr(self):
        """
        filters all the attr before attempting to load them in to a nc file
        """
        for v in ['lon_20_ku','lat_20_ku']+self.varadd:
            variable = getattr(self,v)
            if hasattr(variable,'attr'):
                attr_check = getattr(variable,'attr')
                try:
                    del attr_check['_FillValue']
                except KeyError:
                    pass
                try:
                    del attr_check['scale_factor']
                except KeyError:
                    pass
                try:
                    del attr_check['add_offset']
                except KeyError:
                    pass
    
    def save_nc(self,save_dir,info_name,merge_list = [''],add_attr = False,verbos=False):
        """
        saves the virtual track list as a new track
        adds attributes too if they are on each variable.
        before writing, add the attribute with:
            CS_track_list.variable.attr = {
                    "label":"info in text form",
                    "more_labels":"more info in text form",
                    }
        """
        ### file name
        ### get time vector end points
        t0str = CSt2dt(self.time_20_ku.vlist[0][0]).strftime('%Y-%m-%dT%H-%M_')
        t1str = CSt2dt(self.time_20_ku.vlist[-1][-1]).strftime('%Y-%m-%dT%H-%M')
        ### check if new dir exists, if not, make it 
        save_long_dir = os.path.dirname(save_dir+
                            CSt2dt(self.time_20_ku.vlist[0][0]).strftime('%Y'))
        # check it exists
        if not os.path.exists(save_long_dir):
            # make if it doesn't
            os.makedirs(save_long_dir)
            print('Creating directory: ',save_long_dir)
        ### use exact same L2 names but in the save_dir
        file = save_long_dir+'/'+info_name+t0str+t1str+'.nc'
        
        ### make a netcdf file
        if verbos: print('Saving NetCDF: '+file)
        f_nc = Dataset(file, mode='w')

        # Create the dimensions of the file
        # this is time_20_ku
        f_length = self.shape()
        f_nc.createDimension('time_20_k', f_length)

        # Copy the global attributes
        Fattr = {
            'Description':'CSAO merged data list',
            'CSAO_code_by':'HDBS Heorton',
            'Data_Partners':' '.join(merge_list),
            'No_of combined_tracks':str(len(self.time_20_ku.vlist)),
        }
        f_nc.setncatts(Fattr)

        # Create the variables in the file
        ## time
        v = 'time_20_ku'
        fill_var = getattr(self,v)
        # Fill it
        f_nc.createVariable(v, fill_var.vlist[0][0].dtype, ('time_20_k'))

        ### fill it
        f_nc.variables[v][:] = fill_var.ravel()
        ## lon/lat and all
        for v in ['lon_20_ku','lat_20_ku']+self.varadd:
            fill_var = getattr(self,v)
            f_nc.createVariable(v, fill_var.vlist[0][0].dtype, ('time_20_k'))
            if add_attr and hasattr(fill_var,'attr'):
                fill_var_attr = getattr(fill_var,'attr')
                f_nc.variables[v].setncatts(fill_var_attr)
            # Fill it
            f_nc.variables[v][:] = fill_var.ravel()
            
        # Save the file
        f_nc.close()


class NSIDC_nt():
    """
    forcing class for the budget
    lets the forcing load efficiently
    
    """
    def __init__(self,ppath):
        self.name = 'NSIDC_n'
        self.path = ppath
# next function will take a list of dates and return an appropriately orientated arrays
# give a 
    def get_aice(self,dates_u,verbos=False):
        # does dates_u cover one year or more
        #daily files
        dimY = 316
        dimX = 332
        d_no = np.shape(dates_u)[0]
        data =  np.empty([d_no, dimX, dimY])
        for n,d in enumerate(dates_u):
#             infile = self.path+d.strftime('/%Y/')+"nt_"+d.strftime('%Y%m%d')+"_f17_v1.1_n.bin"
            if d.year<2020:
#                 infile = self.path+"/nt_"+d.strftime('%Y%m%d')+"_f17_v1.1_n.bin"
                infile = self.path+d.strftime('/%Y/')+"nt_"+d.strftime('%Y%m%d')+"_f17_v1.1_s.bin"
            if d.year>2019:
#                 infile = self.path+"/nt_"+d.strftime('%Y%m%d')+"_f18_nrt_n.bin"
                infile = self.path+d.strftime('/%Y/')+"nt_"+d.strftime('%Y%m%d')+"_f18_nrt_s.bin"
#             if d.year<2019:
#                 infile = self.path+"nt_"+d.strftime('%Y%m%d')+"_f17_v1.1_n.bin"
#             if d.year>2018:
#                 infile = self.path+"nt_"+d.strftime('%Y%m%d')+"_f18_nrt_n.bin"
            with open(infile, 'rb') as fr:
                hdr = fr.read(300)
                ice = np.fromfile(fr, dtype=np.uint8)

            ice = ice.reshape(dimX,dimY)
            ice = np.flipud(ice)
            data[n] = ice / 250.
        data[data>1.0] = np.nan
        return data

    def get_dates(self,time_start,time_end):
        # does dates_u cover one year or more
        #daily files
        dates_u = []
        d_no = (time_end-time_start).days +3 
        # make sure we get the bracket points
        for dn in range(d_no):
            d = time_start+ relativedelta(days = dn - 1)
            if d.year<2020:
                infile = self.path+d.strftime('/%Y/')+"nt_"+d.strftime('%Y%m%d')+"_f17_v1.1_n.bin"
            if d.year>2019:
                infile = self.path+d.strftime('/%Y/')+"nt_"+d.strftime('%Y%m%d')+"_f18_nrt_n.bin"
            # check infile exists 
            if exists(infile):
                dates_u.append(d)
            #if it does append dates_u
        self.dates= dates_u
        print(self.name+' Found '+str(np.shape(dates_u)[0])+' dates')


class add_file:
    def __init__(self,name,path,dir_form,multi_file = True,name_form = '%Y%m%d',name_blurb = ''):
        #### set the format for how we access files
        #### paths, extension, blurbs, etc
        self.name = name
        self.path = path
        self.dir_form = dir_form
        self.name_form = name_form
        self.name_blurb = name_blurb
        self.multi_file = multi_file
        self.nc_loaded = False
        self.tshift = -1e-7
        self.F_no = 0
        self.F_pt = slice(0,None)
        self.F_pt_end = 0
        
    def get_files_times(self,time,verbos=False):
        #### uses the attr to load all the files we care about
        #### resolution of tw (default to a day)
        if self.multi_file:
            self.files = []
            self.times = []
            tstring = time.strftime(self.name_form)
            ystring = time.strftime('%Y')
            temp_files = glob.glob(self.path+time.strftime(self.dir_form)+
                        self.name_blurb+'*.nc')
#                         '*.nc')
            temp_files.sort()
            for file in temp_files:
                if tstring in file:
                    try:
                        ftime = dt.datetime.strptime(
                            tstring+file.split(self.path)[1].split(tstring)[1],
                            '%Y%m%dT%H%M%S_')
                        self.times.append(dt2CSt(ftime))
                        self.files.append(file)
                    except ValueError:
                        try:
                            ### second time is next day
                            ftime = dt.datetime.strptime(tstring+
                                file.split(self.path)[1].split(tstring)[1].split(ystring)[0],
                                '%Y%m%dT%H%M%S_')
                            self.times.append(dt2CSt(ftime))
                            self.files.append(file)
                        except ValueError:
                            pass
            self.times = np.array(self.times)
            if verbos: print(self.name+' found '+str(self.times.shape[0])+' files')
        elif not self.nc_loaded:
            ### there 'should; be a new file the day in it
            self.files = glob.glob(self.path+
                time.strftime(self.dir_form)+'*'+time.strftime(self.name_form)+'*.nc')
            ### open if there
            try:
                self.f_nc = Dataset(self.files[self.F_no])
                Ftime = self.f_nc['time_20_ku']
                ### leave f_nc open
                gaps = np.diff(Ftime)
                self.track_points = np.hstack([[0],np.where(gaps>1)[0]+1,[Ftime.size-1]] )
                self.times = np.array([Ftime[t] for t in self.track_points])
                self.Fmax = self.times.shape[0]
                self.nc_loaded = True
                if verbos: print(self.name+' found '+str(self.times.shape[0])+' file slices')
            except IndexError:
                ### leave it unloaded
                self.times = np.array([])
                self.Fmax = 0
        ### end up with self.times to search
        ### plus the way to find the file list or merged
                
                

    def get_tp(self,ftime,tolerance,verbos = False):
        #### set whether we do each of these
        tp_ar = np.where(np.abs(self.times-ftime)<tolerance)[0]
        if tp_ar.shape[0]>0:
            ### files exists
            #### looks for the time/files corresponding to ftime
            #### matches file time start points
            if self.multi_file:
                self.F_no = tp_ar[0] 
                if verbos: print(self.name+' opens file no:'+str(self.F_no))
                return True
            #### find array point to access from
            #### work for merged files
            else:
                if tp_ar[0] == self.Fmax-1:
                    ### end of the file so don't
                    self.f_nc.close()
                    self.nc_loaded = False
                    self.times = np.array([])
                    self.Fmax = 0
                    return False
                else:
                    self.F_pt = slice(self.track_points[tp_ar[0]],
                                      self.track_points[tp_ar[0]+1])
                    self.F_pt_end = tp_ar[0]+1
                    if verbos: print(self.name+' slice no:'+str(self.F_pt_end-1)+
                                     ' of '+str(self.Fmax))
                    return True
        else:
            return False

    def load_nc(self):
        #### checks whether we need to load a new nc file or not
        if self.multi_file:
            self.f_nc = Dataset(self.files[self.F_no])

    def load_clean(self):
        #### tidies the nc_file if we need to
        if self.multi_file:
            self.f_nc.close()
            #### end of file
        elif self.F_pt_end == self.Fmax:
            self.f_nc.close()
            self.nc_loaded = False
            self.times = np.array([])
            self.Fmax = 0
            

def add_atm_cor(track):
    cor = np.ma.zeros([track.n_u])
    cor_list = [
 'mod_wet_tropo_cor_to_20_ku',
 'mod_dry_tropo_cor_to_20_ku',
 'iono_cor_gim_to_20_ku',
 'inv_bar_cor_to_20_ku',
 'hf_fluct_total_cor_to_20_ku']
    for cl in cor_list:
        cor += getattr(track,cl)
    return cor

def add_geo_cor(track):
    cor = np.ma.zeros([track.n_u])
    cor_list = [
 'pole_tide_to_20_ku',
 'solid_earth_tide_to_20_ku',
 'load_tide_to_20_ku',
 'ocean_tide_to_20_ku',
 'ocean_tide_eq_to_20_ku',
 ]
    for cl in cor_list:
        cor += getattr(track,cl)
    return cor
    