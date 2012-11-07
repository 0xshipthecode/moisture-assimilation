

import netCDF4
import pytz
from datetime import datetime, timedelta
import os
import numpy as np



class WRFModelData:
    """
    This class contains aggregate information loaded from a WRF model, methods for loading data from a WRF simulation
    are provided.
    """
    
    def __init__(self, file_name, fields = None, tz_name = None):
        """
        Load data from a file file_name. See load_wrf_data for standard fields that
        are loaded.  The fields can be overridden by passing a new list in the fields
        argument.  The model simulation times can be moved into a different time zone
        by passing in a time zone descriptor in tz_name (must be recognizable for pytz).
        If no time zone is given, the get_times() function assumes GMT is local time. 
        """
        self.file_name = file_name
        self.load_data(file_name, fields)
        if tz_name:
            self.change_time_zone(tz_name)
    

    def load_data(self, data_file,
                  var_names):
        """
        Load required variables from the file data_file.  A list of variables
        is either supplied or the default list is used which contains the following
        variables: 'T2', 'Q2', 'PSFC', 'XLAT', 'XLONG', 'RAINNC'.  The fields
        'Times', 'XLAT', 'XLONG' are always loaded.
        """
        
        # replace empty array by default
        if var_names is None:
            var_names = ['T2', 'Q2', 'PSFC', 'RAINNC']
            
        self.fields = {}
        
        d = netCDF4.Dataset(os.path.join(data_file))
        for vname in var_names:
            self.fields[vname] = d.variables[vname][:,...]
            
        # we assume the domain does not move in time
        self.fields['lat'] = d.variables['XLAT'][0,:,:]
        self.fields['lon'] = d.variables['XLONG'][0,:,:]
            
        # time is always loaded and encoded as a list of python datetime objects
        gmt_tz = pytz.timezone('GMT')
        tm = d.variables['Times'][:,...]
        tp = []
        for t in tm:
            dt = datetime.strptime(''.join(t), '%Y-%m-%d_%H:%M:%S')
            dt = dt.replace(tzinfo = gmt_tz)
            tp.append(dt)
            
        self.fields['GMT'] = tp
        self.fields['LT'] = tp
 
        d.close()
        
        # precompute the equilibrium fields needed everywhere
        self.equilibrium_moisture()

    
    def change_time_zone(self, tz_name):
        """
        Changes the local_time variable to a new timezone.  The local time
        can be accessed using the get_times() function.
        """
        # recode times into datetime objects
        gt = self['GMT']
        self.tz = pytz.timezone(tz_name)
        lt = [ dt.astimezone(self.tz) for dt in gt]
        self.fields['LT'] = lt
        
        
    def get_times(self):
        """
        Returns the local time (depends on time zone set).
        """
        return self['LT']
    
    
    def get_lons(self):
        """
        Return longitude of grid points.
        """
        return self['lon']
    
    
    def get_lats(self):
        """
        Return lattitute of grid points.
        """
        return self['lat']
    
    def get_field(self, field_name):
        """
        Return the field with the name field_name.
        """
        return self.field[field_name]
    
    
    def get_domain_extent(self):
        """
        Return smallest enclosing aligned rectangle of domain.
        return is a tuple (min(lon), min(lat), max(lon), max(lat)). 
        """
        lat = self['lat']
        lon = self['lon']
        return (np.min(lon), np.min(lat), np.max(lon), np.max(lat))
    

    def __getitem__(self, name):
        """
        Access a variable from the fields dictionary.
        """
        return self.fields[name]


    def equilibrium_moisture(self):
        """
        Uses the fields of the WRF model to compute the equilibrium
        field.
        """
        
        # load the standard fields
        P = self['PSFC']
        Q = self['Q2']
        T = self['T2'] 
    
        # saturated vapor pressure (at each location, size n x 1)
        Pws = np.exp(54.842763 - 6763.22/T - 4.210 * np.log(T) + 0.000367*T + np.tanh(0.0415*(T - 218.8)) 
            * (53.878 - 1331.22/T - 9.44523 * np.log(T) + 0.014025*T))
          
        # water vapor pressure (at each location, size n x 1)
        Pw = P * Q / (0.622 + (1 - 0.622) * Q)
        
        # relative humidity (percent, at each location, size n x 1)
        H = 100 * Pw / Pws;
        
        # drying/wetting fuel equilibrium moisture contents (location specific,
        # n x 1)
        Ed = 0.924*H**0.679 + 0.000499*np.exp(0.1*H) + 0.18*(21.1 + 273.15 - T)*(1 - np.exp(-0.115*H))
        Ew = 0.618*H**0.753 + 0.000454*np.exp(0.1*H) + 0.18*(21.1 + 273.15 - T)*(1 - np.exp(-0.115*H))
    
        Ed *= 0.01
        Ew *= 0.01
        
        self.fields['Ed'] = Ed
        self.fields['Ew'] = Ew


    def get_moisture_equilibria(self):
        """
        Return the drying and wetting equilibrium.
        """
        return self['Ed'], self['Ew']

