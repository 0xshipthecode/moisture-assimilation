


from time_series_utilities import match_sample_times
from spatial_model_utilities import find_closest_grid_point, great_circle_distance

import re
import pytz
import codecs
from datetime import datetime, timedelta
import xlrd


class Observation:
    """
    An observation of a field value at a certain time.
    """
    
    def __init__(self, s, tm, fm_obs, fm_var, field_name):
        """
        Constructs an observation packet from the info dictionary.
        """
        self.s = s
        self.tm = tm
        self.obs_val = fm_obs
        self.field_name = field_name
        self.obs_var = fm_var

        
    def get_time(self):
        """
        Return time of observation.
        """
        return self.tm
        
        
    def get_value(self):
        """
        Return the observation value.
        """
        return self.obs_val
    

    def get_measurement_variance(self):
        """
        Return the variance of the measurement of the given field.
        """
        return self.obs_var
        

    def get_position(self):
        """
        Longitude and lattitude of the originating station (shortcut).
        """
        return self.s.get_position()
    

    def get_nearest_grid_point(self):
        """
        Return the indices that identify the nearest grid point.
        """
        return self.s.get_nearest_grid_point()
    

    def get_station(self):
        """
        Return the station from which this observation originates 
        """
        return self.s



class Station:
    """
    An observation station which stores and yields observations.
    All times must be in GMT.
    """
    def __init__(self, wrf_data=None):
        """
        Load a station from data file.
        """
        # array of observation times and dictionary of obs variables
        self.tm = []
        self.obs_vars = {}
        
        # only co-register to grid if required
        if wrf_data is not None:
            mlon, mlat = wrf_data.get_lons(), wrf_data.get_lats()
            self.grid_pt = find_closest_grid_point(self.lon, self.lat, mlon, mlat)
            self.dist_grid_pt =  great_circle_distance(self.lon, self.lat,
                                                       mlon[self.grid_pt], mlat[self.grid_pt])
        else:
            self.grid_pt = None
            self.dist_grid_pt = None
    
        # the measurement variance of different observations is unknown
        self.measurement_variance = {}
        

    def get_name(self):
        """
        Returns the name of the station.
        """
        return self.name
    
    
    def get_position(self):
        """
        Get geographical position of the observation station as a (lat, lon) tuple.
        """
        return self.lon, self.lat

    
    def get_nearest_grid_point(self):
        """
        Returns the indices of the nearest grid point.
        """
        return self.grid_pt
    
    
    def get_dist_to_grid(self):
        """
        Returns the distance in kilometers to the nearest grid point.
        """
        return self.dist_grid_pt
    
        
    def get_obs_times(self):
        """
        Return observatuion times.
        """
        return self.tm
    
    
    def get_elevation(self):
        """
        Return the elevation in meters above sea level.
        """
        return self.elevation
    
    
    def get_measurement_variance(self, field_name):
        """
        Return the variance of the measurement from this station.
        """
        return self.measurement_variance[field_name]
    
    
    def set_measurement_variance(self, field_name, v):
        """
        Set the variance of the measurements from this station.
        """
        self.measurement_variance[field_name] = v


    def get_observations_raw(self, obs_name):
        """
        Return an entire time series of given observation type.  Returns
        the raw observation values.  Call must manage time stamps.
        """
        return self.obs_vars[obs_name]


    def get_observations(self, obs_name):
        """
        Returns a list of Observations for given observation type (var name).
        """
        obs = self.get_observations_raw(obs_name)
        mv = self.get_measurement_variance(obs_name)
        return [Observation(self, self.tm[i], obs[i], mv, obs_type) for i in range(len(obs))]

        
    def get_observations_for_times(self, obs_type, tm):
        """
        Get the observations of the field which match the times passed in tm.
        Also returns an index array which shows which times in the argument tm
        match the observation times returned.
        Returns a set of Observations.
        """
        _, indx_me, _ = match_sample_times(self.tm, tm)
        if len(indx_me) < len(tm):
            raise ValueError('Observations for some times not available.')
        ts = self.get_observations_raw(obs_type) 
        mv = self.get_measurement_variance(obs_type)
        return [Observation(self, self.tm[i], ts[i], mv, obs_type) for i in indx_me]



class StationAdam(Station):
    """
    An specialization of observation station which loads data from
    the format Adam sent to me for the Witch Creek data.
    """
    
    def __init__(self, wrf_data):
        """
        Load a station from data file.
        """
        Station.__init__(self, wrf_data)
        

    def load_station_data(self, station_file, tz):
        """
        Load all available fuel moisture data from the station information text file.
        These files have been supplied by A. Kochanski.
        """
        f = codecs.open(station_file, 'r', encoding = 'utf-8')
        s = f.readline().strip()
        self.name = str(s[0:min(8, len(s))])
        
        # next line is location string, not interesting    
        f.readline()
        
        # next 2 lines are lattitude & longitude
        lat_lon_re = re.compile("\D+(\d{2,3})\D+(\d{2})\D+(\d{2})")
    
        l = f.readline()
        mo = lat_lon_re.match(l)
        lat_info = map(lambda x: int(x), mo.groups())
        self.lat = lat_info[0] + lat_info[1] / 60.0 + lat_info[2] / 3600.0
    
        l = f.readline()
        mo = lat_lon_re.match(l)
        lon_info = map(lambda x: int(x), mo.groups())
        self.lon = -(lon_info[0] + lon_info[1] / 60.0 + lon_info[2] / 3600.0)
        
        # read off the elevation in feet and convert to meters
        l = f.readline()
        self.elevation = int(l[9:l.find(' ft')]) * 0.3048
        
        # bypass lines 6 through 8
        for i in range(6,8):
            f.readline()
    
        while True:
    
            # read in and parse date
            l = f.readline()
            date = datetime.strptime(l.strip(), '%B %d, %Y').replace(tzinfo = tz)
            
            # read lines until a line starts with daily
            while l[0] < '0' or l[0] > '9' and len(l) > 0:
                l = f.readline()
                
            if len(l) == 0:
                break
            
            while l[0] >= '0' and l[0] <= '9' and len(l) > 0:
                fields = filter(lambda x: len(x) > 0, l.split('\t'))
                time = datetime.strptime(fields[0], "%I %p")
                timed = timedelta(0, time.hour * 3600)
                mtime = date + timed
                if len(fields) != 12:
                    print fields
                    
                self.tm.append(mtime)
                self.air_temp.append(float(fields[5]))
                self.fuel_temp.append(float(fields[6]))
                self.fm10.append(float(fields[7]) / 100.0)
                self.rh.append(float(fields[8]))
                self.rain.append(float(fields[11]))
                l = f.readline()
            
            while l[:5] != 'Daily' and len(l) > 0:
                l = f.readline()
                
            if len(l) == 0:
                break
                    
        f.close()
        
        
class MesoWestStation(Station):
    """
    An observation station with data downloaded from the MesoWest website in xls format.
    """
    
    def __init__(self, info_string, wrf_data = None):
        """
        Initialize the station using an info_string that is written by the scrape_stations.py
        script into the 'station_infos' file.
        """
        # initialize a flag that will be set to False if any data is missing or invalid
        # when loading
        self.data_loaded_ok = True

        # parse the info_string
        tokens = info_string.split(',')
        self.name = tokens[0]
        self.lat = float(tokens[1])
        self.lon = float(tokens[2])
        self.elevation = float(tokens[3]) * 0.3048

        Station.__init__(self, wrf_data)
        
        # construct empty variable lists for what we wish to load
        self.obs_vars['air_temp'] = []
        self.obs_vars['rh'] = []
        self.obs_vars['fm10'] = []


    def load_station_data(self, station_file):
        """
        Load all available fuel moisture data from the station measurement file
        in xls format.
        """
        # load the worksheet
        x = xlrd.open_workbook(station_file)
        s = x.sheet_by_index(0)

        # find the order of the variables
        var_ord = []
        cell_ord = []
        for i in range(1,5):
            cv = s.cell_value(0,i)
            if 'TMP' in cv:
                var_ord.append(self.obs_vars['air_temp'])
                cell_ord.append(i)
            elif 'RELH' in cv:
                var_ord.append(self.obs_vars['rh'])
                cell_ord.append(i)
            elif 'FM' in cv:
                var_ord.append(self.obs_vars['fm10'])
                cell_ord.append(i)

        if len(var_ord) != 3:
            self.data_loaded_ok = False

        # now read 24 entries starting at 
        i = 1
        gmt_tz = pytz.timezone('GMT')
        while True:
            # parse the time stamp string
            try:
                tstamp = gmt_tz.localize(datetime.strptime(s.cell_value(i,0), '%m-%d-%Y %H:%M %Z'))
                self.tm.append(tstamp.replace(minute=0)) 
            except ValueError:
                break

            # parse the variables in order
            for j in range(len(cell_ord)):
                try:
                    var_ord[j].append(float(s.cell_value(i,cell_ord[j])))
                except ValueError:
                    var_ord[j].append(float('nan'))
                    self.data_loaded_ok = False

            i += 1

        if i != 25:
            self.data_loaded_ok = False


    def data_ok(self):
        """
        Check if data loaded without any errors.
        """
        return self.data_loaded_ok


if __name__ == '__main__':

    o = MesoWestStation('../real_data/colorado_stations/BAWC2.xls', 'BAWC2,39.3794,-105.3383,2432.9136') 
    print(o.get_observations('fm10'))

