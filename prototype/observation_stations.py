

from time_series_utilities import match_sample_times
from spatial_model_utilities import find_closest_grid_point, great_circle_distance

import re
import pytz
import codecs
from datetime import datetime, timedelta
import string


def readline_skip_comments(f):
    """
    Read a new line while skipping comments.
    """
    l = f.readline().strip()
    while len(l) > 0 and l[0] == '#':
        l = f.readline().strip()
    return l


class Observation:
    """
    An observation of a field value at a certain time.
    """
    
    def __init__(self, s, tm, obs, var, field_name):
        """
        Constructs an observation packet from the info dictionary.
        """
        self.s = s
        self.tm = tm
        self.obs_val = obs
        self.field_name = field_name
        self.obs_var = var

        
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
    def __init__(self):
        """
        Load a station from data file.
        """
        # array of observation times and dictionary of obs variables
        self.tm = []
        self.obs_vars = {}

        # no co-registration by default
        self.grid_pt = None
        self.dist_grid_pt = None



    def register_to_grid(self, wrf_data):
        """
        Find the nearest grid point to the current location.
        """
        # only co-register to grid if required
        mlon, mlat = wrf_data.get_lons(), wrf_data.get_lats()
        self.grid_pt = find_closest_grid_point(self.lon, self.lat, mlon, mlat)
        self.dist_grid_pt =  great_circle_distance(self.lon, self.lat,
                                                   mlon[self.grid_pt], mlat[self.grid_pt])


    def get_id(self):
        """
        Returns the id of the station.
        """
        return self.id


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
    
    

    def get_observations(self, obs_type):
        """
        Returns a list of Observations for given observation type (var name).
        """
        return self.obs[obs_type]

        

class StationAdam(Station):
    """
    An specialization of observation station which loads data from
    the format Adam sent to me for the Witch Creek data.
    """
    
    def __init__(self):
        """
        Load a station from data file.
        """
        Station.__init__(self)
        
        # construct empty variable lists for what we wish to load
        for var in [ 'air_temp', 'rh', 'fm10', 'fuel_temp', 'rain']:
            self.obs_vars[var] = []


    def load_station_data(self, station_file, tz):
        """
        Load all available fuel moisture data from the station information text file.
        These files have been supplied by A. Kochanski.
        """
        f = codecs.open(station_file, 'r', encoding = 'utf-8')
        s = f.readline().strip()
        self.name = str(s[0:min(8, len(s))])
        self.id = self.name
        
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
                self.obs_vars['air_temp'].append(float(fields[5]))
                self.obs_vars['fuel_temp'].append(float(fields[6]))
                self.obs_vars['fm10'].append(float(fields[7]) / 100.0)
                self.obs_vars['rh'].append(float(fields[8]))
                self.obs_vars['rain'].append(float(fields[11]))
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
    
    def __init__(self, name):
        """
        Initialize the station using an info_string that is written by the scrape_stations.py
        script into the 'station_infos' file.
        """
        # initialize a flag that will be set to False if any data is missing or invalid
        # when loading
        self.data_loaded_ok = True

        # parse the info_string
        self.name = name

        Station.__init__(self)


    def load_station_info(self, station_info):
        """
        Load station information from an .info file.                                                                  """
        with open(station_info, "r") as f:

            # read station id
            self.id = readline_skip_comments(f)

            # read station name
            self.name = readline_skip_comments(f)

            # read station geo location
            loc_str = readline_skip_comments(f).split(",")
            self.lat, self.lon = float(loc_str[0]), float(loc_str[1])

            # read elevation
            self.elevations = float(readline_skip_comments(f))

            # read sensor types
            self.sensors = map(lambda x: x.strip(), readline_skip_comments(f).split(","))

            # create empty lists for observations
            self.obs = {}
            for s in self.sensors:
                self.obs[s] = []



    def load_station_data(self, station_file):
        """
        Load all available fuel moisture data from the station measurement file
        in an obs file.
        """
        gmt_tz = pytz.timezone('GMT')

        with open(station_file, "r") as f:

            while True:

                # read in the date or exit if another packet is not found
                tm_str = readline_skip_comments(f)
                if len(tm_str) == 0:
                    break

                tstamp = gmt_tz.localize(datetime.strptime(tm_str, '%Y-%m-%d_%H:%M %Z'))
            
                # read in the variable names
                var_str = map(string.strip, readline_skip_comments(f).split(","))
            
                # read in observations
                vals = map(lambda x: float(x), readline_skip_comments(f).split(","))

                # read in variances
                variances = map(lambda x: float(x), readline_skip_comments(f).split(","))

                # construct observations
                for vn,val,var in zip(var_str, vals, variances):
                    self.obs[vn].append(Observation(self, tstamp, val, var, vn))


    def data_ok(self):
        """
        Check if data loaded without any errors.
        """
        return self.data_loaded_ok


if __name__ == '__main__':

    o = MesoWestStation('../real_data/colorado_stations/BAWC2.xls', 'BAWC2,39.3794,-105.3383,2432.9136') 
    print(o.get_observations('fm10'))

