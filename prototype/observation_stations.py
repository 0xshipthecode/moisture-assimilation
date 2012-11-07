


from time_series_utilities import match_sample_times
from spatial_model_utilities import find_closest_grid_point, great_circle_distance

import re
import pytz
import codecs
from datetime import datetime, timedelta


class Observation:
    """
    An observation of a field value at a certain time.
    """
    
    def __init__(self, s, tm, fm_obs):
        """
        Constructs an observation packet from the info dictionary.
        """
        self.s = s
        self.tm = tm
        self.obs_val = fm_obs
        
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
    An observation station which yields observations.
    """
    
    def __init__(self, file_name, time_zone, wrf_data):
        """
        Load a station from data file.
        """
        self.tz = time_zone
        self.file = file_name

        # array of values
        self.tm = []
        self.fuel_moisture = []
        self.rh = []
        self.rain = []
        self.air_temp = []
        self.fuel_temp = []
        
        self.load_station_data(file_name, time_zone)
        
        mlon, mlat = wrf_data.get_lons(), wrf_data.get_lats()
        self.grid_pt = find_closest_grid_point(self.lon, self.lat, mlon, mlat)
        self.dist_grid_pt =  great_circle_distance(self.lon, self.lat, mlon[self.grid_pt], mlat[self.grid_pt])
        

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
    
    
    def get_fuel_moisture(self):
        """
        Return the fuel moisture time series.
        """
        return self.fuel_moisture
    
    
    def get_elevation(self):
        """
        Return the elevation in meters above sea level.
        """
        return self.elevation

        
    def load_station_data(self, station_file, tz):
        """
        Load all available fuel moisture data from the station information.
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
                self.fuel_moisture.append(float(fields[7]) / 100.0)
                self.rh.append(float(fields[8]))
                self.rain.append(float(fields[11]))
                l = f.readline()
            
            while l[:5] != 'Daily' and len(l) > 0:
                l = f.readline()
                
            if len(l) == 0:
                break
                    
        f.close()
        
        
    def get_observations_for_times(self, obs_type, tm):
        """
        Get the observations of the field which match the times passed in tm.
        """
        _, indx_me, indx_other = match_sample_times(self.tm, tm)
        ts = vars(self)[obs_type]
        return indx_other, [Observation(self, self.tm[i], ts[i]) for i in indx_me]
