#! /usr/bin/env python

import xlrd
import sys
import string
import glob
import pytz
from datetime import datetime


skip_vars = [ 'QFLG' ]
process_vars = { 'FM' : lambda x: 0.01 * x,  # fuel moisture comes in gm water /100 gm pine wood
                 'TMPF' : lambda x: 273.15 + x,  # we actually download TMPF in degrees C in the scraper
                 'RELH' : lambda x : 0.01 * x,  # relative humidity comes in percent
                 'SKNT' : lambda x : 273.15 + x # skin temperature is also in deg C
               }


def load_station_data(station_file):
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
    obs = {}
    for i in range(1,s.ncols):
        cv = s.cell_value(0,i).strip().split(" ")[0]
        if len(cv) > 0 and cv not in skip_vars:
            var_ord.append(cv)
            cell_ord.append(i)
            if cv not in process_vars:
                process_vars[cv] = lambda x: x

    # now read 24 entries starting at 
    i = 1
    gmt_tz = pytz.timezone('GMT')
    while i < s.nrows:
        # parse the time stamp string
        try:
            tstamp = gmt_tz.localize(datetime.strptime(s.cell_value(i,0), '%m-%d-%Y %H:%M %Z'))
        except ValueError:
            break

        # parse the variables in order
        obs_i = []
        for j in range(len(cell_ord)):
            try:
                val = process_vars[var_ord[j]](float(s.cell_value(i,cell_ord[j])))
                obs_i.append((var_ord[j], val))
            except ValueError:
                pass

        obs[tstamp] = obs_i

        i += 1

    return obs


if __name__ == "__main__":

    gmt_tz = pytz.timezone("GMT")

    if len(sys.argv) < 3:
        print("Usage: extract_observations.py station_list obs_var_table")
        sys.exit(1)

    # read in station codes
    with open(sys.argv[1], "r") as f:
        sids = filter(lambda x: len(x) > 0 and x[0] != '#', map(string.strip, f.readlines()))

    # read in standard observation variances
    with open(sys.argv[2], "r") as f:
        obs_var_tbl = dict(map(lambda x: x.split(","), map(string.strip, f.readlines())))
        for k,v in obs_var_tbl.iteritems():
            obs_var_tbl[k] = float(v)

    # for each station find all xls files
    for sid in sids:
        print("Processing station %s" % sid)
        fnames = glob.glob("%s*.xls" % sid)

        stored_obs = {}

        # read in each file and parse the contents
        obs = {}
        for fname in fnames:
            obs_i = load_station_data(fname)
            obs.update(obs_i)

        # write out the data in a neat csv file
        with open("%s.obs" % sid, "w") as f:
            f.write("# Data file generated on %s by extract_observations.py\n" % str(datetime.now()))
            f.write("# Format is time, observed_vars, observations, variances\n")
            for tm in sorted(obs.keys()):
                # sometimes stations need not supply a single measurement at a given time, then zip() fails
                if len(obs[tm]) > 0:
                    # remove unavailable measurments
                    var_list, obs_list = zip(*obs[tm])
                    obs_var = [ obs_var_tbl[x] if x in obs_var_tbl else float("nan") for x in var_list ]
                    f.write(tm.strftime('%Y-%m-%d_%H:%M %Z'))
                    f.write('\n')
                    f.write(string.join(map(str, var_list), ", "))
                    f.write('\n')
                    f.write(string.join(map(str, obs_list), ", "))
                    f.write('\n')
                    f.write(string.join(map(str, obs_var), ", "))
                    f.write('\n')


