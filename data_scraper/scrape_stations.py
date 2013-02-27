#!/usr/bin/env python

from bs4 import BeautifulSoup
import sys
import urllib
import time
import argparse
from datetime import datetime
import string
import os

# mesowest list of stations
mesowest_station_url = 'http://mesowest.utah.edu/cgi-bin/droman/download_ndb.cgi?stn=%s'
mesowest_net_url = 'http://mesowest.utah.edu/cgi-bin/droman/stn_mnet.cgi?mnet=%d'
mesowest_dl_url = 'http://mesowest.utah.edu/cgi-bin/droman/meso_download_mesowest_ndb.cgi'
mesowest_station_pos = 'http://mesowest.utah.edu/cgi-bin/droman/side_mesowest.cgi?stn=%s'

# format of timestamp
tstamp_fmt = '%Y-%m-%d_%H:%M'

# renaming dictionary handling differences between variable names in xls files from Mesowest
# and in the variable listing
renames = { 'TMPF' : 'TMP', 'DWPF' : 'DWP' }


def retrieve_web_page(addr):
    """
    Retrieve the web page and return it as a string.
    """
    f = urllib.urlopen(addr)
    content = f.read()
    f.close()
    return content


def extract_stations(table):
    """
    Go through rows of a table found on the HTML page listing stations
    on the MesoWest website and extract station data.  The station data
    is returned as a list of dicts, which contain the station info. 
    """
    stations = []

    name = [ 'code', 'name', 'state', 'rep24', 'status' ]

    # find all rows in the table
    rows = table.find_all('tr')
    for r in rows:
        cells = r.find_all('td')
    	if len(cells) == 5:
	        # we have five cells, retrieve the text and add to list
            val = [ c.getText() for c in cells ]
    	    stations.append(dict(zip(name,val)))

    	if len(cells) == 6:
    	    # we have six cells, second is WIMS id, we throw that out
            val = [ c.getText() for c in cells ]
            del val[1]
            stations.append(dict(zip(name,val)))
	
    return stations


def get_station_info(station_info):
    """
    Return a list all the variables that are measured by a station with the station
    code given, store it in the station_info dictionary and also return it.
    """
    # retrieve web page, parse it and pause slightly
    p = retrieve_web_page(mesowest_station_url % station_info['code'])
    soup = BeautifulSoup(p)
    table = soup.find_all('table')[-2]
    varc = table.find_all('td')[2]
    vlist = [ inp.get('value') for inp in varc.find_all('input') ]
    station_info['vlist'] = vlist

    # retrieve web page with position info for station
    p = retrieve_web_page(mesowest_station_pos % station_info['code'])
    soup = BeautifulSoup(p)
    data = filter(lambda x: x.find(':') > 0, map(string.strip, soup.div.getText().split('\n')))
    d = dict([ s.split(':') for s in data ])
    station_info['elevation'] = int(d['ELEVATION'][:-3]) * 0.3048
    station_info['lat'] = float(d['LATITUDE'])
    station_info['lon'] = float(d['LONGITUDE'])
    station_info['wims'] = d['WIMS ID']
    station_info['mnet'] = d['MNET']
    station_info['name'] = d['NAME']
	   

def observes_variables(station_info, req_vars):
    """
    Checks if a station observs all variables req_vars.
    Downloads the variable list if not already available.
    """
    if station_info.has_key('vlist'):
        vlist = station_info['vlist']
    else:
        get_station_info(station_info)
        vlist = station_info['vlist']
#        print('# vars for %s : %s' % (station_info['code'], str(station_info['vlist'])))
        time.sleep(0.1)

    return all([v in vlist for v in req_vars])


def find_and_list_stations(args):
    """
    Find stations that satisfy search criteria and return a list of
    dicts with their information.
    """
    # dictionary to convert station names to MesoWest network ids
    network_to_mesowest_id = { 'RAWS' : 2, 'NWS' : 1 }
    net_id = network_to_mesowest_id[args.network]
    
    # retrieve the network page
    content = retrieve_web_page(mesowest_net_url % net_id)
    soup = BeautifulSoup(content)

    # find the second-to-last table, which contains the station list
    # & extract rows into a list of lists
    stations = extract_stations(soup.find_all('table')[-2])

    print('# total stations: %d' % len(stations))

    # if requested, filter stations by name
    if len(args.state) > 0:
        stations = filter(lambda x: x['state'] == args.state, stations)
        print('# stations in %s : %d' % (args.state, len(stations)))

    # only report active stations
    stations = filter(lambda x: x['status'] == 'ACTIVE', stations)
    print('# ACTIVE stations : %d' % len(stations))

    # if requested, filter by variables
    if args.vlist is not None:
        stations = filter(lambda x: observes_variables(x, args.vlist), stations)
        print('# stations yielding %s: %d' % (str(args.vlist), len(stations)))

    return stations


def download_station_data(station_info, out_fmt, tstmp, length_hrs, vlist = None):
    """
    Downloads station observations for the given station and parameters.
    """
    output_map = { 'xls' : 'Excel', 'csv' : 'csv', 'xml' : 'xml' }
    
    params = [ ['product', ''],                 # empty in the form on the website
               ['stn', station_info['code']],   # the station code
               ['unit', '1'],                   # metric units
               ['time', 'GMT'],                 # tstamp will be in GMT
               ['day1', tstmp.day],             # the end timestamp of the measurements
               ['month1', '%02d' % tstmp.month],
               ['year1', tstmp.year],
               ['hour1', tstmp.hour],
               ['hours', length_hrs],
               ['daycalendar', 1],
               ['output', output_map[out_fmt]], # output format (excel/csv/xml)
               ['order', '0']                   # order is always ascending
               ]

    # append all variables
    var_list = vlist if vlist is not None else station_info['vlist']
    if not observes_variables(station_info, var_list):
        return None

    for v in var_list: 
        params.append([v, v])

    # join al internal parameters
    get_rq = mesowest_dl_url + '?' + string.join(map(lambda x: x[0] + '=' + str(x[1]), params), '&')

    # download the observed variables
    content = retrieve_web_page(get_rq)
    return content 


def parse_dt(dt):
    tstmp = datetime.strptime(dt, tstamp_fmt)
    return tstmp

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Downloads and processes fuel moisture data from the MesoWest website.')
    parser.add_argument('command', type=str, help='the command to execute, one of [list|info|dl]')
    parser.add_argument('-n', '--network', metavar='N', type=str, dest='network', help='the name of the station network to query (one of RAWS, NWS), list command only')
    parser.add_argument('-c', '--code', metavar='D', type=str, dest='station_code', help='the code of the station to process, dl/vars command only')
    parser.add_argument('-s', '--state', metavar='S', type=str, dest='state', help='only return stations from a particular state (use state code), list command only')
    parser.add_argument('-v', '--vlist', metavar='V', type=str, dest='vlist', default=None, help='(dl/list) download the var(s) from a station if command is dl or list stations which provide var(s) if list; separate variables by a comma, no spaces')
    parser.add_argument('-i', '--interval', metavar='H', type=int, dest='hours', help='the timespan for which to download measurements (number of hours, dl only)')
    parser.add_argument('-t', '--tstamp', metavar='T', type=parse_dt, dest='tstamp', help='end date/time for the downloaded range, dl only (fmt=Y-m-d_H:M)')
    parser.add_argument('-f', '--fmt', metavar='M', type=str, default='xls', dest='fmt', help='(dl only) download in what format? (xls/csv/xml)')
    parser.add_argument('-l', '--line', action='store_const', const='terse', default='loose', dest='info_fmt', help='(info only) print information in terse format')
    
    args = parser.parse_args()

    if args.vlist is not None:
        args.vlist = args.vlist.split(',')

    if args.command == 'list':
        stations = find_and_list_stations(args)
        map(lambda x: sys.stdout.write(x['code'] + '\n'), stations)

    elif args.command == 'info':
        si = { 'code' : args.station_code }
        get_station_info(si)
        if args.info_fmt == 'loose':
            print("# Info created by scrape_stations.py on %s" % str(datetime.now()))
            print(args.station_code)
            print("# Station name")
            print(string.strip(si['name']))
            print("# Station geo location")
            print("%g, %g" % (si['lat'], si['lon']))
            print("# Elevation (meters)")
            print("%g" % si['elevation'])
            print("# Station sensors")
            nnv = [ renames[x] if x in renames else x for x in si['vlist'] ]
            print(string.join(nnv, ", "))
        else:
            print(si['code'] + ',' + str(si['lat']) + ',' + str(si['lon']) + ',' + str(si['elevation']))

    elif args.command == 'dl':
        station_info = { 'code' : args.station_code }
        get_station_info(station_info)
        doc = download_station_data(station_info, args.fmt, args.tstamp, args.hours, args.vlist)
        
        with open('%s_%s.%s' % (args.station_code, args.tstamp.strftime(tstamp_fmt), args.fmt), 'w') as f:
            f.write(doc)
    else:
        sys.exit(1)
        
