from bs4 import BeautifulSoup
import sys
import urllib
import time
import argparse
from datetime import datetime

# mesowest list of stations
mesowest_dl_url = 'http://mesowest.utah.edu/cgi-bin/droman/download_ndb.cgi?stn=%s'
mesowest_net_url = 'http://mesowest.utah.edu/cgi-bin/droman/stn_mnet.cgi?mnet=%d'

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
    	    stations.append({key: value for (key, value) in zip(name,val)})

    	if len(cells) == 6:
    	    # we have six cells, second is WIMS id, we throw that out
            val = [ c.getText() for c in cells ]
            del val[1]
    	    stations.append({key: value for (key, value) in zip(name,val)})
	
    return stations


def get_station_variables(station_info):
    """
    Return a list all the variables that are measured by a station with the station
    code given, store it in the station_info dictionary and also return it.
    """
    # retrieve web page, parse it and pause slightly
    p = retrieve_web_page(mesowest_dl_url % station_info['code'])
    soup2 = BeautifulSoup(p)
    table = soup2.find_all('table')[-2]
    varc = table.find_all('td')[2]
    vlist = [ inp.get('value') for inp in varc.find_all('input') ]
    station_info['vlist'] = vlist
    return vlist
	   

def observes_variables(station_info, req_vars):
    """
    Checks if a station observs all variables req_vars.
    Downloads the variable list if not already available.
    """
    if station_info.has_key('vlist'):
        vlist = station_info['vlist']
    else:
        vlist = get_station_variables(station_info)
        print('# vars for %s : %s' % (station_info['code'], str(station_info['vlist'])))
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

    # if requested, filter by variables
    if args.vlist is not None:
        stations = filter(lambda x: observes_variables(x, args.vlist), stations)
        time.sleep(0.1)
        print('# stations yielding %s: %d' % (str(args.vlist), len(stations)))

    return stations


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Downloads and processes fuel moisture data from the MesoWest website.')
    parser.add_argument('command', type=str, help='the command to execute, one of [list|vars|dl]')
    parser.add_argument('--network', metavar='N', type=str, dest='network', help='the name of the station network to query (one of RAWS, NWS), list command only')
    parser.add_argument('--code', metavar='D', type=str, dest='station_code', help='the code of the station to process, dl/vars command only')
    parser.add_argument('--state', metavar='S', type=str, dest='state', help='only return stations from a particular state (use state code), list command only')
    parser.add_argument('--vlist', metavar='V', type=str, dest='vlist', help='download the var(s) from a station if command is dl or list stations which provide var(s) in list command, separate stations by a comma, no spaces')
    parser.add_argument('--from', metavar='F', type=datetime, dest='from_date', help='starting date/time for the downloaded range, dl only')
    parser.add_argument('--to', metavar='T', type=datetime, dest='to_date', help='end date/time for the downloaded range, dl only')
    args = parser.parse_args()

    if args.vlist is not None:
        args.vlist = args.vlist.split(',')
        print args.vlist

    if args.command == 'list':
        stations = find_and_list_stations(args)
        map(lambda x: sys.stdout.write(x['code'] + '\n'), stations)
    elif args.command == 'vars':
        vlist = list_station_variables(args)
        map(lambda x: sys.stdout.write(x + '\n'), vlist)
    elif args.command == 'dl':
        doc = download_station_vars(args)
        with open('%s.xls' % args.code, 'w') as f:
            f.write(doc)
    else:
        sys.exit(1)
        
	
