from bs4 import BeautifulSoup
import sys
import urllib
import time

mesowest_dl_url = 'http://mesowest.utah.edu/cgi-bin/droman/download_ndb.cgi?stn=%s'

def retrieve_web_page(addr):
    f = urllib.urlopen(addr)
    content = f.read()
    f.close()
    return content

def extract_variable_list(table):
    varc = table.find_all('td')[2]
    return [ inp.get('value') for inp in varc.find_all('input') ]


def extract_stations(table):
    """
    Go through rows of a table found on the HTML page listing stations
    on the MesoWest website.
    """
    stations = []

    # find all rows in the table
    rows = table.find_all('tr')
    for r in rows:

    	# if the row does not have five cells, it's not in the right
	# format for us, skip it
    	cells = r.find_all('td')
	if len(cells) == 5:
	    # we have five cells, retrieve the text and add to list
	    stations.append([ c.getText() for c in cells ])

	if len(cells) == 6:
	    # we have six cells, second is WIMS id, we throw that out
	    s = [c.getText() for c in cells]
	    del s[2]
	    stations.append(s)
	
    return stations


if __name__ == '__main__':

    # open and read in the file which contains the station list
    with open(sys.argv[1], 'r') as f:
    	soup = BeautifulSoup(f.read())

    # find the second-to-last table, which contains the station list
    # & extract rows into a list of lists
    stations = extract_stations(soup.find_all('table')[-2])

    # extract colorado stations and print them
    colorado_stations = filter(lambda x: x[2] == 'CO', stations)

    print([x[0] for x in colorado_stations])
    print("Found %d stations in Colorado in this network." % len(colorado_stations))

    fm_stations = 0
    for sname in [x[0] for x in colorado_stations]:

        # retrieve web page, parse it and pause slightly
    	p = retrieve_web_page(mesowest_dl_url % sname)
	soup2 = BeautifulSoup(p)

	svars = extract_variable_list(soup2.find_all('table')[-2])
	if 'FM' in svars:
	    print('Station %s yields FM data.' % sname)
            fm_stations += 1

	# 0.5s pause
	time.sleep(0.5)

    print('Found %d stations yielding FM data.' % fm_stations)
	
