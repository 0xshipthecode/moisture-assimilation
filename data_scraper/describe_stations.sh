#!/usr/bin/env bash

STATIONS=`grep -v "^#" $1`

for S in $STATIONS ;
do
    echo "Retrieving info for $S ..."
    python scrape_stations.py -c $S info >> $S.info
done
