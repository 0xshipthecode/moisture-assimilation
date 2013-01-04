#!/usr/bin/env bash

STATIONS=`cat $1`
cat /dev/null > station_infos

for S in $STATIONS ;
do
    echo "Processing $S ..."
    python scrape_stations.py -c $S -l info >> station_infos
    python scrape_stations.py -c $S -v FM,RELH,TMPF,QFLG -i 24 -t $2 dl
done
