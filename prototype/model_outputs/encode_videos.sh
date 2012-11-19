#!/usr/bin/env bash

if [ ! -d videos ]; then
    echo "Directory does not exist, creating."
    mkdir videos
fi

for dir in $(find -name "r*" -type d -printf %f\\n)
do
    avconv -qscale 1 -r 20 -b 9600 -i $dir/moisture_model_t%03d.png videos/$dir.mp4
done
