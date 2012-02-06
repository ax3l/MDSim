#!/bin/sh

rm -rf video/*

for i in `seq 1 $1`
do
  echo $i > timestep
  ./gnuplot_video.sh
  convert -density 300 a.eps -flatten -quality 95 -resize 800x560 video/$i.png
done

# convert to video
# ...
