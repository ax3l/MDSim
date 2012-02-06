#!/bin/sh

rm -rf video/*

for i in `seq 1 $1`
do
  echo $i > timestep
  ./gnuplot_video.sh
  timestep=`printf "%.5i" $i`
  convert -density 300 a.eps -flatten -quality 95 -resize 800x560 video/$timestep.jpg
done

# convert to video
ffmpeg -f image2 -b 4000k -qscale 1 -i video/%05d.jpg output.mpg
