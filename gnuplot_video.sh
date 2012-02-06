#!/usr/bin/gnuplot -persist
# 
# File:   gnuplot_video.sh
# Author: ahuebl
#
# Created on 06.02.2012, 10:22
#

set terminal postscript eps color "Helvetica" 20
set grid
set out 'a.eps'

#set xlabel "t in s"
#set ylabel "p_y/p_y_soll - 1"

#MyYRange=3.0e-3

set xrange [0.0:5.0]
set yrange [0.0:8.0]
#set format y "%11.1e"

#set ytics MyYRange/5.0

# with lines
#plot './pos.dat' using 1:2 with points
plot "< grep \"time(`cat ./timestep`)\" ./pos.dat | awk '{if($4 == \"ghost(0)\") print}'" u 1:2 t "Stern" w p pt 2
#plot "< awk '{if($4 == \"ghost(0)\") print}' ./pos.dat" u 1:2 t "Normal" w p pt 2, \
#     "< awk '{if($4 == \"ghost(1)\") print}' ./pos.dat" u 1:2 t "Ghost" w p pt 2

set xlabel "x"
set ylabel "y"
