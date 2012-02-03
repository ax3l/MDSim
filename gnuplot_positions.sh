#!/usr/bin/gnuplot -persist
# 
# File:   gnuplot_positions.sh
# Author: ahuebl
#
# Created on 24.01.2012, 17:31
#

#set terminal postscript eps color "Helvetica" 20
#set grid
#set out 'a.eps'

#set xlabel "t in s"
#set ylabel "p_y/p_y_soll - 1"

#MyYRange=3.0e-3

#set yrange [-MyYRange:MyYRange]
#set format y "%11.1e"

#set ytics MyYRange/5.0

# with lines
#plot './pos.dat' using 1:2 with points
plot "< awk '{if($4 == \"ghost(0)\") print}' ./pos.dat" u 1:2 t "Normal" w p pt 2
#plot "< awk '{if($4 == \"ghost(0)\") print}' ./pos.dat" u 1:2 t "Normal" w p pt 2, \
#     "< awk '{if($4 == \"ghost(1)\") print}' ./pos.dat" u 1:2 t "Ghost" w p pt 2

set xlabel "x"
set ylabel "y"

