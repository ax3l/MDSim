#!/usr/bin/gnuplot -persist
# 
# File:   gnuplot_positions.sh
# Author: ahuebl
#
# Created on 24.01.2012, 17:31
#
#set xlabel "t in s"
#set ylabel "p_y/p_y_soll - 1"

#MyYRange=3.0e-3

#set yrange [-MyYRange:MyYRange]
#set format y "%11.1e"

#set ytics MyYRange/5.0

# with lines
plot './pos.dat' using 1:2 with points

set xlabel "x"
set ylabel "y"

