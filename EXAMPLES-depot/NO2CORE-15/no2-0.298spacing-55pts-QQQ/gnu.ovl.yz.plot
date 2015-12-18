#!/bin/bash

echo '

set style data lines

set border 31 lw 2

set logscale y

set ytics 2

set key bottom right
set xlabel "Time (as)"

set term unknown

set title "yz"

plot "Dat/Overlaps.yz.dat" using 1:3 lw 4 title "2"
replot "Dat/Overlaps.yz.dat" using 1:4 lw 4 title "3"
replot "Dat/Overlaps.yz.dat" using 1:5 lw 4 title "4"

set yrange [0.00001:*]
set yrange [0.0000001:0.1]

set term x11
replot

set key top left
set term post solid color 20
set out "overlaps.yz.ps"
replot

' |gnuplot -persist


