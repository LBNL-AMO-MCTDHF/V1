#!/bin/bash

echo '

set style data lines

set border 31 lw 2

set logscale y

set ytics 2

set key bottom right
set xlabel "Time (as)"

set term unknown

#  lebedev quadrature weights.
#  50 point lebedev quadrature (10 points w/o redundancy)
#  Overlaps.some26.dat should be ordered xy,yz,xz, xyz, xxyz,xyyz,xyzz  
#  50 points = 6 + 12 + 8 + 24 (multiply by degeneracy... here c2v)

# NOT USED x,y,z  A = 0.012698412698413 * 2

B = 0.022574955908289 * 4
C = 0.021093750000000 * 8
D = 0.020173335537919 * 8

# print 3 * A + 3 * B + C + 3 * D

plot "Dat/Overlaps.some26.dat" using 1:($3 * B + $8 * B + $13 * B + $18 * C + $23 * D + $28 * D + $33 * D) lw 3 title "2"
replot "Dat/Overlaps.some26.dat" using 1:($4 * B + $9 * B + $14 * B + $19 * C + $24 * D + $29 * D + $34 * D) lw 3 title "3"
replot "Dat/Overlaps.some26.dat" using 1:($5 * B + $10 * B + $15 * B + $20 * C + $25 * D + $30 * D + $35 * D) lw 3 title "4"

set yrange [0.00001:*]
set yrange [0.0000001:0.1]

set term x11
replot

set key top left
set term post solid color 20
set out "overlaps.some50.ps"
replot

set out
set term unknown
set xrange [77.5:*]
set yrange [*:*]

set table
set out "overlaps.some50.table"
replot

' |gnuplot -persist

line=`grep " i" overlaps.some50.table |perl -pi -e "s/ i/ /" |perl -pi -e "s/78 / /"`
echo $line > overlaps.some50.table

