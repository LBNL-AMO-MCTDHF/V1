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
#  26 point lebedev quadrature (7 points w/o redundancy, 4 points w/o x,y,z)
#  Overlaps.some26.dat should be ordered xy,yz,xz, xyz.  
#  26 points = 6 + 12 + 8 (multiply by degeneracy... here c2v)

# NOT USED x,y,z   A = 0.047619047619048 * 2
B = 0.038095238095238 * 4
C = 0.032142857142857 * 8

plot "Dat/Overlaps.some26.dat" using 1:($3 * B + $8 * B + $13 * B + $18 * C) lw 3 title "2"
replot "Dat/Overlaps.some26.dat" using 1:($4 * B + $9 * B + $14 * B + $19 * C) lw 3 title "3"
replot "Dat/Overlaps.some26.dat" using 1:($5 * B + $10 * B + $15 * B + $20 * C) lw 3 title "4"

set yrange [0.00001:*]
set yrange [0.0000001:0.1]

set term x11
replot

set key top left
set term post solid color 20
set out "overlaps.some26.ps"
replot

set out
set term unknown
set xrange [77.5:*]
set yrange [*:*]

set table
set out "overlaps.some26.table"
replot

' |gnuplot -persist

line=`grep " i" overlaps.some26.table |perl -pi -e "s/ i/ /" |perl -pi -e "s/78 / /"`
echo $line > overlaps.some26.table

