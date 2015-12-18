#!/bin/bash

echo '

set style data lines

set border 31 lw 2

set logscale y

set ytics 2

set key bottom right
set xlabel "Time (as)"

set term unknown

# lebedev quadrature weights.
#  38 point lebedev quadrature (7 points w/o redundancy)
#  Overlaps.some38.dat should be ordered xxy,xxz,xyy,yyz,xzz,yzz, xyz.
#  38 points = 6 (not used) + 24 + 8 (multiply by degeneracy... here c2v)

## not used A = 0.009523809523810 * 2
B = 0.028571428571429 * 4
C = 0.032142857142857 * 8

## print 3 * A + 6 * B + C

plot "Dat/Overlaps.some38.dat" using 1:($3 * B + $8 * B + $13 * B + $18 * B + $23 * B + $28 * B + $33 * C) lw 3 title "2"
replot "Dat/Overlaps.some38.dat" using 1:($4 * B + $9 * B + $14 * B + $19 * B + $24 * B + $29 * B + $34 * C) lw 3 title "3"
replot "Dat/Overlaps.some38.dat" using 1:($5 * B + $10 * B + $15 * B + $20 * B + $25 * B + $30 * B + $35 * C) lw 3 title "4"

set yrange [0.00001:*]
set yrange [0.0000001:0.1]

set term x11
replot

set key top left
set term post solid color 20
set out "overlaps.some38.ps"
replot


set out
set term unknown
set xrange [77.5:*]
set yrange [*:*]

set table
set out "overlaps.some38.table"
replot

' |gnuplot -persist

line=`grep " i" overlaps.some38.table |perl -pi -e "s/ i/ /" |perl -pi -e "s/78 / /"`
echo $line > overlaps.some38.table


