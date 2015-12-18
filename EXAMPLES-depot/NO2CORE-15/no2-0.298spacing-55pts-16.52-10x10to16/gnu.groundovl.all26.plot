#!/bin/bash

echo '

set style data lines

set border 31 lw 2

set key bottom right
set xlabel "Time (as)"

set term unknown

# lebedev quadrature weights.  Overlaps.all.dat should be ordered x,y,z, xy,yz,xz, xyz.  26 point lebedev quadrature (7 points w/o redundancy)
#  26 points = 6 + 8 + 12 (multiply by degeneracy... here c2v)

A = 0.047619047619048 * 2
B = 0.038095238095238 * 4
C = 0.032142857142857 * 8

print 3*A+3*B+C

plot "Dat/Overlaps.all26.dat" using 1:($2 * A + $7 * A + $12 * A + $17 * B + $22 * B + $27 * B + $32 * C) lw 3 title "1"

set term x11
replot

' |gnuplot -persist


