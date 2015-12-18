#!/bin/bash

echo '

set style data lines

set border 31 lw 2

set key bottom right
set xlabel "Time (as)"

set term unknown

# lebedev quadrature weights.
#  38 point lebedev quadrature (7 points w/o redundancy)
#  Overlaps.all.dat should be ordered x,y,z, xxy,xxz,xyy,yyz,xzz,yzz, xyz.
#  38 points = 6 + 24 + 8 (multiply by degeneracy... here c2v)

A = 0.009523809523810 * 2
B = 0.028571428571429 * 4
C = 0.032142857142857 * 8

print 3 * A + 6 * B + C

plot "Dat/Overlaps.all38.dat" using 1:($2 * A + $7 * A + $12 * A + $17 * B + $22 * B + $27 * B + $32 * B + $37 * B + $42 * B + $47 * C) lw 3 title "1"

set term x11
replot


' |gnuplot -persist

