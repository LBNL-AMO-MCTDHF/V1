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
#  Overlaps.all.dat should be ordered x,y,z, xxy,xxz,xyy,yyz,xzz,yzz, xyz.
#  38 points = 6 + 24 + 8 (multiply by degeneracy... here c2v)

A = 0.009523809523810 * 2
B = 0.028571428571429 * 4
C = 0.032142857142857 * 8

print 3 * A + 6 * B + C

plot "Dat/Overlaps.all38.dat" using 1:($3 * A + $8 * A + $13 * A + $18 * B + $23 * B + $28 * B + $33 * B + $38 * B + $43 * B + $48 * C) lw 3 title "2"
replot "Dat/Overlaps.all38.dat" using 1:($4 * A + $9 * A + $14 * A + $19 * B + $24 * B + $29 * B + $34 * B + $39 * B + $44 * B + $49 * C) lw 3 title "3"
replot "Dat/Overlaps.all38.dat" using 1:($5 * A + $10 * A + $15 * A + $20 * B + $25 * B + $30 * B + $35 * B + $40 * B + $45 * B + $50 * C) lw 3 title "4"

set yrange [0.00001:*]
set yrange [0.0000001:0.1]

set term x11
replot

set key top left
set term post solid color 20
set out "overlaps.all38.ps"
replot


set out
set term unknown
set xrange [77.5:*]
set yrange [*:*]

plot "Dat/Overlaps.all38.dat" using 1:($2 * A + $7 * A + $12 * A + $17 * B + $22 * B + $27 * B + $32 * B + $37 * B + $42 * B + $47 * C) lw 3 title "1"
replot "Dat/Overlaps.all38.dat" using 1:($3 * A + $8 * A + $13 * A + $18 * B + $23 * B + $28 * B + $33 * B + $38 * B + $43 * B + $48 * C) lw 3 title "2"
replot "Dat/Overlaps.all38.dat" using 1:($4 * A + $9 * A + $14 * A + $19 * B + $24 * B + $29 * B + $34 * B + $39 * B + $44 * B + $49 * C) lw 3 title "3"
replot "Dat/Overlaps.all38.dat" using 1:($5 * A + $10 * A + $15 * A + $20 * B + $25 * B + $30 * B + $35 * B + $40 * B + $45 * B + $50 * C) lw 3 title "4"

set table
set out "overlaps.all38.table"
replot

' |gnuplot -persist

line=`grep " i" overlaps.all38.table |perl -pi -e "s/ i/ /" |perl -pi -e "s/78 / /"`
echo $line > overlaps.all38.table


