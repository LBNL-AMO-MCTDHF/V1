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
#  Overlaps.all.dat should be ordered x,y,z, xy,yz,xz, xyz, xxyz,xyyz,xyzz  
#  50 points = 6 + 12 + 8 + 24 (multiply by degeneracy... here c2v)

A = 0.012698412698413 * 2
B = 0.022574955908289 * 4
C = 0.021093750000000 * 8
D = 0.020173335537919 * 8

print 3 * A + 3 * B + C + 3 * D

plot "Dat/Overlaps.all26.dat" using 1:($3 * A + $8 * A + $13 * A + $18 * B + $23 * B + $28 * B + $33 * C + $38 * D + $43 * D + $48 * D) lw 3 title "2"
replot "Dat/Overlaps.all26.dat" using 1:($4 * A + $9 * A + $14 * A + $19 * B + $24 * B + $29 * B + $34 * C + $39 * D + $44 * D + $49 * D) lw 3 title "3"
replot "Dat/Overlaps.all26.dat" using 1:($5 * A + $10 * A + $15 * A + $20 * B + $25 * B + $30 * B + $35 * C + $40 * D + $45 * D + $50 * D) lw 3 title "4"

set yrange [0.00001:*]
set yrange [0.0000001:0.1]

set term x11
replot

set key top left
set term post solid color 20
set out "overlaps.all50.ps"
replot

set out
set term unknown
set xrange [77.5:*]
set yrange [*:*]

plot "Dat/Overlaps.all26.dat" using 1:($2 * A + $7 * A + $12 * A + $17 * B + $22 * B + $27 * B + $32 * C + $37 * D + $42 * D + $47 * D) lw 3 title "1"
replot "Dat/Overlaps.all26.dat" using 1:($3 * A + $8 * A + $13 * A + $18 * B + $23 * B + $28 * B + $33 * C + $38 * D + $43 * D + $48 * D) lw 3 title "2"
replot "Dat/Overlaps.all26.dat" using 1:($4 * A + $9 * A + $14 * A + $19 * B + $24 * B + $29 * B + $34 * C + $39 * D + $44 * D + $49 * D) lw 3 title "3"
replot "Dat/Overlaps.all26.dat" using 1:($5 * A + $10 * A + $15 * A + $20 * B + $25 * B + $30 * B + $35 * C + $40 * D + $45 * D + $50 * D) lw 3 title "4"

set table
set out "overlaps.all50.table"
replot

' |gnuplot -persist

line=`grep " i" overlaps.all50.table |perl -pi -e "s/ i/ /" |perl -pi -e "s/78 / /"`
echo $line > overlaps.all50.table

