#!/bin/bash

eth=`grep -A 2 CONV Outs/Out.shiftneon|tail -n 1 |tr -s " " |cut -f 2 -d " "`

but=`grep -A 2 CONV Outs/Out.shifthelium|tail -n 1 |tr -s " " |cut -f 2 -d " "`

combine=`grep -A 2 CONV Outs/Out.combine|tail -n 1 |tr -s " " |cut -f 2 -d " "`

sum=`echo "$eth + $but " |bc -l`

echo $sum

echo $combine



