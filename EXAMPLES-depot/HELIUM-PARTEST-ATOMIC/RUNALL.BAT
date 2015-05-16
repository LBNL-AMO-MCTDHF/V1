#!/bin/bash

for opt in .par0 ""
do

 echo
 echo  "XXX OPTION $opt XXX"
 echo

 RUN0.BAT                 $opt
 RUN0.BAT.NOPAR           $opt
 RUN0.BAT.NOPAR.SPARSE    $opt
 RUN0.BAT.SPARSE          $opt

done
