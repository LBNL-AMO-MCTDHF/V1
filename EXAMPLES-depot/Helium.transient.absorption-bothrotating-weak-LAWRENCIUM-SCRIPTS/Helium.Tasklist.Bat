#!/bin/bash

rm TASKS.LIST

numphase=9

iline=1;  

mkdir Inps

delax=-100

while [[ $delax -lt 1400 ]]
do

   (( delax+=100 ))

   if [[ $delax -lt 0 ]]
   then
      echo "POSITIVE DELAY ONLY"
      exit
   fi

   delay=`expr $delax + 100000|colrm 1 1`
   echo "DELAY $delay"

   iphase=0
   while [[ $iphase -lt $numphase ]]
   do

   jphase=0
   while [[ $jphase -lt $numphase ]]
   do

      shift0=`echo "$iphase * 2 * 3.14159265 / $numphase" |bc -l`
      shift00=`echo "$jphase * 2 * 3.14159265 / $numphase" |bc -l`

      shift1=$shift0
      shift2=$shift00
      shift3=`echo "$shift0 + 1.570796326795" | bc -l`
      shift4=`echo "$shift00 + 1.570796326795" | bc -l`

      shift=`echo $shift0 + 10 |bc -l |colrm 1 1 | colrm 7 10000`
      shiftx=`echo $shift00 + 10 |bc -l |colrm 1 1 | colrm 7 10000`

      ext="${delay}_${shift}_${shiftx}"

      cp Input.Inp.Pulse.QQQ Inps/Input.Inp.Pulse.$ext
      perl -pi -e "s/QQQ/$delay/" Inps/Input.Inp.Pulse.$ext

      perl -pi -e "s/RR1/$shift1/" Inps/Input.Inp.Pulse.$ext
      perl -pi -e "s/RR2/$shift2/" Inps/Input.Inp.Pulse.$ext
      perl -pi -e "s/II3/$shift3/" Inps/Input.Inp.Pulse.$ext
      perl -pi -e "s/II4/$shift4/" Inps/Input.Inp.Pulse.$ext
      perl -pi -e "s/EEE/$ext/" Inps/Input.Inp.Pulse.$ext

      echo "chmctdhf_atom Inp=Inps/Input.Inp.Pulse.$ext > Outs/Out.pulse.$ext" >> TASKS.LIST

     (( jphase++ ))
   done
     (( iphase++ ))
   done

done













