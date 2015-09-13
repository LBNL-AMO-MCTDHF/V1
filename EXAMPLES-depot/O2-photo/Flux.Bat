#!/bin/bash

if [[ "$1 " == " " ]]
then
	echo "need time"
	exit
fi

energy=`grep -A 2 CONV Outs/Out.relax |tail -n 1|tr -s " "|cut -f 2 -d " "`

echo "ENERGY $energy"

chmctdhf_diatom Inp=Input.Inp.Fourier$2 Act=16   Eground=$energy T=$1 |tee Outs/Out.flux$2.t$1

cp -p Dat/xsec.spi$2.dat Dat/xsec.spi$2.dat.t$1




