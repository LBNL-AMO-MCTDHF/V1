#!/bin/bash

if [[ "$1 " == " " ]]
then
	echo "need time"
	exit
fi

energy=`grep -A 2 CONV Outs/Out.relax |tail -n 1|tr -s " "|cut -f 2 -d " "`

echo "ENERGY $energy"

export OMP_NUM_THREADS=4

chmctdhf_diatom Inp=Input.Inp.Fourier Act=16   Eground=$energy T=$1 |tee Outs/Out.flux.t$1

cp Dat/xsec.spi.dat Dat/xsec.spi.dat.t$1

