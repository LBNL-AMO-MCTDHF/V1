#!/bin/bash

if [[ "$1 " == " " ]]
then
	echo "need time"
	exit
fi

energy=`grep T= Outs/Out*ax |tail -n 1 |colrm 1 28 |colrm 21 10000 |cut -f 1 -d "E"`

e2=`echo "$energy * 10" |bc -l`
energy=$e2

echo "ENERGY $energy"

export OMP_NUM_THREADS=4



chmctdhf_diatom Inp=Input.Inp.Fourier.perp Act=16   Eground=$energy T=$1 |tee Outs/Out.flux.perp.t$1


cp Dat/xsec.spi.perp.dat Dat/xsec.spi.perp.dat.t$1




