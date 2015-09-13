#!/bin/bash

if [[ "$2 " == " " ]]
then
	echo "need time and extension"
	exit
fi

rm WALKS/cation.configlist.BIN
rm Bin/cation.spfs.bin Bin/cation.avector.bin

if [[ ! -e "Bin/spfs.bin.cation.$2" ]]
then
     echo "No exist Bin/spfs.bin.cation.$2"
     exit
fi
if  [[ ! -e "Bin/avector.bin.cation.$2" ]]
then
     echo "No exist Bin/avector.bin.cation.$2"
     exit
fi
if [[ ! -e "WALKS/cation.configlist.BIN.$2" ]]
then
     echo "No exist WALKS/cation.configlist.BIN.$2"
     exit
fi

ln -s cation.configlist.BIN.$2 WALKS/cation.configlist.BIN
ln -s avector.bin.cation.$2 Bin/cation.avector.bin
ln -s spfs.bin.cation.$2 Bin/cation.spfs.bin

energy=`grep -A 2 CONV Outs/Out.relax |tail -n 1|tr -s " "|cut -f 2 -d " "`

echo "ENERGY $energy"

chmctdhf_diatom Inp=Input.Inp.Fourier$3 Act=17   Eground=$energy T=$1 |tee Outs/Out.projflux$3.$2.t$1

for file in Dat/xsec.proj.spi${3}*_*dat
do
	mv $file $file.$2.t$1
done

