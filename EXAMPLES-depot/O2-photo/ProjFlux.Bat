
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



energy=`grep T= Outs/Out*ax |tail -n 1 |colrm 1 28 |colrm 21 10000 |cut -f 1 -d "E"`

e2=`echo "$energy * 1000" |bc -l`
energy=$e2

echo "ENERGY $energy"


chmctdhf_diatom Inp=Input.Inp.Fourier Act=17   Eground=$energy T=$1 |tee Outs/Out.projflux.t$1

for file in Dat/FTGTau*_*Dat Dat/KVLsum*_*dat Dat/xsec.proj*_*dat
do
	mv $file $file.$2
done

# 



