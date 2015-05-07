
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


#energy=`echo "$energy + 1.5 / 1.4" |bc -l|colrm 14 1000`
#energy=0.000000001
#echo "ENERGY $energy"
#exit


chmctdhf_diatom Inp=Input.Inp.Fourier Act=16   Eground=$energy T=$1 |tee Outs/Out.flux.t$1

cp Dat/xsec.spi.dat Dat/xsec.spi.dat.t$1




