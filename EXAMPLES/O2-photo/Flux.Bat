
if [[ "$1 " == " " ]]
then
	echo "need time"
	exit
fi

energy=`grep T= Outs/Out*ax |tail -n 1 |colrm 1 28 |colrm 21 10000 |cut -f 1 -d "E"`

e2=`echo "$energy * 1000" |bc -l`
energy=$e2

echo "ENERGY $energy"

chmctdhf_diatom Inp=Input.Inp.Fourier Act=16   Eground=$energy T=$1 |tee Outs/Out.flux.t$1

cp Dat/xsec.spi.dat Dat/xsec.spi.dat.t$1




