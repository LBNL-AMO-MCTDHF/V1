#!/bin/bash

evals=9

phase=( xx 0 0.6981 1.3962 2.0943 2.7925 3.4906 4.1887 4.8869 5.5850 )

totevals=`expr $evals \* $evals`
echo "TOT $totevals"

rm BigFourier.Dat

blah=( `wc Dips/ZDipoleft.Dat.01000_0_0_005000` )

num=${blah[0]}

#echo "$num LINES"
#exit

delax=-100

while [[ $delax -lt 1400 ]]
do

(( delax+=100 ))
delay=`expr $delax + 100000|colrm 1 1`
echo "DELAY $delay"


qlist=" "
ii=0
while [[ $ii -lt $evals ]]
do
	(( ii++ ))
	jj=0
	while [[ $jj -lt $evals ]]
	do
		(( jj++ ))
		qlist="$qlist Dips/ZDipoleft.Dat.${delay}_${phase[ii]}_${phase[jj]}_005000"
	done
done

#ls $qlist |more
#echo $qlist
#exit


iline=580

while [[ $iline -lt 780 ]]
do
	(( iline++ ))

	i=0


	for file in $qlist
	do

		(( i++ ))
		line=( `head -n $iline $file |tail -n 1` )

		energy=${line[0]}

#		abs[$i]=`echo ${line[6]} |perl -pi -e "s/E/\*10^/"|perl -pi -e "s/\+//"`

		redip[$i]=`echo ${line[1]} |perl -pi -e "s/E/\*10^/"|perl -pi -e "s/\+//"|bc -l`
		imdip[$i]=`echo ${line[2]} |perl -pi -e "s/E/\*10^/"|perl -pi -e "s/\+//"|bc -l`
		refield[$i]=`echo ${line[3]} |perl -pi -e "s/E/\*10^/"|perl -pi -e "s/\+//" |bc -l |colrm 14 100`
		imfield[$i]=`echo ${line[4]} |perl -pi -e "s/E/\*10^/"|perl -pi -e "s/\+//" |bc -l |colrm 14 100`


	done

	if [[ $i -ne $totevals ]]
	then
		echo "ERROR $i $totevals"
		exit
	fi

	outline=""

	ifour=0
	while [[ $ifour -lt $evals ]]
	do
		(( ifour++ ))
	jfour=0
	while [[ $jfour -lt $evals ]]
	do
		(( jfour++ ))

		realline=" 0 "
		imagline=" 0 "

		ff=0
		i=0
		while [[ $i -lt $evals ]]
		do
			(( i++ ))
		j=0
		while [[ $j -lt $evals ]]
		do
			(( j++ ))
			(( ff++ ))

#			echo ${imdip[ff]} ${redip[ff]}
			realline="$realline + (${redip[ff]} * c(${phase[i]}*($ifour-1) + ${phase[j]}*($jfour-1)) - ${imdip[ff]} * s(${phase[i]}*($ifour-1) + ${phase[j]}*($jfour-1) ))/$totevals*1000"
			imagline="$imagline + (${redip[ff]} * s(${phase[i]}*($ifour-1) + ${phase[j]}*($jfour-1)) + ${imdip[ff]} * c(${phase[i]}*($ifour-1) + ${phase[j]}*($jfour-1) ))/$totevals*1000"

		done
		done

#		echo "REALLINE"
#		echo "$realline"
#		echo "IMAGLINE"
#		echo "$imagline"
#		exit

		realpart=`echo "( $realline )  " | bc -l`
		imagpart=`echo "( $imagline )  " | bc -l`

		outline="$outline  ${realpart:0:12} ${imagpart:0:12}"

	done
	done

	echo " $delay $energy     ${refield[1]} ${imfield[1]}  $outline" >> BigFourier.Dat


done

echo >> BigFourier.Dat

done


	
