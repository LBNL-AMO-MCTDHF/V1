


blah=( `wc Dips/ZDipoleft.Dat.01000_0_0_005000`)

num=${blah[0]}

rm ZDipoleft.all.Dat

line1=`ls Dips/ZDipoleft.Dat.*_0_0_005000`

echo $line1 
echo "hit enter"
read var

for file in $line1 
do

	ext=`echo $file |cut -f 1 -d _ |cut -f 3 -d .`

	echo "EXT $ext"
	echo "FILE $file"

	rm temp.bb
	i=0
	while [[ $i -lt $num ]]
	do
		(( i++ ))
		echo "$ext   " >> temp.bb
	done
	cat $file |colrm 20 125 > temp.cc
	paste temp.bb temp.cc >> ZDipoleft.all.Dat
	echo >> ZDipoleft.all.Dat

done
