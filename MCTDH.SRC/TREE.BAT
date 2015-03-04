
line=""

for file in H2PROJECT/*90 SINCDVR/*90 HEPROJECT/*90
do

	aa=`grep -i "subroutine" $file |grep -v -i "end " |cut -f 2 -d " " |cut -f 1 -d \(`

	bb=`grep -i "function" $file |grep -v -i "end " |cut -f 2 -d " " |cut -f 1 -d \(`

	line="$line $aa $bb"
done

line2=`echo $line |sort`

#echo "$line2"

for word in $line2
do 

grep -i "call $word" *90 |cut -f 1 -d \( | cut -f 2 -d : |perl -pi -e "s/call//"| tr -s " "

done
