
line=`cat tree.out.edited`

for word in $line
do 

echo 
echo " ** $word ** "
echo
grep -i "call $word" *90 |cut -f 1 -d \(

done





