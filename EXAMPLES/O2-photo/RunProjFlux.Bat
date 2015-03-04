
if [[ "$1 " == " " ]]
then
	echo "Need time"
	exit
fi

 for ext in doub quart; do for sym in sig pi; do for ger in u g
 do 

./ProjFlux.Bat $1 $ext.${sym}_${ger} 

done; done; done

 wait

