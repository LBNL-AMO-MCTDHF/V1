
call system(' echo \'  \
# \
# $1 is the directory
# for one file : 
#    $2 (positive real density file input) 
#    $3 (output file) 
#    $4 (xmax) 
#    $5 (time)
# otherwise 
#    $2-$4 (input) 
#    $5 (output file)
#    $6 (xmax)
#    $7 (time)

#echo "DENSITY.BAT : options are XX $1 XX $2 XX $3 XX $4 XX $5 XX $6 XX $7 XX $8 XX $9 XX $10 "

if [[ "$5 " == " " ]]
then
	echo "Need input"
	exit
fi

cd $1

if [[ "$7 " == " " ]]
then

myext=$3
xmax=$4
thistime=$5

#echo "GOING SMALL $1 x $2 x $3 x $4 x $5 x $6 x $7 xx"
#read var


denline=\'

   media {
      emission <1,1,1> / 5
      density {
        	density_file df3 "\'df3/$2.df3\'" 
			interpolate 1
			color_map {
   			[0.000 rgb <0.0,0.0,0>]
   			[0.01 rgb <0.0,0.1,0>]
   			[0.02 rgb <0.0,0.2,0>]
   			[0.04 rgb <0.0,0.3,0>]
   			[0.08 rgb <0.0,0.4,0>]
   			[0.16 rgb <0.0,0.5,0>]
   			[0.25 rgb <0.0,0.6,0>]
   			[0.36 rgb <0.0,0.7,0>]
   			[0.5 rgb <0.0,0.8,0>]
   			[0.7 rgb <0.0,0.9,0>]
   			[0.9  rgb <0.0,1.0,0>]
   			[1.0   rgb <0.0,1.0,0>]
			}
		}
	}	
\'

else

myext=$5
xmax=$6
thistime=$7


#echo "GOING BIG $1 x $2 x $3 x $4 x $5 x $6 x $7 xx"
#read var


denline=\'
   media {
      emission <0,0,1> / 1
      absorption <0,1,0> * 5
      density {
        	density_file df3 "\'df3/$2.df3\'"
			interpolate 1
			color_map {
   			[0.000 rgb <0.0,0.0,0>]
   			[0.01 rgb <0.1,0.1,0.1>]
   			[0.02 rgb <0.2,0.2,0.2>]
   			[0.04 rgb <0.3,0.3,0.3>]
   			[0.08 rgb <0.4,0.4,0.4>]
   			[0.16 rgb <0.5,0.5,0.5>]
   			[0.25 rgb <0.6,0.6,0.6>]
   			[0.36 rgb <0.7,0.7,0.7>]
   			[0.5 rgb <0.8,0.8,0.8>]
   			[0.7 rgb <0.9,0.9,0.9>]
   			[0.9  rgb <1.0,1.0,1.0>]
   			[1.0   rgb <1.0,1.0,1.0>]
			}
		}
       }	
   media {
      emission <0,1,0> / 1
      absorption <1,0,0> * 5
      density {
        	density_file df3 "\'df3/$3.df3\'"
			interpolate 1
			color_map {
   			[0.000 rgb <0.0,0.0,0>]
   			[0.01 rgb <0.1,0.1,0.1>]
   			[0.02 rgb <0.2,0.2,0.2>]
   			[0.04 rgb <0.3,0.3,0.3>]
   			[0.08 rgb <0.4,0.4,0.4>]
   			[0.16 rgb <0.5,0.5,0.5>]
   			[0.25 rgb <0.6,0.6,0.6>]
   			[0.36 rgb <0.7,0.7,0.7>]
   			[0.5 rgb <0.8,0.8,0.8>]
   			[0.7 rgb <0.9,0.9,0.9>]
   			[0.9  rgb <1.0,1.0,1.0>]
   			[1.0   rgb <1.0,1.0,1.0>]
			}
		}
}
   media {
      emission <1,0,0> / 1
      absorption <0,0,1> * 5
      density {
        	density_file df3 "\'df3/$4.df3\'"
			interpolate 1
			color_map {
   			[0.000 rgb <0.0,0.0,0>]
   			[0.01 rgb <0.1,0.1,0.1>]
   			[0.02 rgb <0.2,0.2,0.2>]
   			[0.04 rgb <0.3,0.3,0.3>]
   			[0.08 rgb <0.4,0.4,0.4>]
   			[0.16 rgb <0.5,0.5,0.5>]
   			[0.25 rgb <0.6,0.6,0.6>]
   			[0.36 rgb <0.7,0.7,0.7>]
   			[0.5 rgb <0.8,0.8,0.8>]
   			[0.7 rgb <0.9,0.9,0.9>]
   			[0.9  rgb <1.0,1.0,1.0>]
   			[1.0   rgb <1.0,1.0,1.0>]
			}
		}
}

\'


fi


echo "
Output_File_Type=C
Width=600
Height=600
Quality=9
Antialias=on
Antialias_Threshold=0.1
Display=off
Library_Path=/home/cluster/pbourke/rendering/povinclude

Initial_Frame = 0
Final_Frame = 0
Initial_Clock = 0
Final_Clock = 1
" > temp.ini

echo \'
#include "metals.inc"
#include "colors.inc"

#declare XMAX=\'"$xmax"\';
#declare TWOPI = 6.283185307179586476925287;
#declare RADIUS = 1;
#declare NX = 100;
#declare NY = 100;
#declare SPHRAD=1;
#declare THICK=20;
#declare NZ = 100;
#declare VVPP = NY * 4/4*12/13;
#declare VVTT = NY * 4/4*5/13;
#declare VVCC = NY * 1/25;
#declare DD = <NX,NY,NZ>;
#declare TEXTTRANS = <-NX/50, NY/10, -NZ/5>;
#declare CC = DD / 2;
#declare VXXP = <VVCC,VVTT,VVPP>;
#declare VP = <0.0,0.0,VVPP>;
#declare NUCONE = <NX*(0.5/XMAX)-SPHRAD/2,-SPHRAD/2,-SPHRAD/2>;
#declare NUCTWO = <NX*(-0.5/XMAX)-SPHRAD/2,-SPHRAD/2,-SPHRAD/2>;


global_settings { 
	ambient_light <500,500,500> 
	assumed_gamma 0.68
}


camera {
   location VP
   up y
   right x
   angle 60
   sky <0,0,-1>
   look_at <0,0,0>
}

light_source {
   VP + <0,0,NZ/2>
   color rgb <1,1,1>/5
   media_interaction on
   media_attenuation on
   shadowless
}

light_source {
   VP - <0,0,NZ/2>
   color rgb <1,1,1>/5
   media_interaction on
   media_attenuation on
   shadowless
}

text {
ttf "timrom.ttf" "T=\'"$thistime"\' fs" 3, 0
texture {
pigment { White } 
finish {F_MetalB}
}
rotate <205,0,0>
translate <-3,11,-20> 
scale DD
}


box {
   <0,0,0>, <1,1,1>
   pigment { rgbf 1 }
   interior { /'"$denline"/' }
   hollow
	translate <-0.5,-0.5,-0.5>
	scale DD
}



sphere { 
NUCONE, SPHRAD 
texture {
pigment { White } 
finish {F_MetalB}
}
}

sphere { 
NUCTWO, SPHRAD 
texture {
pigment { White } 
finish {F_MetalB}
}
}

/' > temp.pov


povray +V +Itemp.pov temp.ini
mv temp.tga $myext.tga


cd ..

\' > Density.Bat')

