#!/bin/sh
# 
YYYY=2014
mm=09
iday=04
here='/home/metno/carrasco/MYWAVE_WAM/WAM10/pre'
wdir='/work/carrasco/MyWaveWAM/WAM10/'
ecspedir='/work/carrasco/MyWaveWAM/ECspect/Sep4/'
ecname='BWP'

idate=$mm$iday'0000'
cd $here
for dd in 04 05 06 
    do
    for each in 00 03 06 09 12 15 18 21
        do 
        edate=$mm${dd}${each}'00'
     
	./make_boundary.sh ${ecspedir}${ecname}${idate}${edate}'1'  $YYYY${edate}'00'
        #mv Boundary/TMP $wdir/TMP$YYYY$edate'00'   #TMP* are ascii

        mv  'BXX'$YYYY$edate'00' $wdir             #BXX* are binary files
	done
done

 for dd in 07
    do
    for each in 00 03 06
        do 
        edate=$mm${dd}${each}'00'
     
	./make_boundary.sh ${ecspedir}${ecname}${idate}${edate}'1'  $YYYY${edate}'00'
        #mv Boundary/TMP $wdir/TMP$YYYY$edate'00'   #TMP* are ascii

        mv  'BXX'$YYYY$edate'00' $wdir              #BXX* are binary files
	done
done

     
