#!/bin/bash

# turn a atomic_proj.xml file into a plain digit format
infile="atomic_proj.xml"
outfile="atomic_proj_digitized.out"

if [ $# == 2 ]
then
	infile=$1
	outfile=$2
elif [ $# != 0 ]
then
	echo "usage: $0"
	echo "usage: $0 infile outfile"
	exit
fi

nbnd=`awk '/NUMBER_OF_BANDS/{getline;print $1;exit}' $infile`

# #row = #local orbitals
# #col = #nbnd

date
echo "Start processing $infile ... "
awk -v nbnd=$nbnd 'BEGIN{FS = "," }
/<ATMWFC/{
	for (i = 1; i <= nbnd; i ++ ) {
		getline
		printf "%15.8e %15.8e ", $1, $2
	}
	printf "\n"
}' $infile > $outfile

date
echo "Finish reading $infile"
nl=`cat $outfile | wc -l`
nlh=`echo "scale = 0; $nl / 2" | bc`

# 
echo "Start pasting spin up and down files ..."
sed -n "1, $nlh p" $outfile > tmp.up
sed -n "$((nlh + 1)), $nl p" $outfile > tmp.down
paste tmp.up tmp.down > $outfile
rm tmp.up tmp.down
date
echo "Finish !"

