#!/bin/bash

# Keep only the digital part of each state projection
# so as to facilitate postprocessing. Each line contains
# only one state and has the format:
# e c_1 index_1 c_2 index_2 ... c_n index_n
# where e is the state's energy, index_n is index of the nth most
# important local orbital in this state, and c_n is its |coefficient| ^ 2
#

# default file names
infile="pdos.out"
outfile="pdos_digitized.out"

if [ $# == 2 ]
then
	infile=$1
	outfile=$2
elif [ $# != 0 ]
then
	echo "usage: $0"
	echo "       $0 infile outfile"
	exit
fi


# use awk to isolate the digits from each text output

sed "s:\*\[#:  :g" $infile | \
sed "s:\]+:  :g" | \
sed "s:    +:     :g" | \
sed "s:psi = ::g" | \
awk '
BEGIN{ns = 0}
/   e = /{
	ns += 1
	printf "%8.5f ", $3
	getline
	while ( $1 != "|psi|^2") {
		printf "%s ", $0
		getline
	}
	printf "%s\n", $3
}
' | \
awk '{
	printf "%8.5f %5.3f %5d ", $1, $NF, (NF-2) / 2
	for ( i = 2; i < NF; i ++ ) if ( i % 2 == 1 ) printf "%5d ", $i
	for ( i = 2; i < NF; i ++ ) if ( i % 2 == 0 ) printf "%5.3f ", $i
	printf "\n"
}' > $outfile


