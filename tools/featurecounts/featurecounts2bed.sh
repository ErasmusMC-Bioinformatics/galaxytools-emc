#!/bin/bash

# featurecounts2bed - converts featureCounts output to BED format

# Copyright 2013-2014, Youri Hoogstrate

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License at <http://www.gnu.org/licenses/> for
# more details.

# This tool has been written by Youri Hoogstrate from the Erasmus
# Medical Center (Rotterdam, Netherlands) on behalf of the Translational
# Research IT (TraIT) project:
# http://www.ctmm.nl/en/programmas/infrastructuren/traitprojecttranslationeleresearch


exon_level="true"
filename=""

# Parse parameters
while getopts e:f: option
do
	case "${option}"
	in
		e) exon_level=${OPTARG};;
		f) filename=$OPTARG;;
	esac
done

# Convert the file
if [ $filename == "" ]; then
	echo "Usage:"
	echo "  -e [true, false]   true = entry for every exon; false = line for genes first exon"
	echo "  -f                 FILENAME from featureCounts"
else
	while read line; do
		first=${line:0:1}
		if [ $first != "#" ]; then
			columns=($line)
			uid=${columns[@]:0:1}
			if [ $uid != "Geneid" ]; then
				chr=${columns[@]:1:1}
				start=${columns[@]:2:1}
				stop=${columns[@]:3:1}
				direction=${columns[@]:4:1}
				length=${columns[@]:5:1}
				count=${columns[@]:6:1}
				
				chr_splitted=($(echo $chr | tr ";" "\n"))
				start_splitted=($(echo $start | tr ";" "\n"))
				stop_splitted=($(echo $stop | tr ";" "\n"))
				strand_splitted=($(echo $direction | tr ";" "\n"))
				
				if [ $exon_level == "true" ]; then
					n=${#chr_splitted[@]}
				else
					n=1
				fi
				
				for (( i=0; i<$n; i++ ))
				do
					echo ${chr_splitted[@]:$i:1}"	"${start_splitted[@]:$i:1}"	"${stop_splitted[@]:$i:1}"	"$uid" ("$((${stop_splitted[@]:$i:1}-${start_splitted[@]:$i:1}))"/"$length"nt)	"$count"	"${strand_splitted[@]:$i:1}
				done
			fi
		fi
	done < $filename
fi
