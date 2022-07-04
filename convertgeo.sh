#!/bin/bash

# Script to convert the geometry CASE file to be suitable for ParaView
# Milos Ivanovic, 2007.

# Count how much time steps we have counting lines in _time.txt file
COUNT=`wc -l $1_time.txt | awk '{print $1}'`

touch $1_tmpgeo.txt

# Write the contents of .geo file into temporary file $COUNT times
for i in `seq 1 $COUNT`;
do
	cat $1.geo >> $1_tmpgeo.txt
done

# Rename the temporary file file
mv -f $1_tmpgeo.txt $1.geo

