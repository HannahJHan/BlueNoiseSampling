#!/bin/bash

if [ -z $1 ] || [ -z $2 ] || [ -z $3 ]
then
echo "Usage ./cat_files.sh <input_file0> <input_file2> <output_file>"
exit
fi

echo "Joining $1 (`cat $1 | wc -l` lines) and $2 (`cat $2 | wc -l` lines) ..."

if [ -e $3 ]
then
	echo "Overwriting $3 ..."
	rm $3
fi

paste $1 $2 >> $3