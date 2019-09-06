#!/bin/bash

# AA 06-06-2008

  file=$1
  echo " Input file :: \"$file\" "

# Get the total number of atoms
  set -- $( sed '7q;d' $file )
  sum=0
  while [ $# -gt 0 ]; do
    i=$1
    sum=$( echo "$sum + $i" | bc )
    shift
  done

# Cut out the header info
  nh=$((10+sum))
  rm -f haus
  sed "$nh"q $file > haus

# Find the FFT grid line and use it to split the file
  set -- $( grep -n "$( head -$nh $file | tail -1 )" $file | \
            sed 's/:/ /' | awk '{print $1}' )
  l1=$1 ; l2=$2         # line numbers for the FFT grid info line
  l3=$( wc -l $file | awk '{print $1}' )
  

# Split, and paste the header on top again
  rm -f cf1 cf2
  cp haus cf1
  cp haus cf2
  echo " Total spin (file \"cf1\") "
  sed -n "$((l1+1)),$((l2-1)) p" $file >> cf1
#  head -$((l2-1)) $file | tail -n +$((l1+1)) >> cf1
  echo " Magnetization (file \"cf2\") "
  sed -n "$((l2+1)),$ p" $file >> cf2
#  tail -$((l3-l2)) $file >> cf2

  rm -f haus
  
