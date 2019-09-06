#!/bin/sh
if [ "$1" == "-h" ]; then
    echo "usage: mdplot.sh [OUTCAR] [REFERENCE ENERGY]"
    echo "       outputs a gnuplot script that plots the energy and"
    echo "       temperature vs time"
    echo 
    echo "example: mdplot.sh | gnuplot -p"
    echo
    exit 0
fi
if [ -z "$1" ]
then
    fn=OUTCAR
else
    fn=$1
fi

if [ -z "$2" ]
then
    eref=0
else
    eref=$2
fi

timestep=`awk '/POTIM/ && $3 ~ /[0-9].*/ {print $3}' $fn`
nions=$(awk '$10~/NIONS/{print $NF}' $fn)

title=$(basename $(pwd))
echo "set title \"$title\""
echo 'set ytics nomirror'
echo 'set xlabel "Time (fs)"'
echo 'set y2tics auto'
echo 'set ylabel "Temperature (K)"'
echo 'set y2label "Energy (eV)"'
echo 'plot "-" w l t "Temperature (K)", "-" w l axis x1y2 t "Total Energy per Atom", "-" w l axes x1y1 ls -1 t "Average Temperature"'
awk "/temperature/ {n+=1*$timestep;print n, \$6}" $fn
echo 'e'
awk "/total energy   ETOTAL =/ {n+=1*$timestep;print n, (\$5 - $eref)/$nions}" $fn
echo 'e'
awk "/temperature/ {t+=1*$timestep;T+= \$6;n+=1}END{print 0,T/n;print t,T/n}" $fn
