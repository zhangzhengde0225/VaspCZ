set grid
#set nokey
set pointsize 1.2
#set xlabel "Step"

set terminal postscript eps color

set pointsize 0.8
set ytics nomirror
set y2tics
set xlabel "Number of Ionic Step"
set ylabel "Force [eV/A]"
set y2label "Energy [eV]"
set terminal postscript eps color
set output "dimer.eps"
set log y
plot "dimer.dat" u 1:2 axis x1y1 t "Max force" w lp lt 1 lw 2.4 pt 7 ps 0.9, \
     "dimer.dat" u 1:4 axis x1y2 t "Energy" w lp lt 3 lw 2.4 pt 7 ps 0.9


#set output "force.eps"
#set ylabel "Max ( Abs ( Force ) ) [eV/Ang]"
#plot "outtmp.dat"  u 2 w lp lt 1 lw 3 pt 7

#set output "energy.eps"
#set ylabel "Energy [eV]"
#plot "outtmp.dat"  u 4 w lp lt 1 lw 3 pt 7

#set output "curvature.eps"
#set ylabel "Curvature"
#plot "outtmp.dat"  u 5 w lp lt 1 lw 3 pt 7

#set terminal postscript eps color enhanced
#replot
#set output
#set terminal aqua
