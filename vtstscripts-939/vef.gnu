set grid
set pointsize 0.8
set ytics nomirror
set y2tics 
set xlabel "Number of Ionic Step"
set ylabel "Force [eV/A]"
set y2label "Energy [eV]"
set terminal postscript eps color 
set log y
set output "vaspout.eps"
plot "fe.dat" u 1:2 axis x1y1 t "Max force" w lp lt 1 lw 2.4 pt 7 ps 0.9, \
     "fe.dat" u 1:3 axis x1y2 t "Energy" w lp lt 3 lw 2.4 pt 7 ps 0.9


#set output
#set terminal x11
