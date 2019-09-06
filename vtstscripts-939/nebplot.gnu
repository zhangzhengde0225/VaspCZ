set grid
set pointsize 2
set xlabel "Reaction Coordinate [A]"
set ylabel "Energy [eV]"
set nokey
set terminal postscript eps color
set output "mep.eps"
plot "spline.dat"  u 2:3 w l lt 1 lw 2.4 , \
     "neb.dat" u 2:3 w p lt 3 lw 3.0 pt 7 ps 1.3
