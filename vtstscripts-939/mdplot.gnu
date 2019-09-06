set grid
set pointsize 2
set xlabel "Ionic Time Steps"
set ylabel "Energy [eV]"
set nokey
set terminal postscript eps color
set output "energy.eps"
plot "md.totale"  u 1 w l lt 1 lw 2.4

set output "pressure.eps"
set ylabel "Pressure [kB]"
plot "md.press"  u 1 w l lt 1 lw 2.4
