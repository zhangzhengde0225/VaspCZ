  set terminal postscript enhanced "Helvetica" 18
  set output "res_ins.eps"
  set xlabel "Force calls"
  set nokey
  set grid
#  set style line 1 lt 2 lw 1.5 pt 7 ps 0.8

  set ylabel "Potential energy [eV]"
  plot "o.u.t.t.e.m.p" u 2:3 w lp lt 3 lw 1.5 pt 6 ps 0.8

  set ylabel "Spring energy [eV]"
  plot "o.u.t.t.e.m.p" u 2:4 w lp lt 2 lw 1.5 pt 6 ps 0.8

  set ylabel "|F_{max}| [eV/A]"
  plot "o.u.t.t.e.m.p" u 2:5 w lp lt 2 lw 1.5 pt 6 ps 0.8

  set ylabel "Imaginary mode [eV/A^{2}]"
  plot "o.u.t.t.e.m.p" u 2:6 w lp lt 2 lw 1.5 pt 6 ps 0.8

  set ylabel "S_{o} [eV{/Symbol t}]"
  plot "o.u.t.t.e.m.p" u 2:7 w lp lt 2 lw 1.5 pt 6 ps 0.8

  exit
