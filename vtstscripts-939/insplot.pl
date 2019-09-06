#!/usr/bin/env perl
#;-*- Perl -*-

use FindBin qw($Bin);

system "grep ut insout.dat | grep -v itr | cut -c 5-100 > o.u.t.t.e.m.p";
system "gnuplot $Bin/insplot.gnu";
system "rm o.u.t.t.e.m.p";
