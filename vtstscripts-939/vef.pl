#!/usr/bin/env perl

use FindBin qw($Bin);

$plotflag = 1;
if(@ARGV > 0) { $plotflag = $ARGV[0]; }

$zip = $ENV{'VTST_ZIP'};
if($zip eq ''){ $zip = 'gzip'; }

$outzipped = 0;
if(-e "OUTCAR.gz") { $outzipped = 1; system "gunzip OUTCAR.gz"; }
if(-e "OUTCAR.bz2") { $outzipped = 1; system "bunzip2 OUTCAR.bz2"; }

$forces = `grep 'FORCES: max atom, RMS' OUTCAR`;
$energy = `grep 'energy  without entropy' OUTCAR`;

@forces = split /\n/, $forces;
@energy = split /\n/, $energy;

$num = @forces;

open OUT, ">fe.dat";

for($i=0; $i<$num; $i++){
    $line = $forces[$i] ;
    chomp($line) ;
    $line =~ s/^\s+//g;
    @line = split /\s+/,$line;
    $f = $line[4],"\n";

    $line = $energy[$i]; 
    chomp($line);
    $line =~ s/^\s+//g;
    @line = split /\s+/,$line;
    $e = $line[6],"\n";
    if(!$i) { $e0 = $e; }
    printf OUT "%5i %20.8f %20.6f %20.6g \n",$i,$f,$e,$e-$e0;
    printf "%5i %20.8f %20.6f %20.6g \n",$i,$f,$e,$e-$e0;
}

close OUT;

if($num>1 and $plotflag) {
    system "gnuplot $Bin/vef.gnu";
}

if($outzipped) { system "$zip -9 OUTCAR \&"; }
if($outzipped) { system "$zip -9 OUTCAR"; }

