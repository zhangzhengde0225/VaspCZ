#!/usr/bin/env perl
#;-*- Perl -*-

use FindBin qw($Bin);
use lib "$Bin";
use Vasp;

# This program prints the difference in each coordinate for two POSCAR files.

# variables needed by the script

@args = @ARGV;
@args == 2 || die "usage: diffcon.pl <POSCAR 1 filename> <POSCAR 2 filename>";

$poscarfile1 = $args[0];
$poscarfile2 = $args[1];

# convert to poscar files if we have .con files

$conflag = 0;
if(($poscarfile1 =~ /.con\b/) and ($poscarfile2 =~ /.con\b/)) {
    print "Found con files, converting to poscar files\n";
    $conflag = 1;
    $confile1 = $poscarfile1;
    $confile2 = $poscarfile2;
    $poscarfile1 = "tmp_poscarfile1";
    $poscarfile2 = "tmp_poscarfile2";;
    system "$Bin/pos2con.pl $confile1 $poscarfile1";
    system "$Bin/pos2con.pl $confile2 $poscarfile2";
}

($coordinates1,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective)
  = read_poscar($poscarfile1);

($coordinates2,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective)
  = read_poscar($poscarfile2);

if($conflag) {
    system "rm $poscarfile1";
    system "rm $poscarfile2";
}

$dirdiff = pbc_difference($coordinates1,$coordinates2,$total_atoms);
$cardiff = dirkar($dirdiff,$basis,$lattice,$total_atoms);
$distance = magnitude($cardiff,$total_atoms);

# Print distances between atoms

$dist_tot = 0;
$dist_sq = 0;

$tmp = 0;

for ($atm=0; $atm<$total_atoms; $atm++) {
    $dist = 0;
    for ($i=0; $i<3; $i++) {
        $r = $cardiff->[$atm][$i];
        $dist += ($r*$r);
        $vdist[$i] += $r;
        if($i==1) { $tmp += $r; }
    }
    $dist_sq += $dist;
    $dist = sqrt($dist);

   printf("%9.6f   %9.6f   %9.6f :  %9.6f  %3s\n",@{$cardiff->[$atm]},$dist,($atm+1));
}
$dist_sq = sqrt($dist_sq);

printf("----------------------------------------------------\n");
printf("Displacement:  Sum: %9.6f   Vector: %9.6f\n",$dist,$dist_sq);
printf("                 V: %9.6f %9.6f %9.6f\n",@vdist);
printf("----------------------------------------------------\n");
