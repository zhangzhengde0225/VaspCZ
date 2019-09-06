#!/usr/bin/env perl
#;-*- Perl -*-

# This program interpolates between two POSCAR files by the given fraction.

use FindBin qw($Bin);
use lib "$Bin";
use Vasp;

@args = @ARGV;
@args == 3 || die "usage: posinterp.pl <POSCAR 1> <POSCAR 2> <fractional distance between>\n";

$poscarfile1 = $args[0];
$poscarfile2 = $args[1];
$fraction = $args[2];

($coordinates1,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective)
  =read_poscar($poscarfile1);

print "Read $poscarfile1...\n";

($coordinates2,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective)
  =read_poscar($poscarfile2);

$description =`head -n 1 $poscarfile1`;
chop($description);

print "Read $poscarfile2...\n";
print "Total atoms: $total_atoms...\n";
print "Lattice: $lattice...\n";

for ($i=0; $i<$total_atoms; $i++) {
    for ($j=0; $j<3; $j++) {
        $interpolated->[$i][$j] = pbc($coordinates1->[$i][$j]+$fraction*
        pbc($coordinates2->[$i][$j]-$coordinates1->[$i][$j]));
    }
}

write_poscar($interpolated,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective,$description,"POSCAR.out");

