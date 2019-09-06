#!/usr/bin/env perl
#;-*- Perl -*-

# This program selects the degrees of freedom to use in the dynamical
# matrix calculation. It requires two POSCARS (the minimum and the
# transition state), the number of atoms to include (degrees of 
# freedom = 3 * Number of atoms), and the displacement size.

use FindBin qw($Bin);
use lib "$Bin";
use Vasp;

# variables needed by the script

@args = @ARGV;
@args==4 || die "usage: dymseldsp.pl <POSCAR 1 filename> <POSCAR 2 filename> <number of atoms to include> <displacement>";

$poscarfile1 = $args[0];
$poscarfile2 = $args[1];
$number_atoms = $args[2];
$displacement = $args[3];

($coordinates1,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective) = read_poscar($poscarfile1);
($coordinates2,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective) = read_poscar($poscarfile2);

for ($i=0; $i<$total_atoms; $i++) {
    for ($j=0; $j<3; $j++) {
        $difference->[0][$j]=pbc($coordinates1->[$i][$j]-$coordinates2->[$i][$j]);
    }
    $difference = dirkar($difference,$basis,$lattice,"1");
    $mag_difference->[$i] = magnitude($difference,"1");
    $displacecar->[$i] = "0 0 0 ".($i+1)." ".$mag_difference->[$i];
}

@{$indices} = (0..$total_atoms);
@{$indices} = sort {$mag_difference->[$b] <=> $mag_difference->[$a]} @{$indices};

$disp_string = "$displacement $displacement $displacement";
for ($i=0; $i<$number_atoms; $i++) {
    $displacecar->[$indices->[$i]] =~ s/0 0 0/$disp_string/;
    $ii = $indices->[$i]+1;
    print $ii." ".$mag_difference->[$indices->[$i]]."\n";
}

open (OUT,">DISPLACECAR");
for ($i=0;$i<$total_atoms;$i++) {
    # print $indices->[$i]." ".$mag_difference->[$indices->[$i]]."\n";
    print OUT $displacecar->[$i]."\n";
}
close (OUT);

