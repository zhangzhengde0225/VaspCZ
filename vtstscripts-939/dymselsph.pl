#!/usr/bin/env perl
#;-*- Perl -*-

# This program selects the degrees of freedom around a given atom.
# It requires one POSCAR, the central atom number (starting from 1),
# the radius in angstroms of atoms to include, and the size of the
# displacement.

use FindBin qw($Bin);
use lib "$Bin";
use Vasp;

print "----------------------------------------------------------------------\n";
@args = @ARGV;
(@args==4 || @args==5) || die "usage: dymselsph.pl <POSCAR filename> <central atom(s)>  <radius to include>\n                                        <displacement in angstroms>";
$inputfilename = $args[0];
if (@args==4) {
    print "Using 1 central atom\n";
    $centralatom = $args[1];
    $radius = $args[2];
    $displacement = $args[3];
} elsif (@args==5) {
    print "Using 2 central atoms\n";
    $centralatom1 = $args[1];
    $centralatom2 = $args[2];
    $radius = $args[3];
    $displacement = $args[4];
}

($coordinates,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective) = read_poscar($inputfilename);

$scale = $lattice;
#$displacement=$displacement/$scale**2;

$s = "Central Coordinates: ";
$c1 = "Central atom 1: ";
$c2 = "Central atom 2: ";
for ($i=0;$i<3;$i++) {
    if (@args==4) {
        $centralcoords->[$i] = $coordinates->[$centralatom-1][$i];
        $c1 .= $coordinates->[$centralatom-1][$i]." ";
    } elsif (@args==5) {
        $centralcoords->[$i] = pbc(.5*$coordinates->[$centralatom1-1][$i] + .5*$coordinates->[$centralatom2-1][$i]);
        $c1 .= $coordinates->[$centralatom1-1][$i]." ";
        $c2 .= $coordinates->[$centralatom2-1][$i]." ";
    }
    $s .= $centralcoords->[$i]." ";
}
$s .= "\n";
$c1 .= "\n";
$c2 .= "\n";
print $s;
print $c1;
if (@args==5) {
    print $c2;
}

$displacefile = "";
$dof = 0;

for ($i=0; $i<$total_atoms; $i++) {
    for ($j=0; $j<3; $j++) {
        $distance->[0][$j] = pbc($coordinates->[$i][$j] - $centralcoords->[$j]);
    }

    $distance = dirkar($distance,$basis,$lattice,"1");

    $magdistance = 0;
    for ($j=0; $j<3; $j++) {
        $magdistance += ($distance->[0][$j])**2;
    }
    $magdistance = sqrt($magdistance);
    if ($magdistance<=$radius) {
        $displacefile .= "$displacement $displacement $displacement ".($i+1)." $magdistance\n";
        $dof++;
    } else {
        $displacefile .= "0 0 0 ".($i+1)." $magdistance\n";
    }
}

print "----------------------------------------------------------------------\n";
print "$dof atoms were found within a radius of $radius of atom $centralatom, \n";
print "leading to ".($dof*3)." degrees of freedom selected.\n";
print "----------------------------------------------------------------------\n";

open (OUT,">DISPLACECAR");
print OUT $displacefile;
close (OUT);

