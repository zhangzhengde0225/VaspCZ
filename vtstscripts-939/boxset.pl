#!/usr/bin/env perl

use FindBin qw($Bin);
use lib "$Bin";

$Options::usageString = "[POSCAR] [New Lattice Constant]";
$Options::docString = 'Modifies the lattice constant of the given POSCAR';
use Options;
my %options = &Options::parseOptions("h");
@args=@ARGV;
(@args==2) || &Options::usage("Must give [POSCAR] and [New Lattice Constant]");

use Vasp;

$poscarfile1 = $args[0];
($coordinates,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective,$description)
  =read_poscar($poscarfile1);

$karcoords = dirkar($coordinates,$basis,$lattice,$total_atoms);
$lattice_old = $lattice;
$lattice = $args[1];
for ($i=0; $i<3; $i++) {
    for ($j=0; $j<3; $j++) {
        $basis->[$i][$j] *= $lattice/$lattice_old;
    }
}
$coordsnew = kardir($karcoords,$basis,$lattice,$total_atoms);

print "Total atoms: $total_atoms...\n";
print "Lattice: $lattice...\n";
print "Shift: ".$shift[0]."  ".$shift[1]."  ".$shift[2]."\n";

write_poscar($coordsnew,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective,$description,"POSCAR.out");

