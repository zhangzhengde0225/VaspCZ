#!/usr/bin/env perl
#;-*- Perl -*-

# This program calculates the radial distribution function of an atom in 
# a POSCAR file.  The bin size is an argument.

use FindBin qw($Bin);
use lib "$Bin";
use Vasp;

# variables needed by the script

@args = @ARGV;
@args==3 || die "usage: pos2rdf.pl <POSCAR filename> <atom to measure from> <bin size in Angstroms> ";

$poscarfile = $args[0];
$central_atom = $args[1];
$bin_size = $args[2];

($coordinates,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective)
  = read_poscar($poscarfile);

for ($i=0; $i<$total_atoms; $i++) {
    if ($i+1!=$central_atom) {
        for ($j=0; $j<3; $j++) {
            $difference->[0][$j] = pbc($coordinates->[$central_atom-1][$j]-$coordinates->[$i][$j]);
        }
# new WS pbc
        $difference = pbc_difference_ws($difference,$basis,"1");
# end WS pbc
        $difference = dirkar($difference,$basis,$lattice,"1");

        $mag_difference = magnitude($difference,"1");
        $index = int($mag_difference/$bin_size);
        $bin->{$index}[@{$bin->{$index}}] = $i;
    }
}

print "----------------------------------------------------------------------\n";
print "RDF of atom ".$central_atom." with bin size of ".$bin_size.".\n";
print "Distance: Neighbors\n";
print "----------------------------------------------------------------------\n";

foreach $index (sort {$a<=>$b} keys %$bin) {
    $number = @{$bin->{$index}};
    printf "%5.3f (%3i) : ", $index*$bin_size,$number;
    for ($i=0; $i<@{$bin->{$index}}; $i++) {
        printf "%4i ",($bin->{$index}[$i]+1);
    }
    print "\n";
}

