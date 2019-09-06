#!/usr/bin/env perl
#;-*- Perl -*-

# Calculates the distance, with periodic boundary conditions,
# from a specific atom to all other atoms in the file.

use FindBin qw($Bin);
use lib "$Bin";
use Vasp;

# variables needed by the script

@args = @ARGV;
@args==2 || die "usage: neighbor.pl <POSCAR filename> <central atom>\n";

# Open up and read the poscar file.

$poscarfile = $args[0];
$central_atom = $args[1];

($coordinates,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective)
  = read_poscar($poscarfile);

# Calculate the distance from atom # $ARGV[1] to all other atoms and
# sort them in ascending order.

open OUT , ">out.tmp";
for ($i=0; $i<$total_atoms; $i++) {
    for ($j=0; $j<3; $j++) {
        $cart_coord->[0][$j] = $coordinates->[$i][$j];
        $difference->[0][$j] = pbc($coordinates->[$central_atom-1][$j]-$coordinates->[$i][$j]);
    }
    $difference = dirkar($difference,$basis,$lattice,"1");
    $mag_difference = magnitude($difference,"1");
    $cart_coord = dirkar($cart_coord,$basis,$lattice,"1");
    printf OUT "%4d %13.7f %13.7f %13.7f    ... %13.7f\n",$i+1,
            $cart_coord->[0][0],$cart_coord->[0][1],$cart_coord->[0][2],$mag_difference;
}
close OUT;

# Here we get the distances and index them.
open IN , "out.tmp" ;
$i = 0 ;
while ($line=<IN>) {
    $file .= $line;
    chomp($line);
    $line =~ s/^\s+//g;
    @line = split /\s+/,$line;
    $d[$i] = $line[5];
    $indx[$i] = $i;
    $i++;
}
close IN;
unlink "out.tmp";

# Sort with selection sort.
for ($i=0; $i<$total_atoms; $i++) {
    $iptr = $i;
    for ($j=$i+1; $j<$total_atoms; $j++) {
        if($d[$j] < $d[$iptr]) {
            $iptr = $j;
        }
    }
    if ($i != $iptr) {
        $tmp = $d[$i];
        $d[$i] = $d[$iptr];
        $d[$iptr] = $tmp;
        $tmp = $indx[$i];
        $indx[$i] = $indx[$iptr];
        $indx[$iptr] = $tmp;
    }
}

# Print out the sorted data.
open OUT, ">neighdist.dat";
@file = split /\n/ , $file;
for ($i=0; $i<$total_atoms; $i++) {
    print OUT $i,"   ",$file[$indx[$i]],"\n";
}
close OUT;

