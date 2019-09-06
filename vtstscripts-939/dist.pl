#!/usr/bin/env perl
#;-*- Perl -*-

use FindBin qw($Bin);
use lib "$Bin";
use Vasp;

# This program prints the difference in each coordinate for two POSCAR files.

# variables needed by the script

@args=@ARGV;
@args==2 || die "usage: dist.pl <POSCAR 1 filename> <POSCAR 2 filename> ";

$poscarfile1 = $args[0];
$poscarfile2 = $args[1];

($coordinates1,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective)
  = read_poscar($poscarfile1);

($coordinates2,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective)
  = read_poscar($poscarfile2);

$difference = pbc_difference($coordinates1,$coordinates2,$total_atoms);
$cartesian = dirkar($difference,$basis,$lattice,$total_atoms);
$distance = magnitude($cartesian,$total_atoms);

print $distance,"\n";
