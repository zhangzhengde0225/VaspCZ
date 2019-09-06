#!/usr/bin/env perl
#;-*- Perl -*-

# Calculates the zero-point energy contribution from the positive curvature modes

$hbar = 6.46538e-2;

die "INPUT THE OMEGA^{2} FILE \n" , if @ARGV != 1;
open OMG , $ARGV[0];
while(<OMG>) { $file .= $_; }
close OMG;
@file = split /\n/ , $file;
$n = @file;

$c = 0;
$sw = 0;
for($i=0; $i<$n; $i++) {
    $w[$i] = $file[$i];
    if( $w[$i] < 0) { $c++; }
    else { $sw += sqrt($w[$i]); }
}
print $c," NEGATIVE EIGENVALUE FOUND ! \n";
print " ZERO-POINT ENERGY CONTRIBUTION FROM THE POSITIVE CURVATURE MODES = \n";
print 0.5 * $hbar * $sw," eV \n";

