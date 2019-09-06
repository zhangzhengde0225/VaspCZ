#!/usr/bin/env perl
#;-*- Perl -*-

# This program takes two different eigenvalue files, which have column labelling
# which eigenvalues are imaginary (1 = imaginary, 0 = real).  It calculates the prefactor
# from them according to Vineyard's formula.

@args = @ARGV;
@args == 2 || die "usage: dymprefactor.pl <minimum freq.dat> <transition freq.dat>\n";
$eigenfile1 = $args[0];
$eigenfile2 = $args[1];

open (IN,$eigenfile1);
while (<IN>) {
    $_ =~ s/^\s+//g;
    $file1 .= $_;
}
close (IN);
@eigenfile1 = split(/\n/,$file1);

open (IN,$eigenfile2);
while (<IN>) {
    $_ =~ s/^\s+//g;
    $file2 .= $_;
}
close (IN);
@eigenfile2 = split(/\n/,$file2);

@eigenfile1 == @eigenfile2 || die "Files do not have the same number of eigenvalues!\n";

$numimg1 = 0;
$numimg2 = 0;
$imgigen = 0;
$prodeigen1 = 1;
$prodeigen2 = 1;
$RatioEig = 1;

for ($i=0; $i<@eigenfile1; $i++) {
    @line1 = split(/\s+/,$eigenfile1[$i]);
    @line2 = split(/\s+/,$eigenfile2[$i]);
  
    if($line1[3] == 1) {
        $numimg1++;
    }
    if($line2[3] == 1) {
        $numimg2++;
        $imgigen = $line2[0];
    }
    $val1 = $line1[0];
    $val2 = $line2[0];
  
    $prodeigen1 *= $line1[0];
    $prodeigen2 *= $line2[0];
    $RatioEig *= ($line1[0]/$line2[0]);
}

$numimg1 == 0 || die "Initial state freq.dat has imaginary frequencies.\n";
$numimg2 == 1 || die "Transition state freq.dat does not have only 1 imaginary frequency.\n";

$prefactor = $RatioEig*$imgigen;

print "The prefactor for this system (in units of inverse cm) is $prefactor\n";

$c = 2.99792458e10; 
$prehertz = $prefactor*$c/1e12; 
print " (in units of TerraHertz) is $prehertz\n"; 

