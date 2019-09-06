#!/usr/bin/env perl
#;-*- Perl -*-

@args = @ARGV;
@args == 2 || die "usage: chgdiff.pl <reference CHGCAR> <CHGCAR2>\n";

open (IN1,$args[0]) || die ("Can't open file $!");
open (IN2,$args[1]) || die ("Can't open file $!");
open (OUT,">CHGCAR_diff");

for ($i=0; $i<5; $i++) {
    $line1 = <IN1>;
    $line2 = <IN2>;
    $header1 .= $line1;
}

# check whether it is a vasp5 format
$line1 = <IN1>;
$header1 .= $line1;
$line1 =~ s/^\s+//;
@line1 = split(/\s+/,$line1);

if ($line1[0] =~ /^\d+$/) {
    @atoms1 = @line1;
} else {
    $atoms1 = <IN1>;
    $header1 .= $atoms1;
    $atoms1 =~ s/^\s+//;
    @atoms1 = split(/\s+/,$atoms1);
}

$line2 = <IN2>;
$line2 =~ s/^\s+//;
@line2 = split(/\s+/,$line2);

if ($line2[0] =~ /^\d+$/) {
    @atoms2 = @line2;
} else {
    $atoms2 = <IN2>;
    $atoms2 =~ s/^\s+//;
    @atoms2 = split(/\s+/,$atoms2);
}

$sum1 += $_ for @atoms1;
$sum2 += $_ for @atoms2;

print "Atoms in file1: ".$sum1.", Atoms in file2: ".$sum2."\n";

for ($i=0; $i<$sum1+2; $i++) {
    $header1 .= <IN1>;
}
for ($i=0; $i<$sum2+2; $i++) {
   $dummy = <IN2>;
}

$points1 = <IN1>;
$header1 .= $points1;
$points2 = <IN2>;

$points1 =~ s/^\s+//;
$points2 =~ s/^\s+//;
@points1 = split(/\s+/,$points1);
@points2 = split(/\s+/,$points2);

$psum1 = $points1[0]*$points1[1]*$points1[2];
$psum2 = $points2[0]*$points2[1]*$points2[2];

print "Points in file1: ".$psum1.", Points in file2: ".$psum2."\n";

if ($psum1 != $psum2) {die ("Number of points not same in two files!");}

print OUT $header1;

for ($i=0; $i<$psum1/5; $i++) {
    $line1 = <IN1>;
    $line1 =~ s/^\s+//;
    $line2 = <IN2>;
    $line2 =~ s/^\s+//;
    @line1 = split(/\s+/,$line1);
    @line2 = split(/\s+/,$line2);
    for ($j=0; $j<@line1; $j++) {
        $line1[$j] = $line2[$j]-$line1[$j];
    }
#    printf OUT " %18.11E %18.11E %18.11E %18.11E %18.11E\n",$line1[0],$line1[1],$line1[2],$line1[3],$line1[4];
    printf OUT " %18.11E" x @line1 . "\n", @line1;
}

close(IN1);
close(IN2);
close(OUT);
