#!/usr/bin/env perl
#;-*- Perl -*-

@ARGV>=1 || die "usage: chgsplit.pl <CHGCAR> \n";

open (IN1,$ARGV[0]) || die ("Can't open file $!");
open (OUT,">CHGCAR_tot");

# read the header of the CHGCAR
for ($i=0; $i<6; $i++) {
    $line1 = <IN1>;
    $header1 .= $line1;
}

$atoms1 = <IN1>;
$header1 .= $atoms1;

@atoms1 = split(/\s+/,$atoms1);

$sum1 += $_ for @atoms1;

print "Atoms in file: ".$sum1."\n";

for ($i=0; $i<$sum1+2; $i++) {
    $header1 .= <IN1>;
}

# start reading the total charge density (spin up plus spin down) part
$points1 = <IN1>;
$header1 .= $points1;

@points1 = split(/\s+/,$points1);
$psum1 = 1;

for ($i=1; $i<@points1; $i++) {
    $psum1 *= $points1[$i];
}

print "Points in total charge density: ".$psum1."\n";

print OUT $header1;

for ($i=0; $i<$psum1/5; $i++) {
    $line1 = <IN1>;
    @line1 = split(/\s+/,$line1);
    printf OUT "%18.11E %18.11E %18.11E %18.11E %18.11E\n",$line1[1],$line1[2],$line1[3],$line1[4],$line1[5];  
}
close(OUT);

# start reading the magnetization density (spin up minus spin down) part
open (OUT,">CHGCAR_mag");

$line1 = <IN1>;
while($line1 != $points1){
    $line1 = <IN1>;
}

$points2 = $line1;
@points2 = split(/\s+/,$points2);
$psum2 = 1;

for ($i=1; $i<@points2; $i++) {
   $psum2 *= $points2[$i];
}

print "Points in magnetization density: ".$psum2."\n";
if ($psum1 != $psum2) {die ("Number of points not same in two parts!");}

print OUT $header1;
for ($i=0; $i<$psum1/5; $i++) {
    $line1 = <IN1>;
    @line1 = split(/\s+/,$line1);
    printf OUT "%18.11E %18.11E %18.11E %18.11E %18.11E\n",$line1[1],$line1[2],$line1[3],$line1[4],$line1[5];
}

close(OUT);
close(IN1);
