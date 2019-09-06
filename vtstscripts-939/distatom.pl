#!/usr/bin/env perl
#;-*- Perl -*-

# 10-12-2002

# Calculates the distance between two atoms.
# They can be in the same file or two different files.
# Needs .con file

# @ARGV[0] :: input .con file.
# @ARGV[1] :: second .con file or number of the 1st atom.
# @ARGV[2] :: 1st or 2nd atom number.
# @ARGV[3] :: 2nd atom number.

$ift0 = index(@ARGV[0],".con") + 1;
$ift1 = index(@ARGV[1],".con") + 1;
unless ($ift0) { die "Need at least one .con file"; }
if ($ift0 && $ift1) {
    print "Read 2 .con files\n";
    open CON , @ARGV[1]; 
    while (<CON>) { $cf2 .= $_; }
    close CON;
    @cf2 = split /\n/ , $cf2;
}
if (!$ift1) { print "Read 1 .con file\n"; }
open CON , @ARGV[0];
while (<CON>) { $cf1 .= $_; }
close CON;
@cf1 = split /\n/ , $cf1;

# Read box side lengths, the number of types of atoms and the number of each type.

$line = @cf1[2];
chomp($line);
$line =~ s/^\s+//g;
@line = split /\s+/,$line;
@box = @line[0..2];
$line = @cf1[6];
chomp($line);
$line =~ s/^\s+//g;
@line = split /\s+/,$line;
$ntyp = @line[0];
$line = @cf1[7];
chomp($line);
$line =~ s/^\s+//g;
@line = split /\s+/,$line;
@natoms = @line[0..$ntyp-1];

$cumsum[0] = @natoms[0];
for ($i=1; $i<$ntyp; $i++) {
    $cumsum[$i] += $cumsum[$i-1] + $natoms[$i];
}

# Calculate the distance. Need first to find the coordinates of the atoms.
# Then calculate the distance using periodic boundary conditions.

if (!$ift1) {
    $ca = @ARGV[1];
    @Rca = &coordinates($ca,@cf1);
    $a = @ARGV[2];
    @Ra = &coordinates($a,@cf1);
} elsif ($ift0 && $ift1) {
    $ca = @ARGV[2];
    @Rca = &coordinates($ca,@cf1);
    $a = @ARGV[3];
    @Ra = &coordinates($a,@cf2);
}

@d = &dpc($Rca[0]-$Ra[0],$Rca[1]-$Ra[1],$Rca[2]-$Ra[2]);
$dist = sqrt($d[0]**2 + $d[1]**2 + $d[2]**2);

# Output

printf "Atom # %3d in @ARGV[0] : %10.7g %10.7g %10.7g \n",
      $ca,$Rca[0],$Rca[1],$Rca[2];
if (!$ift1) {
    printf "Atom # %3d in @ARGV[0] : %10.7g %10.7g %10.7g \n",
            $a,$Ra[0],$Ra[1],$Ra[2];
} elsif ($ift0 && $ift1) {
    printf "Atom # %3d in @ARGV[1] : %10.7g %10.7g %10.7g \n",
            $a,$Ra[0],$Ra[1],$Ra[2];
}
printf "The distance between them is:  %13.10g \n",$dist;

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

sub coordinates {
    my($in,@file) = @_;
    for ($i=0; $i<$ntyp; $i++) {
        if ($in <= $cumsum[$i]) {
            $line = @file[$in+10+2*$i] ; chomp($line) ; $line=~s/^\s+//g;
            @line = split /\s+/,$line;
            last;
        }
    }
    @line[0..2];
}  

sub dpc {
    my @d = @_;
    while($d[0] >  $box[0]/2) { $d[0] -= $box[0]; }
    while($d[0] < -$box[0]/2) { $d[0] += $box[0]; }
    while($d[1] >  $box[1]/2) { $d[1] -= $box[1]; }
    while($d[1] < -$box[1]/2) { $d[1] += $box[1]; }
    while($d[2] >  $box[2]/2) { $d[2] -= $box[2]; }
    while($d[2] < -$box[2]/2) { $d[2] += $box[2]; }
    @d;
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
