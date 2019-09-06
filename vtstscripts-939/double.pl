#!/usr/bin/env perl
#;-*- Perl -*-

# double the cell along any vector 
use FindBin qw($Bin);
use lib "$Bin";
@args=@ARGV;
(@args==1 || @args==2) || die "usage: double.pl <inputfile> <outputfile>\n";
$inputfile = @args[0];
if(@args == 2) {
    $outputfile = @args[1];
} else {
    $outputfile = "ciPOSCAR";
}

open (IN,"<$inputfile");
@poscar = <IN>;
close(IN);

# prompt USER for file format
print "SELECT DIRECTION TO DOUBLE CELL  1=v1 or 2=v2 3=v3:  ";
chomp($vec_ = <STDIN>);
#$filetype_ = 1;
if ($vec_ == 1) { $filetype= "vector 1"; }
if ($vec_ == 2) { $filetype= "vector 2"; }
if ($vec_ == 3) { $filetype= "vector 3"; }
if ($vec_ == 1 or $vec_ == 2 or $vec_ == 3) {
  print "     DOUBLING THE CELL ALONG '$vec' \n";
  print "\n";
} else { die "     YOU MUST ENTER 1 or 2 or 2  FOR DOUBLING\n"; }

open (OUT,">$outputfile");

for($index=0; $index<5; $index++){
    $out = $poscar[$index];
    print OUT "$out";
}

$line = $poscar[5];
$line =~ s/^\s+//;
@line = split(/\s+/,$line);
#  if ($line =~ /^s/i) {
if ($line[0] =~ /^\d+$/) {
    $index = 5;
} else {
    $atomtypeflag = 1;
    $out = $poscar[$index];
    print OUT "$out";
    $index = 6;
}

$natoms = $poscar[$index];
$natoms =~ s/^\s+//g;
@natoms = split(/\s+/,$natoms);
$totatoms = 0;
$natomtypes = 0;
for($i=0; $i<@natoms; $i++) {
    $natomtypes++;
    $totatoms += $natoms[$i];
    @natoms[$i] *= 2;
}
$natoms=join("  ",@natoms);
print OUT "$natoms\n";
$index+=1;

for($line=$index; $line<$index+2; $line++) {
    $out = $poscar[$line];
    print OUT "$out";
}
$index = $line;

for($i=$index; $i<$index+$totatoms; $i++) {
    $_ = $poscar[$i];
    $_ =~ s/^\s+//g;

    @out = split(/\s+/,$_);
    $out = join("  ",@out);
    print OUT "$out\n";

    @out = split(/\s+/,$_); 
    @out[$vec_-1] += 1;   # which vector to double along ($vec_-1)
    $out = join("  ",@out);
    print OUT "$out\n";
}

close(OUT);
system "$Bin/pos2con.pl $outputfile tmp.con" ;

open (IN2, "<tmp.con");
@confile = <IN2>;
close(IN2);

open (OUT2, ">tmp2.con");
for($index=0; $index<2; $index++) {
    $out = $confile[$index];
    print OUT2 "$out";
}
# double vector in con format
$_ = $confile[2];
$_ =~ s/^\s+//g;
@out = split(/\s+/,$_);
@out[$vec_-1] *= 2;   # which vector to double along ($vec_-1)
$out = join("  ",@out);
print OUT2 "$out\n";
for($index=3; $index<9+2*$totatoms+2*$natomtypes; $index++) {
    $out = $confile[$index];
    print OUT2 "$out";
}

close(OUT2);
system "$Bin/pos2con.pl tmp2.con $outputfile" ;

system "rm tmp2.con tmp.con" ;
