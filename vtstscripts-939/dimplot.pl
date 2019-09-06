#!/usr/bin/env perl
#;-*- Perl -*-

use FindBin qw($Bin);

# Get a single data point per step

open DIM,"DIMCAR"
  or die "No DIMCAR to open\n";

# Strip first line, and read the second, which is the first data point
<DIM>;
$oldline=<DIM>;
$line=$oldline;
$line=~s/^\s+//g;
@line=split(/\s+/,$line);
$line[0]==1 || die "First data line in DIMCAR is not step 1\n";
$step=1;

open(OUT,">dimer.dat");
while($newline=<DIM>){
  $line=$newline;
  $line=~s/^\s+//g;
  @line=split(/\s+/,$line);
  if($line[0]!=$step){
    print OUT $oldline;
    $step++;
    $step==$line[0] || die "Error: non-sequential steps in the DIMCAR\n";
  }
  $oldline=$newline;
}
print OUT $oldline;

system "gnuplot $Bin/dimplot.gnu" ;
