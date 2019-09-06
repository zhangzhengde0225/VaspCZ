#!/usr/bin/env perl
#;-*- Perl -*-

# Created by dano 10-27-09
# LAST MODIFIED by dano: 10-27-09
# Reads OUTCAR for electrostatic potential information
# prints to screen and file
######################################################
######################################################


use FindBin qw($Bin);
use lib "$Bin";
use Vasp;
use Math::Trig;
$fact=180/pi;

print "\n";
@args=@ARGV;
#@args>=1 || die "usage: alcesp.pl <OUTCAR> <output>\n";
if(@args>=1){$inputfilename=$args[0];}
else{$inputfilename="OUTCAR";}
if(@args>=2){$outputfilename=$args[1];}
else{$outputfilename="esp.out";}


$NIONS=`grep 'NIONS' -m 1 "$inputfilename"| awk '{print \$12}'`;
#print "$NIONS \n";
$lines=int($NIONS/5.+1);
#$string = "\"%i %s\\n%i %s\\n%i %s\\n%i %s\\n%i %s\\n\"";
#`grep "the norm of the test charge is" -A $lines OUTCAR |tail -n $lines > temp`;
#print "$string \n";
#`awk '{printf($string,\$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10)}' < temp > temp2`;
#`grep "the norm of the test charge is" -A $lines OUTCAR |tail -n $lines | awk '{printf($string,\$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10)}' > $outputfilename`;
# this is how it would look from the command line
# awk '{printf("%i %s\n%i %s\n%i %s\n%i %s\n%i %s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10)}' < temp > temp2;

$string = "\"%s\\n%s\\n%s\\n%s\\n%s\\n\"";
`grep "the norm of the test charge is" -A $lines $inputfilename |tail -n $lines | awk '{printf($string,\$2,\$4,\$6,\$8,\$10)}' > $outputfilename`;

open(OUT,$outputfilename);
@line=<OUT>;
close(OUT);
print @line;
