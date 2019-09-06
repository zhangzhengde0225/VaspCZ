#!/usr/bin/env perl
#only works for triclinic cell, each atom position line must have 7 numbers and the 3,4,5(starts from 0) collums are x,y,z

use FindBin qw($Bin);
use lib "$Bin";
use Vasp;
use Math::Trig;
$fact=180/pi;

@args=@ARGV;
@args>=1 || die "usage: cif2pos.pl <POSCAR>\n";
$inputfilename=$args[0];

while (<>) {
	$_ =~ s/^\s+//;
	@F = split(/\s+/,$_);
	if ($F[0] =~ /^_cell_length_a/) {
		@temp = split(/\(/,$F[1]);
		$la = $temp[0];
	}
	if ($F[0] =~ /^_cell_length_b/) {
		@temp = split(/\(/,$F[1]);
		$lb = $temp[0];
	}
	if ($F[0] =~ /^_cell_length_c/) {
		@temp = split(/\(/,$F[1]);
		$lc = $temp[0];
	}
	if ($F[0] =~ /alpha/) {
		$ang[0] = $F[1];
	}
	if ($F[0] =~ /beta/) {
		$ang[1] = $F[1];
	}
	if ($F[0] =~ /gamma/) {
		$ang[2] = $F[1];
	}
	if (@F == 7) {
		push(@atom_symbols, $F[1]);
		push(@coords,$F[3],$F[4],$F[5]);
	}
}
for($i=0;$i<3;$i++){
	$arad[$i]=$ang[$i]/$fact;
}

$v1[0]=1.0;
$v1[1]=0.0;
$v1[2]=0.0;
$v2[0]=cos($arad[2]);
$v2[2]=sin($arad[2]);
$v2[1]=0.0;
$v3[0]=cos($arad[1]);
$v3[2]=(cos($arad[0])-cos($arad[1])*cos($arad[2]))/sin($arad[2]);
$v3[1]=sqrt(1.0-$v3[0]**2-$v3[1]**2);

$v1[0]*=$la;
$v1[1]*=$la;
$v1[2]*=$la;
$v2[0]*=$lb;
$v2[1]*=$lb;
$v2[2]*=$lb;
$v3[0]*=$lc;
$v3[1]*=$lc;
$v3[2]*=$lc;

$var = 0;
$count[$var] = 1;
$atom_types[$var] = $atom_symbols[0];
for ($i = 1; $i < @atom_symbols; $i++){
	if ($atom_symbols[$i] eq $atom_symbols[$i-1]){
		$count[$var] += 1;
	}else{
		$var += 1;
		$count[$var] = 1;
		$atom_types[$var] = $atom_symbols[$i];
	}
}
########################
# Write POSCAR
######################
@outfile = split(/\./,$inputfilename);
open(OUT, ">$outfile[0]");
print OUT ("@atom_types\n");
print OUT ("1.0"."\n");
print OUT ($v1[0]."\t".$v1[1]."\t"." $v1[2]\n");
print OUT ($v2[0]."\t".$v2[1]."\t"." $v2[2]\n");
print OUT ($v3[0]."\t".$v3[1]."\t"." $v3[2]\n");
print OUT ("@count\n");
print OUT ("Selective Dynamics\n");
print OUT ("Direct\n");
$k = @coords / 3;
$j = 0;
for ($i = 0; $i < $k; $i++){
	print OUT ("$coords[$j]   $coords[$j+1]   $coords[$j+2]   ");
	print OUT ("T "."T "."T \n");
	$j = $j + 3;
}

