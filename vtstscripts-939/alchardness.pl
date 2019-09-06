#!/usr/bin/env perl
#;-*- Perl -*-

# Created by dano 10-01-10
# LAST MODIFIED by dano: 10-01-10
# Reads CONTCAR prints alchemical hardness matrix
######################################################
######################################################


use FindBin qw($Bin);
use lib "$Bin";
use Vasp;
use Math::Trig;
use Math::MatrixReal ;
$fact=180/pi;

print "\n";
@args=@ARGV;
#@args>=1 || die "usage: alchardness.pl <CONTCAR> <output>\n";
if(@args>=1){$inputfilename=$args[0];}
else{$inputfilename="CONTCAR";}
if(@args>=2){$outputfilename=$args[1];}
else{$outputfilename="hardness.out";}

($coordinates,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective,$description) = read_poscar($inputfilename);
$basis=rotate_basis($basis);  # this ensures that there is a smooth conversion back from con format 

set_bc($coordinates,$total_atoms);

for ($i=0;$i<$total_atoms;$i++) {
  for ($j=0;$j<$total_atoms;$j++) {
    for ($k=0;$k<3;$k++) {
      #print "$coordinates->[$i][$k]\n";
      $diff->[$j][$k]=pbc($coordinates->[$j][$k]-$coordinates->[$i][$k]);
      #print "$diff->[$j][$k]\n";
    };
  };
  #print "$diff->[0][0]\n";
  $diff=dirkar($diff,$basis,$lattice,$total_atoms);
  #print "$diff->[1][1]\n";
  for ($j=0;$j<$total_atoms;$j++) {
    $dR=0;
    for ($k=0;$k<3;$k++) {
      $dR += $diff->[$j][$k]*$diff->[$j][$k];
    };
    if ($i==$j) {$hard->[$i][$j]=0;}
    else{$hard->[$i][$j] = 1/sqrt($dR);}
  };
};

$B = Math::MatrixReal->new($total_atoms,$total_atoms) ;
open(OUT,">$outputfilename");
for($i=0 ; $i<$total_atoms ; $i++) {
  $B->assign($i+1,$i+1,$hard->[$i][$i]) ;
  for($j=$i+1 ; $j<$total_atoms ; $j++) {

    $B->assign($i+1,$j+1,0.5*($hard->[$i][$j]+$hard->[$j][$i])) ;
    $B->assign($j+1,$i+1,0.5*($hard->[$i][$j]+$hard->[$j][$i])) ;
  }
  for ($j=0;$j<@{$hard->[$i]};$j++) {

    printf OUT "%16.8f", $B->element($i+1,$j+1)."  " ;
    printf  "%16.8f", $B->element($i+1,$j+1)."  " ;
  }
  print OUT "\n" ;
  print  "\n" ;
}
close OUT ;


