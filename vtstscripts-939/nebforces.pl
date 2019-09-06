#!/usr/bin/env perl
#;-*- Perl -*-

  use FindBin qw($Bin);
  use lib "$Bin";
  use Vasp;

# Input the two POSCARs/CONTCARS and the name of the directory
  $file1 = $ARGV[0] ;
  $file2 = $ARGV[1] ;
  $ocar = $ARGV[2] ;
 
#  print $file1,"\n" ;
#  print $file2,"\n" ;
#  print $ocar,"\n" ;
#  $IN = <STDIN> ;

# Read the POSCARs/CONTCARs
  ($coordinates1,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective)
    =read_poscar($file1) ;
  ($coordinates2,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective)
    =read_poscar($file2) ;
  $difference = pbc_difference($coordinates1,$coordinates2,$total_atoms) ;
  $cartesian = dirkar($difference,$basis,$lattice,$total_atoms) ;

# Normalize the difference ... $cartesian is now the local tangent
  $norm = sqrt(dot_product($cartesian,$cartesian,$total_atoms)) ;
  for($i=0 ; $i<$total_atoms ; $i++) {
    for($j=0 ; $j<3 ; $j++) {
      $cartesian->[$i][$j] = $cartesian->[$i][$j]/$norm ;
    }
  }

# Read in the OUTCAR to get the forces from the last step, they will be
# projected onto the local tangent
  open OUT , $ocar ;
  while (<OUT>){$ocar .= $_ ;}
  close OUT ;

# Split the OUTCAR so the forces are close to the top in each segment
  @parts = split /TOTAL-FORCE/ , $ocar ; 

# Extract the forces for the final step
  @line = split /\n/ , $parts[@parts-1] ;
  for($i=2 ; $i<$total_atoms+2 ; $i++) {
    $s = $selective->[$i-2] ; $s=~s/^\s+//g ; @s=split /\s+/ , $s ;
    $l = $line[$i] ; $l=~s/^\s+//g ; @l=split /\s+/ , $l ;
    for($j=0 ; $j<3 ; $j++) {
      if($s[$j] eq "T") {$force->[$i-2][$j] = $l[$j+3] ;}
      else {$force->[$i-2][$j] = 0.0 ;}
    }
  }
  $fp = dot_product($cartesian,$force,$total_atoms) ;

  print $fp,"\n" ;

