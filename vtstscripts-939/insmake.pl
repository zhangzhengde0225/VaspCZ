#!/usr/bin/env perl
#;-*- Perl -*-

#  use lib "/usr/local/vtsttools" ;
use FindBin qw($Bin);
use lib "$Bin";
use Vasp ;

# Input the saddle point (POSCAR), number of images and displacement length.
  die "Input the saddle point (POSCAR), MODECAR, number of images and displacement length. \n" , 
       if @ARGV != 4 ;

  $pos = $ARGV[0] ;
  $mod = $ARGV[1] ;
  $nim = $ARGV[2] ;
  $l   = $ARGV[3] ;

# Before we continue, then make sure that the number of images is mod(2) = 0 
  die "Number of images should be mod(2) = 0 \n" , if $nim%2 ; 

# Read the POSCAR-type file
  ($dir,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective)
    =read_poscar($pos) ;

# Get te header, i.e. the element symbols from the first POSCAR line
  $header=`head -n 1 $pos` ;
  chop($header) ;

# Convert the direct coordinates to cartesian coordinates
  $kar = dirkar($dir,$basis,$lattice,$total_atoms);

# Need to read in the mode-file and make sure that it is normalized
  open IN , $mod ;
  while (<IN>) {$w .= $_}
  close IN ;
  @w = split /\n/ , $w ;
  $ndof = @w ; 

# $ndof should equal the number of "T's" in $selective
  $nc = 0 ;
  for($i = 0 ; $i < $total_atoms ; $i++) {
    $line = $selective->[$i] ; @line = split /\s/ , $line ; 
    for($j = 0 ; $j < 3 ; $j++) {
      if($line[0] eq "T") {
        $mask->[$i][$j] = $w[$nc];
        $nc++ ;
      }else{
        $mask->[$i][$j] = 0.0 ;
      }  
    }
  }
#  print $ndof,"  ",$nc,"\n" ;
  die "Number of degrees of freedom is not the same in $pos and $mod \n" , if $ndof != $nc ;

# And finally normalize the mode.
  $dp = sqrt(dot_product($mask,$mask,$total_atoms)) ;
#  print $dp,"\n" ;
  if($dp != 1.0) {
    for($i = 0 ; $i < $ndof ; $i++) {
      for($j = 0 ; $j < 3 ; $j++) {
        $mask->[$i][$j]/$dp ; 
      }
    }
  }

# Step along the unstable mode and generate images, in pairs
  $n = $nim / 2 ;
  $dl = $l / ($n - 1) * 2 ; 
  for($i = 0 ; $i < $n ; $i++) {
    $k = $i + 1 ;
    for($ii = 0 ; $ii < $total_atoms ; $ii++) { 
      for($jj = 0 ; $jj < 3 ; $jj++) {
        $t_kar->[$ii][$jj] = $kar->[$ii][$jj] + $l * $mask->[$ii][$jj] ;
      }
    }
#
# PBC ---- make sure that the bbc is the correct one to use. 
#    
    $t_dir = kardir($t_kar,$basis,$lattice,$total_atoms) ; 
    for($ii = 0 ; $ii < $total_atoms ; $ii++) {
      for($jj = 0 ; $jj < 3 ; $jj++) {
        $t_dir->[$ii][$jj] = bbc($t_dir->[$ii][$jj]) ;
      }
    }    
    if($k < 10) {$dir = "0$k" ;}
    elsif ($k < 100) {$dir = "$k" ;} 
    mkdir $dir ; 
    write_poscar($t_dir,$basis,$lattice,$num_atoms,$total_atoms,
                 $selectiveflag,$selective,$header,"$dir/POSCAR") ;    
    $l -= $dl ;
  }

