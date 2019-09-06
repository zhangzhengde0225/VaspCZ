#!/usr/bin/env perl
#;-*- Perl -*-

# 06-12-2004

# Generates as MODECAR (normalized) from two POSCARs to be used as an initial
# guess for the lowest mode for a min-mode following saddle point search 

# Input two POSCAR-type files
  die "Input two POSCAR-type files\n" if @ARGV != 2 ;
  $p1 = $ARGV[0] ;
  $p2 = $ARGV[1] ;

# Read the first POSCAR file
  open P1 , $p1 or die "$p1 not found in this directory\n" ;
  while (<P1>){$file .= $_ ;}
  @p1 = split /\n/ , $file ;
  $n1 = @p1 ;
  close P1 ;

# Read the second POSCAR file
  open P2 , $p2 or die "$p2 not found in this directory\n" ; 
  $file = undef ;
  while (<P2>){$file .= $_ ;}
  @p2 = split /\n/ , $file ;
  $n2 = @p2 ;
  close P2 ;

# Make sure that they are at least the same length
  die "The two files are not of the same length\n" if $n1 != $n2 ;

# Get the basis vectors 
  for ($i=1; $i<5; $i++){
    $line = $p1[$i] ; chomp($line) ; $line=~s/^\s+//g ; @line=split /\s+/,$line ;
    if($i == 1){$latt = $line[0] ;}
    if($i == 2){@ax[0..2] = @line[0..2] ;}
    if($i == 3){@ay[0..2] = @line[0..2] ;}
    if($i == 4){@az[0..2] = @line[0..2] ;}
  }  

# Get the number of different kinds of atoms and the number of each kind
  if ($p1[5] =~ /\d+/) {
      $lineno = 5;
	#print "is vasp 4 \n";
      $vasp5=0;
  }else{
      $lineno = 6;
      #print "is vasp 5 \n";
      $vasp5=1;
  }
  @elements = split /\s+/ , $p1[$lineno] ;
  $nel = @elements ;
  $line = $p1[$lineno] ; chomp($line) ; $line=~s/^\s+//g ; @line=split /\s+/,$line ;
  @not[0..$nel-1] = @line[0..$nel-1] ;
  while ($not[$k] != undef){
  $natoms+=$not[$k++] ;
  }
  if ($p1[7] =~ /Selective/){
	#print "is selective dynamic \n";
	$selectiveDynamics=1;
  }else{
	#print "not selective dynamic\n";
	$selectiveDynamics=0;
  }
# Calculate the difference between each coordinate in each of the files and apply
# periodic boundary conditions
  $sh = 7 + $vasp5+$selectiveDynamics ; 
  for ($i=0; $i<$natoms; $i++){
    $j = $i+$sh ;
    $l1 = $p1[$j] ; chomp($l1) ; $l1=~s/^\s+//g ;
    @l1=split /\s+/,$l1 ;
    # print "$l1 \n";
    $l2 = $p2[$j] ; chomp($l2) ; $l2=~s/^\s+//g ; @l2=split /\s+/,$l2 ;
# pbc
    $dx = $l1[0]-$l2[0] ; while ($dx > 0.5){$dx -= 1.0 ;} while ($dx < -0.5){$dx += 1.0 ;}
    $dy = $l1[1]-$l2[1] ; while ($dy > 0.5){$dy -= 1.0 ;} while ($dy < -0.5){$dy += 1.0 ;}
    $dz = $l1[2]-$l2[2] ; while ($dz > 0.5){$dz -= 1.0 ;} while ($dz < -0.5){$dz += 1.0 ;}
#
# Transform from direct coordinates to cart. coordinates.
    $w[$i][0] = ($dx*$ax[0]+$dy*$ay[0]+$dz*$az[0])*$latt ;
    $w[$i][1] = ($dx*$ax[1]+$dy*$ay[1]+$dz*$az[1])*$latt ;
    $w[$i][2] = ($dx*$ax[2]+$dy*$ay[2]+$dz*$az[2])*$latt ;
  }
#  print "after the loop \n";
# Normalize the mode and print out
  for ($i=0; $i<$natoms; $i++){
    for ($j=0; $j<3; $j++){
      $sum += $w[$i][$j] * $w[$i][$j] ; 
    }
  }
  $sum = sqrt($sum) ;
#  print "sum is $sum \n";
  open MOD , ">MODECAR" ;
  for ($i=0; $i<$natoms; $i++){
    for ($j=0; $j<3; $j++){
      $w[$i][$j] = $w[$i][$j]/$sum ;
    }
    printf MOD "%20.10E %20.10E %20.10E \n",$w[$i][0],$w[$i][1],$w[$i][2] ;
  }
  close(MOD) ;

