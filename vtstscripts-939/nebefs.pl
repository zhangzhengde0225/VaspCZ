#!/usr/bin/env perl
#;-*- Perl -*-

# Script prints the force, energy etc of OUTCAR's in immediate subdir
# of present working directory. Older versions of this script had specific
# paths hardcoded.  

  use Cwd ;
  $dir = cwd ;      
  @l1=`ls -la $dir/[0-9][0-9]/OUTCAR`;   #specifies location of outcars
  $i = 0 ;
  print " Image   Force        Stress   Volume    Magnet     Rel Energy  \n"; 
  foreach $_ (@l1) {
    chop() ;
    @t1=split() ;
    $t2=$t1[@t1-1] ;
#    $steps  = `grep 'energy  without' $t2 | wc |cut -c 0-8` ; 
    $energy = `grep 'energy  without' $t2 | tail -n 1 |cut -c 67-78` ; 
    $force =  `grep 'max\ at' $t2 | tail -n 1 |cut -c 27-38` ;
    $stress =  `grep 'Stress total and'  $t2 | tail -n 1 | cut -c 48-55` ;
    $volume =  `grep 'volume'  $t2 | tail -n 1 | cut -c 24-30` ;
    $mag =  `grep 'number of electron' $t2 | tail -n 1 | cut -c 51-59` ;
#    chop($steps) ;
    chop($energy) ;
    chop($force) ;
    chop($stree) ;
    if(!$i) {$e0 = $energy ;}
    $rel = $energy - $e0 ;
    #@f4 = ($i,$force,$energy,$rel) ;
    #printf "%4i %16.8f %16.8f %16.8f \n",@f4 ; 
    @f4 = ($i,$force,$stress,$volume,$mag,$rel) ;
    printf "%4i %12.8f %12.8f %6.2f %11.8f %12.8f \n",@f4 ; 

    $i++ ;

  };

