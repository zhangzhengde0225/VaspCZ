#!/usr/bin/env perl
#;-*- Perl -*-

# 08-31-2008: GH added $Bin

  use Cwd ;
  $home = cwd ; 
  use FindBin qw($Bin);
  use lib "$Bin";
  use Vasp ;

  $blackhole = "/dev/null" ;

# Go over all the images and look for OUTCARs
  $i = 1 ;
  $dir = "0$i" ;
  while (-e $dir) {
    chdir $dir ;
# -> Find the spring energy and So
    if ($dir eq "01") {
      if (-e "insout.dat") {
         $line = `grep ut insout.dat | tail -1` ;
         chomp($line) ; $line=~s/^\s+//g ; @line = split /\s+/ , $line ;
         $sp = $line[4] ;
         $so = $line[7] ;
         $ei = $line[6] ;
      }elsif (-e "insout.dat.gz") {
         system "gunzip insout.dat.gz" ;
         $line = `grep ut insout.dat | tail -1` ;
         chomp($line) ; $line=~s/^\s+//g ; @line = split /\s+/ , $line ;
         $sp = $line[4] ;
         $so = $line[7] ;
         $ei = $line[6] ;
         system "gzip insout.dat" ;
      }else {
        die "No insout.dat or insout.dat.gz file found in $dir \n" ;
      }
    }
    print " --> ",$dir,"\n" ;
# -> Look for an OUTCAR or OUTCAR.gz
    if (-e "OUTCAR") {
      system "$Bin/vef.pl >& $blackhole" ;
      $line = `tail -1 fe.dat` ; 
    }elsif (-e "OUTCAR.gz") {
      system "gunzip OUTCAR.gz" ;
      system "$Bin/vef.pl >& $blackhole" ;
      $line = `tail -1 fe.dat` ;
      system "gzip OUTCAR" ;
    }else {
      die "No OUTCAR or OUTCAR.gz found in $dir \n" ;
    }
    chomp($line) ; $line=~s/^\s+//g ; @line = split /\s+/ , $line ;
    $en[$i-1] = $line[2] ;
    chdir $home ; 
# -> Look for the next image
    $i++ ; 
    if ($i < 10) {
      $dir = "0$i" ;
    }elsif ($i < 100) {
      $dir = "$i" ; 
    }else {
      die "Too many images \n" ;
    }
  }
  $i-- ;
# Calculate the average potential energy
  $ave = 0.0 ;
  for ($k = 0 ; $k < $i ; $k++) {
   $ave += $en[$k] ;     
  }
  $ave /= $i ;
  print "Average images energy ... ",$ave,"\n" ;
  print "Spring energy ........... ",$sp,"\n" ;
  print "So ...................... ",$so,"\n" ;
  print "Lowest eigenvalue ....... ",$ei,"\n" ;
