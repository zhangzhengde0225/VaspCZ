#!/usr/bin/env perl
#;-*- Perl -*-

# 08-31-2008: GH added $Bin

  use Cwd ;
  $home = cwd ;

  use FindBin qw($Bin);
  use lib "$Bin";
  use Vasp ;

  $blackhole = "/dev/null" ;

# Go over all the subdirectories and create the first timeframe.
  $i = 1 ; 
  $dir = "0$i" ;
  $fc = 1 ; 
  print $fc,"\n" ;
  while (-e $dir) {
    chdir $dir ;
    print " --> ",$dir,"\n" ;
    system "$Bin/pos2con.pl POSCAR >& $blackhole ; $Bin/con2xyz.pl POSCAR.con >& $blackhole" ;
    chdir $home ;
    $i++ ;
    if ($i < 10) {$dir = "0$i" ;}
    elsif ($i < 100) {$dir = "$i" ;}
    else {die "Too many images \n" ;}
  }
  $i-- ; 
  $file = "timeframe_0.xyz" ;
  if ($i < 10) {system "cat [0][1-9]/POSCAR.xyz > $file" ;}
  elsif ($i < 100) {system "cat [0][1-9]/POSCAR.xyz [1-9][0-9]/POSCAR.xyz > $file" ;}

# Trying to drop the extra lines
  $nl = `wc -l $file` ;
  $natoms = $nl / $i - 2 ;
  &writeout ;

# Then extract the step index where the system is moved and generate a timeframe.
  $i = 1 ; 
  $dir = "0$i" ;
  if (-e "$dir/insout.dat") {
    $ut = `grep ut $dir/insout.dat` ;
  }
  @ut = split /\n/ , $ut ;
  $nt = @ut ;
  for ($j = 1 ; $j < $nt ; $j++) {
    $line = $ut[$j] ; 
    chomp($line) ; $line=~s/^\s+//g ; @line = split /\s+/ , $line ; 
    $fc = $line[2] ;
    print $fc,"\n" ;
    $i = 1 ; 
    $dir = "0$i" ;
    while (-e $dir) {
      chdir $dir ;
      print " --> ",$dir,"\n" ;
      system "$Bin/xdat2pos.pl $fc+1" ; 
      rename "POSCAR.out" , "ptmp" ; 
      system "$Bin/pos2con.pl ptmp >& $blackhole" ; 
      system "$Bin/con2xyz.pl ptmp.con >& $blackhole" ; 
      chdir $home ;
      $i++ ;
      if ($i < 10) {$dir = "0$i" ;}
      elsif ($i < 100) {$dir = "$i" ;}
      else {die "Too many images \n" ;}
    }
    $i-- ; 
    $file = "timeframe_$j.xyz" ;
    if ($i < 10) {system "cat [0][1-9]/ptmp.xyz > $file" ;}
    elsif ($i < 100) {system "cat [0][1-9]/ptmp.xyz [1-9][0-9]/ptmp.xyz > $file" ;}
    &writeout ; 
  } 

# Finally cat all the timeframes together
  print $j,"\n" ;
  if ($j < 10) {
    system "cat timeframe_[0-9].xyz > ins_series.xyz" ;
  }
  elsif ($j < 100) {
    system "cat timeframe_[0-9].xyz timeframe_[1-9][0-9].xyz > ins_series.xyz" ;
  }
  elsif ($j < 1000) {
    system "cat timeframe_[0-9].xyz timeframe_[1-9][0-9].xyz timeframe_[1-9][0-9][0-9].xyz > ins_series.xyz" ;
  }

# ========================================================================================== #

  sub writeout {
    open FILE , $file ;
    undef $f ;
    while (<FILE>) {$f .= $_ ;}
    close FILE ;
    @f = split /\n/ , $f ;
    unlink $file ;
    print $file,"\n" ;
    $mepfile = "mep.xyz" ;                     # Add the MEP path to each timeframe
    if (-e $mepfile) {
      $add = `head -1 $mepfile` ;
      undef $mep ; 
      open MEP , $mepfile ;
      while (<MEP>) {$mep .= $_ ; } 
      close MEP ;
      @mep = split /\n/ , $mep ; 
    }
    open OUT , ">$file" ;
    print OUT $i*$natoms+$add,"\n" ;
    print OUT "fc = ",$fc,"\n" ;
    $kk = 1 ;
    $t = 1 ;
    for ($k = 2 ; $k <= $nl ; $k++) {
      if ($kk == $natoms+1) {
      }
      elsif ($kk == $natoms+2) {
        $kk = 0 ;
      }
      else {
        $out[++$t] = $f[$k] ;
        print OUT $out[$t],"\n" ;
      }
      $kk++ ;
    }
    if (-e $mepfile) {
      for ($k = 2 ; $k < @mep ; $k++) {
        print OUT $mep[$k],"\n" ;
      }
    }
    close OUT ;
  }
