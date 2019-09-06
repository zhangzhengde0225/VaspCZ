#!/usr/bin/env perl
#;-*- Perl -*-

# we are going to get rid of this script

use FindBin qw($Bin);
use lib "$Bin";
use Vasp;
use Cwd;

$dir = cwd;
$zip = $ENV{'VTST_ZIP'} ;
if($zip eq ''){ $zip = 'gzip'; }

@args = @ARGV;
if (@args==0) {
    opendir(DIR,".") or die "couldn't open . ($!)\n";
    @list=readdir(DIR);
    closedir(DIR);

    @directories = grep {-d && /^[0-9][0-9]$/i} @list;
    @directories = sort {$a<=>$b} @directories;
} else {
    @directories = @args;
}

if(!(-d "vaspgr")) {
    system "mkdir vaspgr";
}
  
for ($i=1; $i<@directories-1; $i++) {

    chdir $directories[$i];
    print "In ",$directories[$i]," ... ","\n";

#    (-e "OUTCAR") || die "No OUTCAR in $directories[$i]!\n" ;

    $zipped = 0;
    if(-e "OUTCAR") { ; }
    elsif(-e "OUTCAR.gz"){ $zipped = 1; system "gunzip OUTCAR.gz"; }
    elsif(-e "OUTCAR.bz2"){ $zipped = 1; system "bunzip2 OUTCAR.bz2"; }
    else{ die "No OUTCAR in $directories[$i]!\n" ;} 

    $energy=`grep 'energy  w' OUTCAR`;
    $forces=`grep 'FORCES: m' OUTCAR`;
    @forces = split /\n/ , $forces;
    @energy = split /\n/ , $energy;
    $n = @forces;
    print "$n \n";

    open OUT , ">fe.dat" ;

    for($j=0; $j<$n; $j++) {
      $f = $forces[$j]; chomp($f); $f=~s/^\s+//g; @f=split /\s+/,$f;
      $e = $energy[$j]; chomp($e); $e=~s/^\s+//g; @e=split /\s+/,$e;

      if(!$j) { $e0 = $e[6]; }
      printf OUT "%5i %20.8f %20.6f %20.6g \n",$j,$f[4],$e[6],$e[6]-$e0;
    }
    close OUT;

    system "gnuplot $Bin/vef.gnu";
    rename "vaspout.eps" , "vaspout$i.eps";
    system "mv vaspout$i.eps ../vaspgr";

    if($zipped) {
       system "$zip -9 OUTCAR &";
    }
    chdir $dir;
}

