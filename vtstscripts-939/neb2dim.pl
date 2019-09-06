#!/usr/bin/env perl
#;-*- Perl -*-

# Set up a dimer saddle point search from a NEB run.
# Assume for now that the run as been stopped and vfin.pl run so the
# configuration is contain in POSCARs and not CONTCAR. 
# If there is an input argument then it the image number to be used.

use FindBin qw($Bin);
use lib "$Bin";

if($ARGV[0]) {
    $i = $ARGV[0];
    $f = 0;
} else {
    # Select the image from the exts.dat file
    open EXT , "exts.dat"
        or die " NO EXTS.DAT IN THIS DIRECTORY; RUN nebresults.pl \n";
    @ext = <EXT>;
    close EXT;
    chomp(@ext);
    $n = @ext;
    $max = -1.0;
    for($i=0; $i<$n; $i++) {
        $t = $ext[$i];
        $t =~ s/^\s+//g;
        @t = split /\s+/,$t;
        if($t[8] > $max) {
            $max = $t[8];
            $im = $t[5];
        }
    }
    $i = int($im);
    $f = $im - $i;
    if($im<=0) {
        $i=0;
        $f=0;
    }
}
print "\n FORMING DIMER BETWEEN IMAGES ",$i," and ",$i+1,"\n\n";

die "NUMBER OF IMAGES MUST BE LESS THAN 99 FOR NOW ! \n" if $i >= 99;

$d0= sprintf "%02d",$i-1;
$d1 = sprintf "%02d",$i;
$d2 = sprintf "%02d",$i+1;
$d3 = sprintf "%02d",$i+2;

# Checks to see whether i-1 and i+2 images exist

if (-e "$d3/POSCAR") {
    $temp1=$d3;
} else {
    $temp1=$d2;
}

if (-e "$d0/POSCAR") {
    $temp2=$d0;
} else {
    $temp2=$d1;
}

$d3=$temp1;
$d0=$temp2;
print $d3, "\n";

open (DIST, "$Bin/dist.pl $d1/POSCAR $d2/POSCAR |");
$d12 = <DIST>;
close DIST;
$d12>0 || die "ZERO DISTANCE BETWEEN IMAGES ",$i1," AND ",$i2,"\n";
open (DIST, "$Bin/dist.pl $d0/POSCAR $d2/POSCAR |");
$d02 = <DIST>;
close DIST;
open (DIST, "$Bin/dist.pl $d1/POSCAR $d3/POSCAR |");
$d13 = <DIST>;
close DIST;


#$NdR = 5e-3;
#$f1 = $f - $NdR/$d12;
#$f2 = $f + $NdR/$d12;

# Generate the dimer POSCARs.

#system "$Bin/posinterp.pl $d1/POSCAR $d2/POSCAR $f1 > /dev/null";
#rename "POSCAR.out" , "POSCAR1";
#system "$Bin/posinterp.pl $d1/POSCAR $d2/POSCAR $f2 > /dev/null";
#rename "POSCAR.out" , "POSCAR2";
system "$Bin/posinterp.pl $d1/POSCAR $d2/POSCAR $f > /dev/null";
rename "POSCAR.out" , "POSCAR";
system "cp $d0/POSCAR POSCAR0";
system "cp $d1/POSCAR POSCAR1";
system "cp $d2/POSCAR POSCAR2";
system "cp $d3/POSCAR POSCAR3";

system "$Bin/modemake.pl POSCAR0 POSCAR2 ; mv MODECAR MODECAR02";
system "$Bin/modemake.pl POSCAR1 POSCAR3 ; mv MODECAR MODECAR13";
$natoms= `wc -l< MODECAR02`;

$p1 = "MODECAR02";
$p2 = "MODECAR13";


# Read the  MODECAR02 file
  open P1 , $p1 or die "$p1 not found in this directory\n" ;
  while (<P1>){$file .= $_ ;}
  @p1 = split /\n/ , $file ;
  $n1 = @p1 ;
  close P1 ;

# Read the MODECAR13 file
  open P2 , $p2 or die "$p2 not found in this directory\n" ;
  $file = undef ;
  while (<P2>){$file .= $_ ;}
  @p2 = split /\n/ , $file ;
  $n2 = @p2 ;
  close P2 ;

  for ($i=0; $i<$natoms; $i++){
    $l1 = $p1[$i] ; chomp($l1) ; $l1=~s/^\s+//g ;
    @l1=split /\s+/,$l1 ;
    # print "$l1 \n";
    $l2 = $p2[$i] ; chomp($l2) ; $l2=~s/^\s+//g ; @l2=split /\s+/,$l2 ;
    
    
    $m[$i][0] =$l1[0]*(1-$f)+$l2[0]*$f;
    $m[$i][1] =$l1[1]*(1-$f)+$l2[1]*$f;
    $m[$i][2] =$l1[2]*(1-$f)+$l2[2]*$f;
}

  for ($i=0; $i<$natoms; $i++){
    for ($j=0; $j<3; $j++){
      $sum += $m[$i][$j] * $m[$i][$j] ;
    }
  }
  $sum = sqrt($sum) ;
#  print "sum is $sum \n";
  open MOD , ">MODECAR" ;
  for ($i=0; $i<$natoms; $i++){
    for ($j=0; $j<3; $j++){
      $m[$i][$j] = $m[$i][$j]/$sum ;
    }
    printf MOD "%20.10E %20.10E %20.10E \n",$m[$i][0],$m[$i][1],$m[$i][2] ;
  }
  close(MOD) ;

$dir = "dim";
mkdir $dir ; # mkdir "$dir/01" ; mkdir "$dir/02";
system "cp INCAR KPOINTS $dir";
if(-e "POTCAR") {
    system "cp POTCAR $dir";
} elsif(-e "../POTCAR") {
    system "cp ../POTCAR $dir";
}
system " mv MODECAR $dir ; mv POSCAR $dir";
#system "$Bin/modemake.pl POSCAR1 POSCAR2 ; mv MODECAR $dir ; mv POSCAR $dir";
#  system "cp POSCAR1 $dir/01/POSCAR" ;
#  system "cp POSCAR2 $dir/02/POSCAR" ;
print " FOR DIMER, REMEMBER TO SET IN THE INCAR: \n\n";
print "ICHAIN = 2 \n";
print "IOPT = 2 \n";
print "IBRION = 3 \n";
print "POTIM  = 0.0 \n";
print "EDIFF  = 1E-7 \n \n";
print "DdR       (DEFAULT = 5E-3) \n";
print "DRotMax   (DEFAULT = 4) \n";
print "DFNMax    (DEFAULT = 1.0) \n";
print "DFNMin    (DEFAULT = 0.01) \n\n";

unlink glob "POSCAR1";
unlink glob "POSCAR2";
unlink glob "POSCAR3";
unlink glob "POSCAR4";
unlink glob "MODECAR02";
unlink glob "MODECAR13";
