#!/usr/bin/env perl
#;-*- Perl -*-

# Set up a lanczos saddle point search from n NEB run.
# Assume for now that the run as been stopped and vfin.pl run so the
# configuration is contain in POSCARs and not CONTCAR. 
# If there is an input argument then it is the image number to be used.

use FindBin qw($Bin);
use lib "$Bin";

if($ARGV[0]) {
    $image = $ARGV[0];
    $f = 0;
} else {
# Select the image from the exts.dat file
    open EXT , "exts.dat"
        or die " NO EXTS.DAT IN THIS DIRECTORY \n";
    @ext = <EXT>;
    close EXT;
    chomp(@ext);
    $n = @ext;
    $max = -1.0;
    for($i=0; $i<$n; $i++) {
        $t = $ext[$i];
        $t=~s/^\s+//g;
        @t=split /\s+/,$t;
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
print "\n LANCZOS CREATED BETWEEN IMAGES ",$i," and ",$i+1,"\n\n";

die "NUMBER OF IMAGES MUST BE LESS THAN 99.\n" if $i >= 99;

$d1 = sprintf "%02d",$i;
$d2 = sprintf "%02d",$i+1;
open (DIST, "$Bin/dist.pl $d1/POSCAR $d2/POSCAR |");
$d12 = <DIST>;
close DIST;
$d12>0 || die "ZERO DISTANCE BETWEEN IMAGES ",$i1," AND ",$i2,"\n";
$NdR = 5e-3;
$f1 = $f - $NdR/$d12;
$f2 = $f + $NdR/$d12;

# Generate the dimer POSCARs.

system "$Bin/posinterp.pl $d1/POSCAR $d2/POSCAR $f1 > /dev/null";
rename "POSCAR.out" , "POSCAR1";
system "$Bin/posinterp.pl $d1/POSCAR $d2/POSCAR $f2 > /dev/null";
rename "POSCAR.out" , "POSCAR2";
system "$Bin/posinterp.pl $d1/POSCAR $d2/POSCAR $f > /dev/null";
rename "POSCAR.out" , "POSCAR";

$dir = "lan";
mkdir $dir;
system "cp INCAR KPOINTS $dir";
if(-e "POTCAR") {
    system "cp POTCAR $dir" ;
} elsif(-e "../POTCAR") {
    system "cp ../POTCAR $dir";
}
system "$Bin/modemake.pl POSCAR1 POSCAR2; mv MODECAR $dir; mv POSCAR $dir";
print " FOR LANCZOS, REMEMBER TO SET IN THE INCAR: \n\n";
print "ICHAIN = 3 \n";
print "IBRION = 3 \n";
print "POTIM  = 0.0 \n";
print "EDIFF  = 1E-8 \n \n";
print "SIFCG     (DEFAULT = .TRUE.) \n";
print "SMAXMOVE  (DEFAULT = 0.3) \n";
print "STOL      (DEFAULT = 1E-2) \n";
print "SDR       (DEFAULT = 1E-3) \n";
print "SDT       (DEFAULT = 0.15.) \n";
print "SNL       (DEFAULT = 20) \n\n";

unlink glob "POSCAR1" ;
unlink glob "POSCAR2" ;

