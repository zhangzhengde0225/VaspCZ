#!/usr/bin/env perl
#;-*- Perl -*-

###########################################################################
#                                                                         #
# Goes though all the image folders in a NEB                              #
# run and forms a vasp-file from the POSCAR.                              #
# Note that the first line in the POSCAR has                              #
# to have the elemental symbols in the same                               #
# order as in the POTCAR, i.e. for a run the                              #
# methyl and a Ni slab                                                    #
#                                                                         #
# Ni C H                                                                  #
#  4.97615083550000037                                                    #
#      1.0000000000000000    0.0000000000000000    0.0000000000000000     #
#      0.0000000000000000    0.8660254037832008    0.0000000000000000     #
#      0.0000000000000000    0.0000000000000000    3.6172536956853296     #
#   16   1   4                                                            #
# Selective dynamics                                                      #
# Direct                                                                  #
#                                                                         #
###########################################################################

use Cwd;
use FindBin qw($Bin);
use lib "$Bin";

if(!$ARGV[0]) {
    print "\nUsing POSCARs to generate movie\n\n";
    $filetype = "POSCAR";
} else {
    print "\nUsing CONTCARs to generate movie\n\n";
    $filetype = "CONTCAR";
}

$dir = cwd;
$zip = $ENV{'VTST_ZIP'};
if($zip eq '') { $zip = 'gzip'; }
$i = 0;
$string = "00";
while(chdir $string) {
    # Grab the forces and energies to add to the xyz files
    $zipped = 0;
    $outcar = 1;
    if(-e "OUTCAR") { ; }
    elsif(-e "OUTCAR.gz"){ $zipped = 1; system "gunzip OUTCAR.gz"; }
    elsif(-e "OUTCAR.bz2"){ $zipped = 1; system "bunzip2 OUTCAR.bz2"; }
    else{ $outcar = 0; }

    if($outcar){
        $for = `grep 'Forces: m' OUTCAR | tail -1`;
        if($for == undef){ $for = `grep 'FORCES: m' OUTCAR | tail -1`; }
        $ene = `grep 'energy  w' OUTCAR | tail -1`;
        $f = $for; chomp($f); $f=~s/^\s+//g; @f=split /\s+/,$f;
        $e = $ene; chomp($e); $e=~s/^\s+//g; @e=split /\s+/,$e;
        if($i == 0){ $e0 = $e[6]; }
        if($zipped){ system "$zip -9 OUTCAR &"; }
    }

    $if = 1;
    if(!(-e $filetype)){print "copying\n"; system "cp POSCAR CONTCAR"; $if = 0;}
    system "cp $filetype TMP > /dev/null";
    system "$Bin/pos2jvasp.pl TMP > /dev/null";
    unlink "TMP";
    if(!$if){ print "unlinking\n"; unlink "CONTCAR"; }
    if($outcar) {
        $e = $e[6] - $e0;
        system "sed s/'SOME INFORMATION WE DONT CARE ABOUT'/'F:$f[4]...E:$e'/g TMP.vasp > POSCAR_TMP.vasp";
    } else {
        system "sed s/'SOME INFORMATION WE DONT CARE ABOUT'/'NO_OUTCAR_FOUND '/g TMP.vasp > POSCAR_TMP.vasp";
    }
    unlink "TMP.vasp";
    system "sed '1,3d'< POSCAR_TMP.vasp >..//p$i.vasp";
    unlink "POSCAR.vasp";
    rename "POSCAR_TMP.vasp" , "POSCAR.vasp";

    $i++;
    if($i < 10) { $string = "0$i"; } 
    elsif($i < 100) { $string = "$i"; }
    else { die "Too many images"; }
    chdir $dir 
}
  
if (-e movie.vasp) { unlink "movie.vasp"; }
system "head -3 00/POSCAR.vasp > movie.vasp";
$i--;

if($i < 10) {
    system "cat p[0-9].vasp >> movie.vasp";
} elsif($i < 100) {
    system "cat p[0-9].vasp p[0-9][0-9].vasp >> movie.vasp";
} else {
    print " TOO MANY vasp FILES ...\n";
}

unlink glob "p*.vasp";

