#!/usr/bin/env perl
#;-*- Perl -*-

# Goes though all the image folders in a NEB run and forms a movie file.

use Cwd;
use FindBin qw($Bin);
use lib "$Bin";
use Vasp;

$xdatflag = 0; # turn off xdatcar writing, by default
$xyzflag = 0;  # turn off xyz writing, by default

if(!$ARGV[0]) {
    print "\nUsing POSCARs to generate movie\n\n";
    $filename = "POSCAR";
} else {
    print "\nUsing CONTCARs to generate movie\n\n";
    $filename = "CONTCAR";
}

$dir = cwd;
$zip = $ENV{'VTST_ZIP'};
if($zip eq ''){ $zip = 'gzip'; }
$i = 0;
$string = "00";
while(-d $string) {
    chdir $string;
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
        if($i == 0) { $e0 = $e[6]; }
        if($zipped) { system "$zip -9 OUTCAR &"; }
    }

    $if = 1;
    if(!(-e $filename)) { print "copying\n"; system "cp POSCAR CONTCAR"; $if = 0; }
    system "$Bin/pos2con.pl $filename POSCAR.con > /dev/null";
    system "$Bin/con2xyz.pl POSCAR.con > /dev/null";
#    system "$Bin/pos2con.pl POSCAR.con POSCAR_TMP> /dev/null";
    ($coordinates,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective,$description,$filetype) = read_poscar($filename);
    write_poscar($coordinates,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective,$description,"POSCAR_TMP",$filetype);
    
    unlink "POSCAR.con";
    if(!$if) { print "unlinking\n"; unlink "CONTCAR"; }
    if($outcar) {
        $e = $e[6] - $e0;
        system "sed s/'Generated with con2xyz'/'F:  $f[4]  ...  E:  $e'/g POSCAR.xyz > POSCAR_TMP.xyz";
    } else {
        system "sed s/'Generated with con2xyz'/' NO OUTCAR FOUND '/g POSCAR.xyz > POSCAR_TMP.xyz";
    }
    system "cp POSCAR_TMP.xyz ..//p$i.xyz";
    system "cp POSCAR_TMP ..//pp$i";
    unlink "POSCAR.xyz";
    unlink "POSCAR_TMP";
    rename "POSCAR_TMP.xyz","POSCAR.xyz";

    $i++;
    if($i < 10) { $string = "0$i"; }
    elsif($i < 100) {$string = "$i"; }
    else {die "Too many images"; }
    chdir $dir;
}

if (-e movie){
    unlink "movie"; } 
    if (-e movie.xyz){
        unlink "movie.xyz"; }
        if (-e movie.xdatcar){
            unlink "movie.xdat"; }
            $i--;

            if($i < 10) {
                system "cat pp[0-9] > movie";
                system "cat p[0-9].xyz > movie.xyz";
            } elsif($i < 100) {
                system "cat pp[0-9] pp[0-9][0-9] > movie";
                system "cat p[0-9].xyz p[0-9][0-9].xyz > movie.xyz";
            } else {
                print " TOO MANY IMAGE FILES ...\n";
            }

# make an XDATCAR movie

open (OUT,">movie.XDATCAR");
for($j=0; $j<=$i; $j++){
    $file="pp".$j;
    ($coordinates,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective,$description,$filetype) = read_poscar($file);
    if($j==0){
        system "cp $file movie.POSCAR";
        print OUT $description."\n";
        print OUT "           1"."\n";
        for ($ii=0; $ii<3; $ii++) {
            for ($jj=0; $jj<3; $jj++) {
                printf OUT "%13.6f", $basis->[$jj][$ii]." ";
            }
            print OUT "\n";
        }
# this right now sets default as vasp5 file type 
#if ($filetype eq "vasp5") {
    if ($filetype ne "vasp4") {
#            print OUT "   ".$description."\n";
    }
    for ($ii=0; $ii<@{$num_atoms}; $ii++) {
        printf OUT "%4i", $num_atoms->[$ii];
    }
#        print OUT "\n";
}
print OUT "\n";
for ($ii=0; $ii<$total_atoms; $ii++) {
    print OUT " ";
    for ($jj=0; $jj<3; $jj++) {
        $coord = $coordinates->[$ii][$jj];
        if ($coord>1) { $coord -= 1; }
        elsif ($coord<0) { $coord += 1; }
        printf OUT "%12.8f", $coord."   ";
    }
    print OUT "\n";
}
}
close(OUT);

#  system "rm -f [0-9][0-9]/POSCAR.xyz";

unlink glob "p*.xyz";
unlink glob "pp*";

# remove the xdatcar and xyz movies (we are not currently using them, typically)

if($xyzflag==0){
    unlink "movie.xyz";
}
if($xdatflag==0){
    unlink "movie.POSCAR";
    unlink "movie.XDATCAR";
}

