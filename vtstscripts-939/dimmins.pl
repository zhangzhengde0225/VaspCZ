#!/usr/bin/env perl
#;-*- Perl -*-

#----------------------------------------------------------------------
# This script gets one POSCAR and one MODECAR file and a distance to generate
# two initial images which can be minimized from a converged dimer.
# These two images are "distance" away from the center of the  converged dimer
# along the dimer axis
# written by Lijun Xu and Graeme Henkelman
# Last Modified Sept. 30, 2011 by Lijun Xu, at UVa
#----------------------------------------------------------------------

use FindBin qw($Bin);
use lib "$Bin";
use Vasp;

(@ARGV<5) || die "USAGE: dimmins.pl <POSCAR MODECAR distance replace>\n" ;

# Defaults
$distance = 0.1;
$poscarfilename = "POSCAR";
$modecar = "MODECAR";
$minsfolder = "mins";
$replace_this_poscar = "";

# Read input variables
if(@ARGV>3) { $replace_this_poscar = $ARGV[3]; }
if(@ARGV>2) { $distance = $ARGV[2]; }
if(@ARGV>1) {
    $poscarfilename = $ARGV[0];
    $modecar = $ARGV[1];
}
# print "Displacement along the dimer: $distance\n";

#----------------------------------------------------------------------
#  Read POSCAR and MODECAR files
#----------------------------------------------------------------------

print " Reading POSCAR and MODECAR files\n";
(-e "$poscarfilename") || die (" POSCAR file $poscarfilename for saddle point does not exist\n");
(-e "$modecar") || die (" MODECAR file $modecar for saddle point  does not exist\n");
($coordinates,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective,$description,$filetype)
 = read_poscar($poscarfilename);
($r_vector,$total_atoms2) = read_othercar($modecar);
($total_atoms==$total_atoms2) || die (" POSCAR and MODECAR files should have the same number of atoms\n");
set_bc($coordinates,$total_atoms);
$coordinates = dirkar($coordinates,$basis,$lattice,$total_atoms);

#----------------------------------------------------------------------
# Generate two images by making disploacements
#----------------------------------------------------------------------

print " Generating displacements\n";
$r_mag = magnitude($r_vector,$total_atoms);
$partial = $distance/$r_mag;
$dr1 = vmult($r_vector,$partial,$total_atoms);
$partial = -1.0*$partial;
$dr2 = vmult($r_vector,$partial,$total_atoms);
$image1 = vsum($coordinates,$dr1,$total_atoms);
$image2 = vsum($coordinates,$dr2,$total_atoms);
$image1 = kardir($image1,$basis,$lattice,$total_atoms);
$image2 = kardir($image2,$basis,$lattice,$total_atoms);
set_bc($image1,$total_atoms);
set_bc($image2,$total_atoms);

#----------------------------------------------------------------------
# Write minimization images
#----------------------------------------------------------------------
$tmpfilename="ciPOSCAR";
if(!$replace_this_poscar) {
    print " Writing the dimer image files\n"; 
    if(-e "$minsfolder") { die " Error: $minsfolder already exists.\n"; }
    system "mkdir $minsfolder";
    system "mkdir $minsfolder/min1";
    system "mkdir $minsfolder/min2";
    write_poscar($image1,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective,$description,$tmpfilename,$filetype);
    system "mv $tmpfilename $minsfolder/min1/POSCAR";
    write_poscar($image2,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective,$description,$tmpfilename,$filetype);
    system "mv $tmpfilename $minsfolder/min2/POSCAR";
    system "cp KPOINTS POTCAR $minsfolder/min1/";
    system "cp KPOINTS POTCAR $minsfolder/min2/";
    if(-e "INCAR_min"){
        system "cp INCAR_min $minsfolder/min1/INCAR; cp INCAR_min $minsfolder/min2/INCAR"; 
    }else{
        system "cp INCAR $minsfolder/min1/INCAR; cp INCAR $minsfolder/min2/INCAR";
        print " Warning: copying INCAR from dimer run -- need to change for minimization\n";
    }
    # fix this later.
    if(-e "akmc_min.sub") { system "cp akmc_min.sub $minsfolder/min1/; cp akmc_min.sub $minsfolder/min2/"; }
}else{
    print "replace the poscar in $replace_this_poscar\n"; # by assuming other necessary files are already in the place
    if($replace_this_poscar eq "min1") {
        write_poscar($image1,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective,$description,$tmpfilename,$filetype);
        system "mv $tmpfilename $minsfolder/min1/POSCAR";
    }elsif($replace_this_poscar eq "min2") {
        write_poscar($image2,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective,$description,$tmpfilename,$filetype);
        system "mv $tmpfilename $minsfolder/min2/POSCAR";
    }else{
        print "Warning in dimmins.pl: $replace_this_poscar is not recogniazable. Nothing will happen. Check it!\n";
    }
}

