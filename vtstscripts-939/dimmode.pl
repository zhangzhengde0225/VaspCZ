#!/usr/bin/env perl
#;-*- Perl -*-
#####################################################################################
# Make movies of the negative mode of the dimer                                     #
# Written by Lijun Xu                                                               #
#####################################################################################
use FindBin qw($Bin) ;
use lib "$Bin" ;
use Vasp ;

$poscar="CENTCAR";
if(@ARGV>0){$poscar=$ARGV[0];
  if($poscar eq "--help"){die "Usage: dimcheck.pl CENTCAR <NEWMODECAR> <numimages> <dist>\n";}
}
$MODECAR="NEWMODECAR";
if(@ARGV>1){$MODECAR=$ARGV[1];}
$numimages=32;
if(@ARGV>2){$numimages=$ARGV[2];}
$dist=0.5;
if(@ARGV>3){$dist=$ARGV[3];}

print "<For usage help, type \"dimcheck.pl --help\">\n";
print "Using $poscar, $MODECAR, $numimages images and $dist Angstrom displacement ... \n";

# read in poscar
($coordinates,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective,$description)=read_poscar($poscar);
set_bc($coordinates,$total_atoms);
$coordinates = dirkar($coordinates,$basis,$lattice,$total_atoms);

# read in MODECAR
($mode, $total_atoms_mode) = read_othercar($MODECAR);
if($total_atoms_mode != $total_atoms){die " $MODECAR is inconsistent with $poscar\n";}
# atom type
chomp($temp=`head -1 $poscar`);
$temp=~s/^\s+//g;
$temp=~s/\s+$//g;
@atype=split(/\s+/,$temp);
if(@atype != @$num_atoms){ die "the number of atom types are messed up\n";}
#check the norm of the mode
$norm=0;
for($i=0;$i<@$mode;$i++){
  for($j=0;$j<3;$j++){
    $norm+=$mode->[$i][$j]*$mode->[$i][$j];
  }
}
$norm=sqrt($norm);
if(abs($norm-1)>0.001){
  print "WARNING! Mode ",$m," has norm ",$norm,"\n";
  if($norm == 0){ die "just found a zero mode. Go check the $MODECAR.\n";}
}
$img_mode_movie="dimmode.xyz";
open(XYZFILE, ">$img_mode_movie") || die "cannot open $img_mode_movie\n";
for($n=0;$n<$numimages;$n++){
  # create images cover one period of vibration
  $stepsize=$dist*sin((2.0*3.1415926*$n)/(1.0*$numimages));
  print XYZFILE $total_atoms."\n";
  print XYZFILE "The current imaginery mode\n";
  for($i=0; $i<$total_atoms; $i++){
    for($j=0;$j<3;$j++){
      $image->[$i][$j]=$coordinates->[$i][$j]+$stepsize*$mode->[$i][$j]; 
    }
  }
  # apply periodic boundary conditions (probably we don't need to do this)
  $image=kardir($image,$basis,$lattice,$total_atoms);
  set_bc($image,$total_atoms);
  $image=dirkar($image,$basis,$lattice,$total_atoms);
  # write to the moviefile    
  $k=-1;
  for($x=0;$x<@$num_atoms;$x++) {
    for($y=0;$y<$num_atoms->[$x];$y++){
       $k++;
       print XYZFILE $atype[$x]."  ".$image->[$k][0]."  ".$image->[$k][1]."  ".$image->[$k][2]."\n";
    }
  }
  $k++;
  if($k != $total_atoms){ die "check POSCAR file: make sure the number of atoms is right\n";}
}
close XYZFILE;
print "mode movie: $img_mode_movie\n";
