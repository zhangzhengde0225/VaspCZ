#!/usr/bin/env perl
#;-*- Perl -*-
#####################################################################################
# Make movies of the vibrational modes calculated from the dynamical matrix         #
# written by Lijun Xu and Graeme Henkelman in the Univeristy of Texas at Austin     #
# Last modified by GH on May 23, 2008                                               #
#####################################################################################
use FindBin qw($Bin) ;
use lib "$Bin" ;
use Vasp ;

#comment
#@ARGV>0 || die "Usage: dymmodes2xyz.pl POSCAR <DISPLACECAR> <modes.dat> <moviedir> <freq.dat> <numimages> <dist>\n";
$poscar="CONTCAR";
if(@ARGV>0){$poscar=$ARGV[0];}
$displacecar="DISPLACECAR";
if(@ARGV>1){$displacecar=$ARGV[1];}
$modesfile="modes.dat";
if(@ARGV>2){$modesfile=$ARGV[2];}
$moviefolder=".";
#if(@ARGV>3){$moviefolder=$ARGV[3];}
$freqfile="freq.dat";
#if(@ARGV>4){$moviefolder=$ARGV[4];}
$numimages=32;
#if(@ARGV>5){$moviefolder=$ARGV[5];}
$dist=0.5;
#if(@ARGV>6){$moviefolder=$ARGV[6];}
if(@ARGV>3){ $moviefolder=$ARGV[3]; if(!(-e "$moviefolder")){system "mkdir $moviefolder";} }
if(@ARGV>4){$freqfile=$ARGV[4];}
if(@ARGV>5){$numimages=$ARGV[5];}
if(@ARGV>6){$dist=$ARGV[6];}

# read in poscar
($coordinates,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective,$description)=read_poscar($poscar);
set_bc($coordinates,$total_atoms);
$coordinates = dirkar($coordinates,$basis,$lattice,$total_atoms);

# figure out the order of the displaced atoms
($displacement,$link)=handle_displacecars($displacecar,$total_atoms);
#($displacement,$totalatoms_disp)=read_othercar($displacecar);
#if($total_atoms != $totalatoms_disp){ die "$poscar and $displacecar should have the same number of atoms.\n";}

# figure out the atom types
@atype=split(/\s+/,$description);
if(@atype < @$num_atoms){
  for($i=@atype;$i<@$num_atoms;$i++){
    push @atype, "Type";
  }
}
print "@atype\n";

# read in modes and frequencies
($modes,$numDOF)=readmodes($modesfile);
if($numDOF !=@$displacement){die "Displacecar $displacecar is inconsistent with modes $modesfile\n";}
$freq=readfreq($freqfile,$numDOF);
#for($i=0;$i<$numDOF;$i++){
#   print "displacement ".$displacement->[$i]." linked with "."@{$link->[$i]}\n";
#}

# for each mode, create a movie of $numimages frames 
for($m=0;$m<@$modes;$m++){
  #check the norm of the mode
  $norm=0;
  for $item (@{$modes->[$m]}){ $norm+=$item*$item; }
  $norm=sqrt($norm);
  if(abs($norm-1)>0.001){
    print "Mode ",$m," has norm ",$norm,"\n";
    if($norm == 0){ die "just found a zero mode. Go check the modefile $modesfile.\n";}
  }

  # figure out the file names
  if($m<9){$modename=$moviefolder."/mode00".($m+1).".xyz";}
  elsif($m<99){$modename=$moviefolder."/mode0".($m+1).".xyz";}
  else{$modename=$moviefolder."/mode".($m+1).".xyz";}

  # add the mode to the equilibrium geometry
  open(XYZFILE, ">$modename") || die "cannot open $modename\n";
  for($n=0;$n<$numimages;$n++){
    # create images cover one period of vibration
    $stepsize=$dist*sin((2.0*3.1415926*$n)/(1.0*$numimages));
    print XYZFILE $total_atoms."\n";
    print XYZFILE $freq->[$m]."\n";
    for($i=0; $i<$total_atoms; $i++){
      for($j=0;$j<3;$j++){
        $image->[$i][$j]=$coordinates->[$i][$j];          
      }
    } 
    for($i=0;$i<$numDOF;$i++){
       $atom=$link->[$i][0];
       $j=$link->[$i][1];
       $image->[$atom][$j]=$image->[$atom][$j]+$stepsize*$modes->[$m][$i];
    }

    # apply periodic boundary conditions (probably we don't need to do this)
    $image=kardir($image,$basis,$lattice,$total_atoms);
    set_bc($image,$total_atoms);
    $image=dirkar($image,$basis,$lattice,$total_atoms);

    # write to the moviefile    
    $k=-1;
    for($i=0;$i<@$num_atoms;$i++) {
      for($j=0;$j<$num_atoms->[$i];$j++){
         $k++;
         print XYZFILE $atype[$i]."  ".$image->[$k][0]."  ".$image->[$k][1]."  ".$image->[$k][2]."\n";
      }
    }
    $k++;
    if($k != $total_atoms){ die "check POSCAR file: make sure the number of atoms is right\n";}
  }
  close XYZFILE;
  print "Done with mode # ".($m+1)."\n";
}

# -------------------------------------
# ($displacement,$link)=handle_displacecars($displacecar,$total_atoms)
# -------------------------------------
sub handle_displacecars{
  my ($displacecar,$total_atoms)=@_;
  my ($displacement,$link,$i,$j,$k,$atom);
  my ($temp,$natoms_chunk);
  ($temp,$natoms_chunk)=read_othercar($displacecar);
  if($natoms_chunk < $total_atoms){ die "the number of atoms in the $displacecar should be same as $total_atoms";}
  $k=-1;
  $atom=-1;
  for($i=0;$i<$natoms_chunk;$i++){
    $atom++;
    if($atom==$total_atoms){ $atom=0;}
    for($j=0;$j<3;$j++){
      if($temp->[$i][$j] > 0){
        $k++;
        $displacement->[$k]=$temp->[$i][$j];
        $link->[$k]=[($atom,$j)];
      }
    }
  }
  if($atom==$total_atoms-1){ return ($displacement,$link);
  }else{ die "The chunky $displcarcar should have integer times $total_atoms atoms.\n";}
}

# -------------------------------------
# ($modes,$numDOF)=readmodes($modesfile);
# -------------------------------------
sub readmodes{ # column vectors
  my $filename=shift;
  my ($m,$n,$line,@lines,$modes,$numDOF);
  open(IN, "<$filename") || die "cannot open $filename.\n";
  $n=0;
  $numDOF=0;
  while($line=<IN>){
    $line=~s/^\s+//;
    if($line eq "\n"){next;}
    chomp($line);
    @lines=split(/\s+/,$line);
    for($m=0;$m<@lines;$m++){
       $modes->[$m][$n]=$lines[$m];
    }
    if(!$numDOF){
      $numDOF=@lines;
    }elsif($numDOF!=@lines){
      die "In modefile $filename, each mode should have the same number of degrees of freedom.\n";
    }
    $n++;
  }
  close IN;
  if($n != $numDOF){ die "Warning: the number of modes is not equal to the number of degrees of freedom.\n";}
  return ($modes,$numDOF);
}

# -------------------------------------
# $freq=readfreq($freqfile);
# -------------------------------------
sub readfreq{
  my ($filename,$numDOF)=@_;
  my ($freq,$i);
  if($filename eq "") {
    for($i=0;$i<$numDOF;$i++){
      $freq->[$i]="Generated by dymmodes2xyz.pl";
    }
  }else{
    open(IN, "<$filename") || die "cannot open $filename.\n";
    $i=0;
    while($line=<IN>){
      $line=~s/^\s+//;
      if($line eq "\n"){next;}
      chomp($freq->[$i]=$line);
      $i++;
    }
    if($i != $numDOF){die "$filename and modefile should have the number of modes\n";}
    close IN;
  }
  return $freq;
}
