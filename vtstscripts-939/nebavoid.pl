#!/usr/bin/env perl
#;-*- Perl -*-
use FindBin qw($Bin);
use lib "$Bin";
use Vasp;

################################################
# Generating neb images by interpolation can create images in which two (or more) atoms 
# are too close to each other. This script applies a repelling force between each pair
# and pushes them apart until the minimum distance is reached. 
# However, this might disturb the equal space between neighboring images. We may fix it later.
# Written by Lijun Xu and Graeme Henkelman in the Univeristy of Texas at Austin on March 6, 2007
# Last modified by Lijun on March 8, 2007
################################################

if(@ARGV>0) {
    $rcut = $ARGV[0];
    print "The minimum distance between two atoms is $rcut Angstrom.\n";
} else {
    die "Usage: nebavoid.pl distance (the minimum distance between two atoms; recommended value is 0.1 Angstrom.\n";
}
$ORTHOGONAL = 0;
if(@ARGV>1) { $ORTHOGONAL = $ARGV[1]; }

chomp($folder=`ls -d [0-9][0-9]`);
$folder =~ s/^\s+//;
@folders = split(/\s+/,$folder);
$num = @folders-1;
for($k=1; $k<$num; $k++) {
    $targetfolder = $folders[$k];
    $dummy = substr($targetfolder, -1 , 1);
    if($dummy eq "\/") { chop($targetfolder); }
    if(!(-e $targetfolder)){ die "folder $targetfolder does not exist.\n"; }
    $POSCAR = $targetfolder."/POSCAR";
    if(!(-e $POSCAR)) { die "file $POSCAR does not exist.\n"; }
    print "working on $targetfolder...\n";
    ($coordinates,$basis,$lattice,$num_atom,$total_atoms,$selectiveflag,$selective,$description) = read_poscar($POSCAR);
    set_bc($coordinates,$total_atoms);
    for($i=0; $i<3; $i++){
        for($j=0; $j<3; $j++){
            $lattice_vec->[$i][$j] = $basis->[$j][$i];
        }
    }
    #----------------------------------
    #Some global variables for orthogonal boxes
    if($ORTHOGONAL) {
        @BOX = (abs($lattice_vec->[0][0]),abs($lattice_vec->[1][1]),abs($lattice_vec->[2][2]));
        @BOX_ = (-1.0*$BOX[0], -1.0*$BOX[1], -1.0*$BOX[2]);
        @HALFBOX = ($BOX[0]*0.5, $BOX[1]*0.5, $BOX[2]*0.5);
        @HALFBOX_ = (-0.5*$BOX[0], -0.5*$BOX[1], -0.5*$BOX[2]);
    }
    $coordinates = dirkar($coordinates,$basis,$lattice,$total_atoms);
    spring_relaxation($coordinates,$basis,$lattice,$lattice_vec,$total_atoms,$selective,$rcut);
    $coordinates = kardir($coordinates,$basis,$lattice,$total_atoms);
    system "cd $targetfolder; mv POSCAR POSCAR_orig";
    write_poscar($coordinates,$basis,$lattice,$num_atom,$total_atoms,$selectiveflag,$selective,$description,$POSCAR);
}

sub spring_relaxation{
    my ($R,$basis,$lattice,$lattice_vec,$totalatoms,$selective,$cutoff) = @_;
    my $epsilon = 0.1;
    my $stepmax = 100000;
    my $drmax = 0.2;
    my @line = ();
    my ($i,$j,$k,$difference,$Rij);
    my ($force,$fij,$phiij,$step,$goodsound,$converged,$frozen,$dr,$rcut);
    $rcut = $cutoff*$cutoff;
    # figure out frozen flags for each degree of freedom
    for($i=0; $i<$totalatoms; $i++) {
        @line = split(/\s+/,$selective->[$i]);
        for($j=0; $j<3; $j++) {
            $frozen->[$i][$j] = $line[$j];
        }
    }
    # start the steepest descent loop
    $converged = 0;
    for($step=0; $step<$stepmax; $step++) {
        if($converged) { last; }
        print "Step No. $step\n";
        for($i=0; $i<$totalatoms; $i++) {
            for($j=0; $j<3; $j++) {
                $force->[$i][$j] = 0.0;
            }
        }
        $goodsound = 1;
        for($i=0; $i<$totalatoms-1; $i++) {
            for($j=$i+1; $j<$totalatoms; $j++) {
                for($k=0; $k<3; $k++){
                    $difference->[$k] = $R->[$i][$k]-$R->[$j][$k];
                }
                ($difference, $Rij) = FindMinimumImage($difference,$lattice_vec);
                if($Rij < $rcut) {
                    #print "($i, $j)"."Rij=$Rij"." rcut=$rcut\n";
                    $goodsound = 0;
                    $Rij = sqrt($Rij);
                    ($fij,$phiij) = linear_repulsion_pot($Rij,$cutoff);
                    #($fij,$phiij)=exp_repulsion_pot($Rij,$cutoff);
                    $fij = $fij/$Rij;
                    for($k=0; $k<3; $k++){
                        $phiij = $fij*$difference->[$k];
                        $force->[$i][$k] -= $phiij;
                        $force->[$j][$k] += $phiij;
                    }
                    #print "@{force->[$i]}\n";
                }
            }
        }
        if($goodsound) {
            $converged = 1;
            next;
        }
        for($i=0; $i<$totalatoms; $i++) {
            for($j=0; $j<3; $j++) {
                if(lc($frozen->[$i][$j]) eq "f") {
                    next;
                }
                $dr = $epsilon*$force->[$i][$j];
                #print "$i: $j -- dr: $dr dmax:$drmax\n"; 
                if(abs($dr) > $drmax) {
                    $dr=$drmax*$dr/abs($dr);
                }
                $R->[$i][$j] += $dr;
            }
        }
        if($ORTHOGONAL) {
            for($i=0; $i<$totalatoms; $i++){
                for($j=0; $j<3; $j++){
                    while($R->[$i][$j] > $BOX[$j]) { $R->[$i][$j] -= $BOX[$j]; }
                    while($R->[$i][$j] < $BOX_[$j]) { $R->[$i][$j] += $BOX[$j]; }
                }
            }
        } else {
            $R = kardir($R,$basis,$lattice,$total_atoms);
            set_bc($R,$total_atoms);
            $R = dirkar($R,$basis,$lattice,$total_atoms);
        }
    }
    if($converged == 1) {
        print "fully relaxed after $step (1: trivial) steps\n";
    } else {
        print "Not fully relaxed after $stepmax steps, but we can;t wait. move on\n";
    }
}

sub linear_repulsion_pot{ # E=kx-b, k< 0
    my ($r,$cutoff) = @_;
    my $slope = -0.3;
    my ($energy, $fij);
    my $intercept = $slope*$cutoff;
    $energy = $slope*$r-$intercept;
    $fij = $slope;
    return ($fij,$energy);
} 

#----------------------------------------------------------------------
# FindMinimumImage: (INPUT) A vector (pointer) (OUTPT) the minimum image vector
#----------------------------------------------------------------------
sub FindMinimumImage{
    my ($dr_car,$lattice_vec) = @_;
    my ($dsqmin,$done,$v1,$v2,$v3,$drt_car,$dsq);
    my ($i,$j,$k);
    if($ORTHOGONAL){
        for($i=0; $i<@$dr_car; $i++){
            if($dr_car->[$i]<$HALFBOX_[$i]) { $dr_car->[$i] += $BOX[$i]; }
            if($dr_car->[$i]>$HALFBOX[$i]) { $dr_car->[$i] -= $BOX[$i]; }
        }
        $dsqmin = VecDot($dr_car,$dr_car);
    } else {
        $dsqmin = VecDot($dr_car,$dr_car);
        $done = 1;
        while($done) {
            $found = 1;
            for($i=-1; $i<2; $i++) {
                $v1 = VecMultiply($lattice_vec->[0],$i);
                for($j=-1; $j<2; $j++) {
                    $v2 = VecMultiply($lattice_vec->[1],$j);
                    for($k=-1; $k<2; $k++) {
                        $v3 = VecMultiply($lattice_vec->[2],$k);
                        $drt_car = VecSum($dr_car,$v1,$v2,$v3);
                        $dsq = VecDot($drt_car,$drt_car);
                        if($dsq<$dsqmin) {
                            $dr_car = $drt_car;
                            $dsqmin = $dsq;
                            $found = 0;
                        }
                    }
                }
            }
            if($found){ $done=0; }
        }
    }
    #print sqrt($dsqmin)."\n";
    return ($dr_car,$dsqmin);
}

#----------------------------------------------------------------------
# VecDot: (INPUT) two vectors (pointers) (OUTPUT) a dot product
#----------------------------------------------------------------------
sub VecDot{
    my ($vec1, $vec2) = @_;
    my ($i, $dot);
    #check dimensions before calling this sub
    $dot = 0.0;
    for($i=0; $i<@$vec1; $i++){
        $dot += $vec1->[$i]*$vec2->[$i];
    }
    return $dot;
}

#----------------------------------------------------------------------
# VecSum: (INPUT) 3-d vectors (pointers) (OUTPUT) a sum vector (pointer)
#----------------------------------------------------------------------
sub VecSum{
    my ($dimension, $i, $j, @sum);
    $dimension = @_;
    @sum = (0.0,0.0,0.0);
    for($i=0; $i<$dimension; $i++){
        for($j=0; $j<3; $j++){
            $sum[$j] += $_[$i]->[$j];
        }
    }
    return \@sum;
}

#----------------------------------------------------------------------
# VecMultiply: (INPUT) one vector (pointer) (OUTPUT) a new vector (pointer)
#----------------------------------------------------------------------
sub VecMultiply{
    my ($vec, $fac) = @_;
    my ($i, $vecout);
    for($i=0; $i<@$vec; $i++){
        $vecout->[$i] = $vec->[$i]*$fac;
    }
    return $vecout;
}


