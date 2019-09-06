#!/usr/bin/env perl
#;-*- Perl -*-

#----------------------------------------------------------------------
# This script gets a POSCAR and a DISPLACECAR_sp to generate single-image 
# dimers for vasp dimer calculations.
#
# Last Modified Sept. 30, 2011 by Lijun Xu and Graeme Henkelman, UT-Austin
# Thank Yafan Zhao for spotting problems in applying periodic boundary conditions
# to non-orthogonal boxes.
#----------------------------------------------------------------------

use FindBin qw($Bin);
use lib "$Bin";
use Vasp;

(@ARGV>0) || die " USAGE: diminit.pl DIRECTORY OR NUMBER (OPTIONAL: DisplaceAlgo DisplaceRange Rcut MaxCord POSCAR DISPLACECAR_sp ORTHOGONAL)\n";

$ndimers = 0;
$DisplaceAlgo = 0; #0 lowest coordianted 1 less than a specified value
$DisplaceRange = 0.01; # default: displace itself only
$rcut = 0.01; # default: no neighbors, just itself
$MaxCord = 8;
$poscarfilename = "POSCAR";
$displacecarfilename = "DISPLACECAR_sp";
$ORTHOGONAL = 0; # Captialized variables are global
# The distance between dimer images is now controlled by INCAR
# $deltaR=0.005 ; #0.01 is the distance between two images

if($ARGV[0] =~ /^\d+$/){
    # It's a number, so we create that many dimers in the folders: pr000[1-n]
    $ndimers = $ARGV[0]; }
if(@ARGV>1){ $DisplaceAlgo = $ARGV[1]; }
if(@ARGV>2){ $DisplaceRange = $ARGV[2]; }
if(@ARGV>3){ $rcut = $ARGV[3]; }
if(@ARGV>4){ $MaxCord = $ARGV[4]; }
if(@ARGV>5){ $poscarfilename = $ARGV[5]; }
if(@ARGV>6){ $displacecarfilename = $ARGV[6]; }
if(@ARGV>6){ $ORTHOGONAL = $ARGV[7]; }

#----------------------------------------------------------------------
#  Read in DISPLACEMENT file to decide which atom should be displaced
#  DISPALCEMENT has the same file structure as DISPLACECAR
#  Example:
#    0   0   0   1  : atom 1 will not be displaced
#    0.1 0.1 0.1 2  : atom 2, x, y, z will be displaced by a
#                     0.1-width Gaussian random number
#----------------------------------------------------------------------

open(DISPLACEMENT, "<$displacecarfilename")
    or die (" Error: cannot open the displacement file: ".$displacecarfilename."\n Note: default is DISPLACECAR_sp\n");
close(DISPLACEMENT);
($sigma, $total_atoms_disp) = read_othercar($displacecarfilename);
$DisplaceList=();
for ($i=0; $i<$total_atoms_disp; $i++) {
    #print $sigma->[$i][0]."  ".$sigma->[$i][1]."  ".$sigma->[$i][2]."  ".($i+1)."\n"; 
    if($sigma->[$i][0] != 0 || $sigma->[$i][1] != 0 || $sigma->[$i][2] != 0) {
        push @$DisplaceList, $i;
    }
}
if(@$DisplaceList == 0){ die "Error from diminit.pl; no atoms will be displaced. what is going on here. Check it.\n"; }
#----------------------------------------------------------------------
#  Read POSCAR files and prepare  the random displacement 
#----------------------------------------------------------------------
print " Reading POSCAR\n";
($coordinates,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective,$description,$filetype)=read_poscar($poscarfilename);
($total_atoms_disp == $total_atoms) ||
    die " Error: Displacement file $displacecarfilename should have the same number of lines as there are atoms in the POSCAR file.\n";
set_bc($coordinates,$total_atoms);#important
# build a list of atoms that be a center of displacement(default: displace all allowed to move) and their neighbors
$coordinates = dirkar($coordinates,$basis,$lattice,$total_atoms);
for($i=0; $i<3; $i++){
    for($j=0; $j<3; $j++){
        $lattice_vec->[$i][$j] = $basis->[$j][$i];
    }
}
#----------------------------------
#Some global variables for ortogonal boxes
if($ORTHOGONAL){
    @BOX = (abs($lattice_vec->[0][0]),abs($lattice_vec->[1][1]),abs($lattice_vec->[2][2]));
    @BOX_ = (-1.0*$BOX[0], -1.0*$BOX[1], -1.0*$BOX[2]);
    @HALFBOX = ($BOX[0]*0.5, $BOX[1]*0.5, $BOX[2]*0.5);
    @HALFBOX_ = (-0.5*$BOX[0], -0.5*$BOX[1], -0.5*$BOX[2]);
}
#---------------------------------
($NcList,$NN_DisplaceList)=BuildDisplaceList($coordinates,$lattice_vec,$DisplaceList,$DisplaceAlgo,$DisplaceRange,$rcut,$MaxCord);
print "**********************************\n";
print "@$NcList\n";
print "**********************************\n";
for $i (@$DisplaceList){
   print "$i: @{$NN_DisplaceList->{$i}}\n";
}
#----------------------------------------------------------------------
#  Do the displacement to generate dimer(s)
#----------------------------------------------------------------------
if($ndimers == 0){
    $NumSearches = 1;
    $folder = $ARGV[0];
}else{
    $NumSearches = $ndimers;
}
print " Generating files for $NumSearches dimer searches\n";
srand();
for($k=1; $k<=$NumSearches; $k++){
    $randomnum = rand();
    print "randomnum=$randomnum\n";
    $atoms2bdisplaced = BuildNewDisplacecar($DisplaceList,$NcList,$NN_DisplaceList,$DisplaceAlgo,$randomnum);
    print "&&&@$atoms2bdisplaced\n";
    # folder name
    if($ndimers > 0) {
        $prnum = sprintf "%04d",$k ;
        $folder = "pr".$prnum;
    } else {
        $folder = $ARGV[0];
    }
    print " Dir: $folder\n";

    # position of dimer
    for($i=0; $i<$total_atoms; $i++) {
        for ($j=0; $j<3; $j++) {
            $coordinatesnew->[$i][$j] = $coordinates->[$i][$j];
        }
    }
    for $i (@$atoms2bdisplaced) {
        for ($j=0; $j<3; $j++) {
            $coordinatesnew->[$i][$j] = $coordinatesnew->[$i][$j] + $sigma->[$i][$j]*gauss;
        }
    }
    $coordinatesnew = kardir($coordinatesnew,$basis,$lattice,$total_atoms);
    set_bc($coordinatesnew,$total_atoms);
    $coordinatesnew = dirkar($coordinatesnew,$basis,$lattice,$total_atoms);
    # spring-relaxation: prevent atoms from being too close 
    print "relaxing the initial guess ... ... \n";
    spring_relaxation($coordinatesnew,$basis,$lattice,$lattice_vec,$total_atoms,$selective);
    print "... done\n";
    # displacement direction and magnitude
    $r_mag = 0.0;
    for($i=0; $i<$total_atoms; $i++) {
        for ($j=0; $j<3; $j++) {
            $r_vector->[$i][$j] = $coordinatesnew->[$i][$j]-$coordinates->[$i][$j];
        }
        ($r_vector->[$i],$mag_difference) = FindMinimumImage($r_vector->[$i],$lattice_vec);
        $r_mag += $mag_difference;
    }
    ($r_mag>0) || die " Error in diminit.pl: No atoms have been displaced.\n";
    $r_mag = 1.0/sqrt($r_mag);
    for($i=0; $i<$total_atoms; $i++) {
        $r_vector->[$i] = VecMultiply($r_vector->[$i],$r_mag);
    }
    # write dimer images
    if(-e $folder) { die " Error from diminit.pl: $folder already exists.\n"; }
    system "mkdir $folder";
    $outputfilename = "ciPOSCAR";
    $coordinatesnew = kardir($coordinatesnew,$basis,$lattice,$total_atoms);
    write_poscar($coordinatesnew,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective,$description,$outputfilename,$filetype);
    system "mv ciPOSCAR $folder/POSCAR";
    write_othercar($r_vector,$total_atoms,$outputfilename);
    system "mv ciPOSCAR $folder/MODECAR";
    system "cp KPOINTS POTCAR $folder/";

    # fix this later.
    if(-e "INCAR") { system "cp INCAR $folder/"; }
    if(-e "akmc.sub") { system "cp akmc.sub $folder/"; }
}

# ---------------------------------------------------------------------------------------------------------
# get all the candidate atoms and the neighborlist between them: ($NcList,$NN_DisplaceList)=
# BuildDisplaceList($R,$lattice_vec,$DisplaceList,$DisplaceAlgo,$DisplaceRange,$rcut,$MaxCord);
# ---------------------------------------------------------------------------------------------------------
sub BuildDisplaceList{
    my ($R,$lattice_vec,$DisplaceList,$DisplaceAlgo,$DisplaceRange,$rcut,$MaxCord)=@_;
    my ($NN_NcordFile,$NN,$NCord,$i,$j,$good,$NcList,$NN_DisplaceList,$dummy);
    print "Rcut=$rcut\n";
    $NcList = {};
    $NN_DisplaceList = {};
    $NN_NcordFile = "neighborlist.dat";
    $good = 0;
    $j = @$R;
    if(-e $NN_NcordFile){
        ($NN,$NCord,$good) = readNN_Ncord($NN_NcordFile,$rcut,$j);
    }
    if(!$good){ # need to make a new list
        ($NN,$NCord) = FindNeighbors($R,$rcut,$lattice_vec);
        print "**********************************\n";
        for $j (sort {$a <=> $b} keys %$NN){
            print $j."  ".$NCord->{$j}."  @{$NN->{$j}}\n";
        }
        writeNN_Ncord($NN_NcordFile,$NN,$NCord,$rcut);
    }
    $NN_DisplaceList = FindDLNeighbors($R,$DisplaceList,$DisplaceRange,$lattice_vec);
    SWITCH: {
        if($DisplaceAlgo == 0) { $NcList = FindLowestCoordAtom($DisplaceList,$NCord); last SWITCH; }
        if($DisplaceAlgo == 1) { $NcList = FindLessMaxcord($DisplaceList,$NCord,$MaxCord);last SWITCH; }
        if($DisplaceAlgo == 2 || $DisplaceAlgo == 3) {
            $NcList = FindAtomsInIslands($R,$DisplaceList,$NN,$NCord,$lattice_vec,$rcut,$DisplaceAlgo,$MaxCord);
            last SWITCH;
        }
        $NcList = $DisplaceList; #default: select a random atom
    }
    return ($NcList,$NN_DisplaceList);
}

sub FindNeighbors{
    my ($R,$rcut,$lattice_vec) = @_;
    my ($i,$j,$k,$total_atoms,$difference,$mag_difference);
    my ($NCord,$NN,@dummy);
    print "rcut=$rcut\n";
    $rcut = $rcut*$rcut;
    $total_atoms = @$R-1; 
    for($i=0; $i<=$total_atoms; $i++){
        $NCord->{$i} = 0;
        $NN->{$i} = [()];
    }
    for($i=0; $i<$total_atoms; $i++){
        for($j=$i+1; $j<=$total_atoms; $j++){
            for($k=0;$k<3;$k++){
                $difference->[$k] = $R->[$i][$k]-$R->[$j][$k];
            }
            ($difference, $mag_difference) = FindMinimumImage($difference, $lattice_vec);
            if($mag_difference < $rcut){
                $NCord->{$i}++;
                push @{$NN->{$i}}, $j;
                $NCord->{$j}++;
                push @{$NN->{$j}}, $i;
            }
        }
    }
    @dummy = sort {$a <=> $b} values %$NCord;
    $NCord->{"mincord"} = $dummy[0];
    return ($NN,$NCord);
}

sub FindDLNeighbors{
    my ($R,$DisplaceList,$rcut,$lattice_vec)=@_;
    my ($i,$j,$k,$m,$n,$total_atoms,$difference,$mag_difference,$NN);
    $total_atoms = @$DisplaceList-1;
    print "rcut=$rcut\n";
    $rcut = $rcut*$rcut;
    $NN = {};
    for $i (@$DisplaceList){
        $NN->{$i} = [()];
    }
    for($m=0; $m<$total_atoms; $m++){
        $i = $DisplaceList->[$m];
        for($n=$m+1; $n<=$total_atoms; $n++){
            $j = $DisplaceList->[$n];
            for($k=0; $k<3; $k++){
                $difference->[$k] = $R->[$i][$k]-$R->[$j][$k];
            }
            ($difference, $mag_difference) = FindMinimumImage($difference, $lattice_vec);
            if($mag_difference < $rcut){
                push @{$NN->{$i}}, $j;
                push @{$NN->{$j}}, $i;
            }
        }
    }
    return $NN;
}

sub writeNN_Ncord{
    my ($NN_NcordFile,$NN,$NCord,$rcut) = @_;
    my ($i,$mincord_total);
    open(ISMIN,">$NN_NcordFile") || die "cannot open $NN_NcordFile\n";
    print ISMIN "Rcut  ".$rcut."\n";
    print ISMIN "mincord  ".$NCord->{"mincord"}."\n";
    print ISMIN "atom\t"."ncoord\t"."neighbors\n";
    for $i (sort {$a <=> $b} keys %$NN){
        print ISMIN $i."\t".$NCord->{$i}."\t"."@{$NN->{$i}}\n";
    }
    close ISMIN;
    if(($i+1) == $total_atoms){
        die "Error in writeNN_Ncord: # of neighborlist elements is not same as # of atoms from POSCAR\n";
    }
}

sub readNN_Ncord{
    my ($NN_NcordFile,$rcut,$total_atoms)=@_;
    my ($NCord,$NN,@line,$line,$i,$j,$numatoms,$good);
    $NN = {};
    $NCord = {};
    $i = $j = $good = 0;
    open(ISMIN,"<$NN_NcordFile") || die "cannot open $NN_NcordFile\n";
    while($line = <ISMIN>){
        $line =~ s/^\s+//;
        if($line eq "\n") { next; }
        @line = split(/\s+/,$line);
        if(lc($line[0]) eq "rcut"){
            if($line[1] =~ /^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$/ && $line[1]==$rcut) {
                $good = 1;
                $i++;
                next;
            } else {
                last;
            }
        }
        if($i==1 && lc($line[0]) eq "mincord"){
            if($line[1] =~ /^\d+$/) { $NCord->{"mincord"}=$line[1]; }
            $i++;
            next;
        }
        if($i==2 && lc($line[0]) eq "atom"){
            $i++;
            next;
        }
        if($i==3){
            $numatoms = @line-1;
            $NCord->{$line[0]} = $line[1];
            $NN->{$line[0]} = [(@line[2 .. $numatoms])];
            #print $line[0]." : "."@{$NN->{$line[0]}}.\n";
            $j++;
        }
    }
    close ISMIN;
    if($j && $j != $total_atoms){$good=0;}
    return ($NN,$NCord,$good);
}

sub FindLowestCoordAtom{
    my ($DisplaceList,$NCord) = @_;
    my $NcList = ();
    my ($Nc,$i,$mincord);
    $mincord = $NCord->{$DisplaceList->[0]};
    for $i (@$DisplaceList) {
        $Nc = $NCord->{$i};
        if($Nc < $mincord){
            $mincord = $Nc;
            $NcList = ();
            push @$NcList, $i;
        }elsif($Nc == $mincord){
            push @$NcList, $i;
        }
    }
    if($NCord->{"mincord"} < $mincord){
        print "warning; the overall lowest coordinated atom(s) not included in the displacement. Make sure that is what you want.\n";
    }
    return $NcList;
}

sub FindLessMaxcord{
    my ($DisplaceList,$NCord,$MaxCord) = @_;
    my $NcList = ();
    my $i = 0;
    for $i (@$DisplaceList) {
        if($NCord->{$i} < $MaxCord){
            push @$NcList, $i;
        }
    }
    if(@$NcList == 0){
        print "MaxCord is too small($MaxCord) and no atoms in DISPLACECAR_sp are selected. Do random selection\n";
        $NcList = $DisplaceList;
    }
    return $NcList;
}

# ---------------------------------------------------------------------------------------------------------
# get a list of atoms to be displaced :  $atoms2bdisplaced=
# BuildNewDisplacecar($DisplaceList,$NcList,$total_atoms,$DisplaceAlgo,$DisplaceRange,$rcut);
# ---------------------------------------------------------------------------------------------------------
sub BuildNewDisplacecar{
    my ($DisplaceList,$NcList,$NN_DisplaceList,$DisplaceAlgo,$randomnum)=@_;
    my ($atoms2bdisplaced,$Nc,$atom);
    #pick out one displacing mechanism, currently it looks trivial
    SWITCH: {
        if($DisplaceAlgo == 0) { $Nc = SelectNcAtom($NcList,$DisplaceAlgo,$randomnum); last SWITCH; }
        if($DisplaceAlgo == 1) { $Nc = SelectNcAtom($NcList,$DisplaceAlgo,$randomnum); last SWITCH; }
        if($DisplaceAlgo == 2) { $Nc = SelectNcAtom($NcList,$DisplaceAlgo,$randomnum); last SWITCH; }
        $Nc = SelectNcAtom($NcList,$DisplaceAlgo,$randomnum);  # default: displace around a random atom in DisplaceList
    }
    $atoms2bdisplaced = [($Nc)];
    for $atom (@{$NN_DisplaceList->{$Nc}}){ push @$atoms2bdisplaced, $atom; }
    print "Nc=$Nc leads to atoms to be displaced:  "."@$atoms2bdisplaced\n";
    return $atoms2bdisplaced;
}

sub SelectNcAtom{
    my ($NcList,$DisplaceAlgo,$randomnum)=@_;
    my ($i,$j,$Nc,$total);
    #print "randomnum=$randomnum\n";
    $j = @$NcList;
    $j = int($j*$randomnum);
    print "j=$j\n";
    for($i=0; $i<@$NcList; $i++){
        if(!$j){
            $Nc = $NcList->[$i];
            last;
        } else {
            $j--;
        }
    }
    return $Nc;
}

# ---------------------------------------------------------------------------------------------------------
# relax an initial saddle point guess to prevent atoms from getting too close to each other
# spring_relaxation($coordinates,$basis,$lattice,$lattice_vec,$totalatoms,$selective);
# ---------------------------------------------------------------------------------------------------------
sub spring_relaxation{
    my ($R,$basis,$lattice,$lattice_vec,$totalatoms,$selective) = @_;
    my $epsilon = 0.1;
    my $cutoff = 0.1;
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
                if($Rij < $rcut){
                    $goodsound = 0;
                    $Rij = sqrt($Rij);
                    ($fij,$phiij) = linear_repulsion_pot($Rij,$cutoff);
                    #($fij,$phiij) = exp_repulsion_pot($Rij,$cutoff);
                    $fij = $fij/$Rij;
                    for($k=0; $k<3; $k++){
                        $phiij = $fij*$difference->[$k];
                        $force->[$i][$k] -= $phiij;
                        $force->[$j][$k] += $phiij;
                    }
                }
            }
        }
        if($goodsound) {
            $converged = 1;
            next;
        }
        for($i=0; $i<$totalatoms; $i++) {
            for($j=0; $j<3; $j++) {
                if(lc($frozen->[$i][$j]) eq "f") { next; }
                $dr = $epsilon*$force->[$i][$j];
                if(abs($dr) > $drmax) { $dr = $drmax*$dr/abs($dr); }
                #print "$i: $j -- dr: $dr\n"; 
                $R->[$i][$j] += $dr;
            }
        }
        if($ORTHOGONAL){
            for($i=0;$ i<$totalatoms; $i++){
                for($j=0; $j<3; $j++){
                    while($R->[$i][$j] > $BOX[$j]){ $R->[$i][$j] -= $BOX[$j]; }
                    while($R->[$i][$j] < $BOX_[$j]){ $R->[$i][$j] += $BOX[$j]; }
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
  my $r=shift;
  my $cutoff=shift;
  my $slope=-0.3;
  my ($energy, $fij);
  my $intercept=$slope*$cutoff;
  $energy=$slope*$r-$intercept;
  $fij=$slope;
  return ($fij,$energy);
} 

sub exp_repulsion_pot{ # E=exp(kx)-b,k<0
    my $r = shift;
    my $cutoff = shift;
    my $slope = -0.3;
    my ($energy, $fij);
    my $intercept = exp($cutoff*$slope);
    my $item = exp($slope*$r);
    $energy = $item-$intercept;
    $fij = $slope*$item;
    return ($fij,$energy);
} 

# ---------------------------------------------------------------------------------------------------------
# divide adsorbates on the surface into islands and choose atoms from each island: $NcList=FindAtomsInIslands
# ($R,$DisplaceList,$NN,$NCord,$lattice_vec,$rcut,$DisplaceAlgo,$MaxCord)
# ---------------------------------------------------------------------------------------------------------
sub FindAtomsInIslands{
    my ($R,$DisplaceList,$NN,$NCord,$lattice_vec,$DisplaceAlgo,$MaxCord) = @_;
    my ($islandsinfo,$islands,$NumIslands,$Nc,$total_atoms,$good,$mincord_island,$NcList);
    my ($i,$j,$rcut2);
    $islandsinfo = "islands.dat";
    $total_atoms = @$DisplaceList;
    if(-e $islandsinfo){
        ($islands,$mincord_island,$good) = readIslands($islandsinfo,$rcut,$total_atoms);
    }
    if(!$good){ # create a new island information
        print "create a new island/coordination list\n";
        $rcut2 = $rcut*$rcut;
        $islands = FindSurfIslands($R,$DisplaceList,$rcut2,$lattice_vec);
        $mincord_island = ();
        $NumIslands = @$islands;
        for($i=0; $i<$NumIslands; $i++){      
            $mincord_island->[$i] = FindMostUndercoordianted($i,$islands,$NCord);
            print "island $i: @{$islands->[$i]}\n";
            print "mincord in island $i is ".$mincord_island->[$i]."\n";
        }
        for($i=0; $i<@$DisplaceList; $i++){
            $j = $DisplaceList->[$i];
            print "atom $j: NCord=".$NCord->{$j}." @{$NN->{$j}}\n";
        }
        writeIslands($islandsinfo,$islands,$mincord_island,$rcut);
    } else {
        print "read in old islands and minimumcord information from file $island_minncord\n";
    }
    $NcList = ();
    $NumIslands = @$islands;
    for($i=0; $i<$NumIslands; $i++){
        for $j (@{$islands->[$i]}){
            if($DisplaceAlgo == 2){
                if($NCord->{$j} == $mincord_island->[$i]){ push @$NcList,$j; }
            } elsif($DisplaceAlgo == 3){
                if($NCord->{$j} < $MaxCord){ push @$NcList,$j; }
            }
        }
    }
    if(@$NcList == 0){
        print "MaxCord is too small($MaxCord) in FindAtomsInIslands. Have to choose atoms randomly\n";
        $NcList = $DisplaceList;
    }
    return $NcList;
}

sub FindSurfIslands{
    my ($R,$DisplaceList,$rcut2,$lattice_vec) = @_;
    my ($i,$domains,$mark,$islands);
    $domains = {};
    $mark = 0;
    for $i (@$DisplaceList){
        if(!(exists($domains->{$i}))){
            $islands->[$mark] = [()];
            $domains->{$i} = $mark;
            MarkMyNeighbors($R,$i,$domains,$DisplaceList,$rcut2,$lattice_vec);
            $mark++;
        }
    }
    for $i (keys %$domains){
        push @{$islands->[$domains->{$i}]},$i;
    }
    return $islands;
}

sub MarkMyNeighbors{ # a recursive subroutine
    my ($R,$j,$domains,$DisplaceList,$rcut2,$lattice_vec)=@_;
    my ($i,$m,$k,$difference,$mag_difference,$rcut2,$lattice_vec);
    for $i (@$DisplaceList){
        if($i != $j){ #my neighbors
            for($k=0; $k<3; $k++){
                $difference->[$k] = $R->[$i][$k]-$R->[$j][$k];
            }
            ($difference, $mag_difference) = FindMinimumImage($difference,$lattice_vec);
            if($mag_difference < $rcut2){
                if(!(exists($domains->{$i}))){
                    $domains->{$i} = $domains->{$j};
                    MarkMyNeighbors($R,$i,$domains,$DisplaceList,$rcut2,$lattice_vec);
                }
            }
        }# end of outer if block
    }
}

sub FindMostUndercoordianted{
    my ($i,$islands,$NCord)=@_;
    my ($mincord,$j);
    $mincord = $NCord->{$islands->[$i][0]};
    for $j (@{$islands->[$i]}){
        if($NCord->{$j} < $mincord){
            $mincord = $NCord->{$j};
        }
    }
    return $mincord;
}

sub readIslands{
    my ($islandsinfo,$rcut,$total_atoms) = @_;
    my ($good,@line,$line,$i,$j,$total,$numatoms,$islands,$mincord_island);
    open(ISMIN,"<$islandsinfo") || die "cannot open $islandsinfo\n";
    $i = $j = $good = $total = 0;
    while($line = <ISMIN>){
        $line =~ s/^\s+//;
        if($line eq "\n"){next;}
        @line = split(/\s+/,$line);
        if(lc($line[0]) eq "rcut"){
            if($line[1]=~/^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$/ && $line[1]==$rcut){
                $good = 1;
                $i++;
                next;
            } else {
                last;
            }
        }
        if($i==1 && lc($line[0]) eq "island"){
            $j=0;
            $i++;
            next;
        }
        if($i==2){
            $numatoms = @line-1;
            $mincord_island->[$j] = $line[1];
            $islands->[$j] = [(@line[2 .. $numatoms])];
            $j++;
            $total = $total + $numatoms-1;
        }
    }
    close ISMIN;
    if($total != $total_atoms){ $good=0; }
    return ($islands,$mincord_island,$good);
}

sub writeIslands{
    my ($islandsinfo,$islands,$mincord_island,$rcut)=@_;
    my ($i,$j);
    open(ISMIN,">$islandsinfo") || die "cannot open $islandsinfo\n";
    print ISMIN "Rcut ".$rcut."\n";
    print ISMIN "Island\tMinCoord\tAtomsInIsland\n";
    for($i=0; $i<@$islands; $i++){
        print ISMIN $i."\t".$mincord_island->[$i]."\t"."@{$islands->[$i]}\n";
    }
    close ISMIN;
}

#----------------------------------------------------------------------
# FindMinimumImage: (INPUT) A vector (pointer) (OUTPT) the minimum image vector
# This piece of code was borrowed from Bader code (F90) by Wenjie and Graeme
#----------------------------------------------------------------------
sub FindMinimumImage{
    my ($dr_car,$lattice_vec) = @_;
    my ($dsqmin,$done,$v1,$v2,$v3,$drt_car,$dsq);
    my ($i,$j,$k);
    if($ORTHOGONAL){
        for($i=0; $i<@$dr_car; $i++){
            if($dr_car->[$i]<$HALFBOX_[$i]){ $dr_car->[$i] += $BOX[$i]; }
            if($dr_car->[$i]>$HALFBOX[$i]){ $dr_car->[$i] -= $BOX[$i]; }
        }
        $dsqmin = VecDot($dr_car,$dr_car);
    } else {
        $dsqmin = VecDot($dr_car,$dr_car);
        $done = 1;
        while($done) {
            $found = 1;
            for($i=-1; $i<2; $i++){
                $v1 = VecMultiply($lattice_vec->[0],$i);
                for($j=-1; $j<2; $j++){
                    $v2 = VecMultiply($lattice_vec->[1],$j);
                    for($k=-1; $k<2; $k++){
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
    print sqrt($dsqmin)."\n";
    return ($dr_car,$dsqmin);
}

#----------------------------------------------------------------------
# VecDot: (INPUT) two vectors (pointers) (OUTPUT) a dot product
#----------------------------------------------------------------------
sub VecDot{
    my ($vec1, $vec2)=@_;
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
