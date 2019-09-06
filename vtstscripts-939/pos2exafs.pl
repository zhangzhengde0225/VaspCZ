#!/usr/bin/env perl
#;-*- Perl -*-

use FindBin qw($Bin);
use lib "$Bin";
use Vasp;

@ARGV>=3 || die "Usage: pos2exafs.pl <POSCAR> <BinSize> <outfile> <choosetype 0 or 1> <atomtype1> <atomtype2>\n";

$poscar = $ARGV[0];
$binsize = 0.05;
if(@ARGV>1) { $binsize = $ARGV[1]; }
$outfile = "exafs.dat";
if(@ARGV>2) { $outfile = $ARGV[2]; }
$choosetype = 0;
if(@ARGV>3) { $choosetype = $ARGV[3]; }
if($choosetype==1 && @ARGV<6) { die "please choose two atom types.\n"; }
if(@ARGV>4) {
    $type1 = $ARGV[4];
    if(!($type1 =~ /^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$/)) {
        die "type1:$type1 is not a number.\n";
    }
}
if(@ARGV>5) {
    $type2 = $ARGV[5];
    if(!($type2 =~ /^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$/)) {
        die "type2:$type2 is not a number.\n";
    }
}

coord_number($poscar, $binsize, $outfile);

sub coord_number{
    my ($poscarfile, $bin_size, $outfile) = @_;
    my ($k,$i,$j, $difference,$mag_difference,$bins,$index,$number);
    my ($coordinates,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective);
    if($bin_size == 0) {
        die "Error in coord_number, binsize is $bin_size.\n";
    }

    # read in the poscar file
    ($coordinates,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective) = read_poscar($poscarfile);

    # count the number of neighbors for each atom
    if($choosetype==0) {
        for ($k=0; $k<$total_atoms-1; $k++) {
            for ($i=$k+1; $i<$total_atoms; $i++) {
                for ($j=0; $j<3; $j++) {
                    $difference->[0][$j] = pbc($coordinates->[$k][$j] - $coordinates->[$i][$j]);
                }
                $difference = dirkar($difference,$basis,$lattice,1);
                $mag_difference = magnitude($difference,1);
                $index = int($mag_difference/$bin_size);
                $bins->{$k}->{$index}[@{$bins->{$k}->{$index}}] = $i;
                $bins->{$i}->{$index}[@{$bins->{$i}->{$index}}] = $k;
            }
        }
    } elsif($choosetype==1) {
        # set the boundary of two types of atoms
        $typenum = @{$num_atoms};
        $t1min = $t1max = $t2min = $t2max = 0.0;
        for ($i=0; $i<=$type1-2; $i++) {
            $t1min = $t1min + $num_atoms->[$i]; }
        $t1max = $t1min + $num_atoms->[$type1-1];
        for ($i=0; $i<=$type2-2; $i++) {
            $t2min = $t2min + $num_atoms->[$i]; }
        $t2max = $t2min + $num_atoms->[$type2-1];
        print "boundary of type1:,$t1min,  $t1max\n";
        print "boundary of type2:,$t2min,  $t2max\n";

        for ($k=$t1min; $k<$t1max; $k++) {
            for ($i=$t2min; $i<$t2max; $i++) {
                for ($j=0; $j<3; $j++) {
                    $difference->[0][$j] = pbc($coordinates->[$k][$j] - $coordinates->[$i][$j]);
                }
                $difference = dirkar($difference,$basis,$lattice,1);
                $mag_difference = magnitude($difference,1);
                $index = int($mag_difference/$bin_size);
                $bins->{$k}->{$index}[@{$bins->{$k}->{$index}}] = $i;
            }
        }
    }
  
    foreach $k (sort {$a<=>$b} keys %$bins) {
        @after_sort = (sort {$a<=>$b} keys %{$bins->{$k}});
        $min_dist = $after_sort[0];
        $len = @after_sort;
        $i = 0;
        while ($i<$len && $after_sort[$i]< 2*$min_dist) {
            $index = $after_sort[$i];
            $natom = @{$bins->{$k}->{$index}};
            if(!(exists($exafs->{$index}))) {
                $exafs->{$index} = 0;
            }
            $exafs->{$index} = $exafs->{$index}+$natom;
            $i = $i+1;
        }
    }

    open(OUT,">$outfile");
    print OUT "Distribution of 1st NN:\n";
    @after_sort = (sort {$a<=>$b} keys %$exafs);
    $max_bin = $after_sort[-1];
    $none = 0;
    for($i=1; $i<=$max_bin; $i++){
        if(exists($exafs->{$i})) {
            printf OUT "%5.3f    %3i  \n", $i*$bin_size,$exafs->{$i};
        } else {
            printf OUT "%5.3f    %3i  \n", $i*$bin_size,$none;
        }
    }
    close(OUT);
}
