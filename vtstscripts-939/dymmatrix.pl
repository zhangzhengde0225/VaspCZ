#!/usr/bin/env perl
#;-*- Perl -*-

# This program extracts the dynamical matrix information (the forces vs displacement)
# from the OUTCAR files.

use FindBin qw($Bin);
use lib "$Bin";
use Vasp;

use Math::MatrixReal ;

$defaultflag = 0;
if(@ARGV == 0) {
    $outcarfiles[0] = "OUTCAR";
    $num_directories = 1;
    $num_displacecars = 1;
    $singleflag = 1; # One DISPLACECAR and OUTCAR
    $defaultflag = 1; # default names: DISPLACECAR and OUTCAR
} elsif(@ARGV == 2) {
    # Just take in one DISPLACECAR and one OUTCAR
    # die "INPUT A DISPLACECAR AND A OUTCAR ! \n" , if @ARGV != 2;
    @outcarfiles = $ARGV[1];
    $num_directories = 1;
    $num_displacecars = 1;
    $singleflag = 1;
} else {
    # Older input !
    @ARGV>=3 || die "usage: dymmatrix.pl <number DISPLACECARS> <DISPLACECAR 1> <...> <OUTCAR 1 for DISPLACECAR 1> <OUTCAR 2 for DISPLACECAR 1> <...>\n";
    $num_displacecars = @ARGV[0];
    @outcarfiles = @ARGV[$num_displacecars+1..@ARGV-1];
    # print "Number of OUTCAR files: ",scalar(@outcarfiles),"\n";
    $num_directories = @outcarfiles/$num_displacecars;
    $singleflag = 0;
}

# Get the masses for each atom from the OUTCAR
#  $pomass = `grep POMASS OUTCAR | grep ZVAL`;
print "Reading masses from OUTCAR\n";
$outfilename=$outcarfiles[0];
$pomass = `grep POMASS $outfilename | grep ZVAL`;
@pomass = split /\n/ , $pomass;
$ntyp = @pomass;

# Get the number of each type of atom 
#  $nions_typ = `grep 'ions per type =' OUTCAR`;
$nions_typ = `grep 'ions per type =' $outfilename`;
$line = $nions_typ;
$line =~ s/^\s+//g;
@line = split /\s+/,$line;
@natoms_typ = @line[4..4+$ntyp];
for($i=0; $i<$ntyp; $i++) {
    $line = $pomass[$i];
    $line =~ s/^\s+//g;
    @line = split /\s+/,$line;
    @m = split /;/ , $line[2];
    $typ_mass[$i] = $m[0];
    $natoms_tot += $natoms_typ[$i];
}
# print masses
print join("\n",@typ_mass),"\n";

$j = 0;
$t = $natoms_typ[0];
for($i=0; $i<$natoms_tot; $i++) {
    if($i >= $t) {
        $j++;
        $t += $natoms_typ[$j];
    }
    $mass[$i] = $typ_mass[$j];
}

#  print "Mass: ",@mass ;

print "Reading forces from OUTCAR\n";
# $splitline = 'POSITION.*TOTAL-FORCE \(eV\/Angst\)';
$splitline = 'HIPREC TOTAL-FORCE \(eV\/A\)';

#----------------------------------------------------------------------
# extract the forces out of the OUTCAR files for this DISPLACECAR
#----------------------------------------------------------------------
$num_dof = 0;
for ($i=0; $i<@outcarfiles; $i++) {
    $outfilename = $outcarfiles[$i];
    # print "OUTCAR: $outfilename\n";

    $outcar = "";
    open (IN,$outfilename);
    while (<IN>) { $outcar .= $_; }
    close (IN);

    @outcar = split(/$splitline/,$outcar);
    shift @outcar;
    $time_steps = @outcar;
#  print "OUTCAR: $outfilename, $time_steps\n";
    for ($j=0; $j<$time_steps; $j++) {
        $outcar[$j] =~ s/--+//g;
        $outcar[$j] =~ s/total drift.*//ms;
        $outcar[$j] =~ s/\n\s+/\n/g;
        $outcar[$j] =~ s/^(\s*\n)+//g;
        $outcar[$j] =~ s/(\n\s*)*$//g;
    }

    @forcelines = split(/\n/,$outcar[$0]);
    for ($k=0; $k<@forcelines; $k++) {
        @forcecomps = split(/\s+/,$forcelines[$k]);
        for ($l=0; $l<3; $l++) {
            # $force_ref->[$k][$l] = $forcecomps[$l+3];
            $force_ref->[$k][$l] = $forcecomps[$l];
        }
    }

    for ($j=1; $j<$time_steps; $j++) {
        @forcelines = split(/\n/,$outcar[$j]);
        for ($k=0; $k<@forcelines; $k++) {
            @forcecomps = split(/\s+/,$forcelines[$k]);
            for ($l=0; $l<3; $l++) {
                # $forces->[$num_dof][$k][$l] =- ($forcecomps[$l+3]-$force_ref->[$k][$l]);
                $forces->[$num_dof][$k][$l] =- ($forcecomps[$l] - $force_ref->[$k][$l]);
            }
        }
        $num_dof++;
    }
}

#----------------------------------------------------------------------
# read displacecars, count displacements
#----------------------------------------------------------------------
print "Reading displacements from DISPLACECAR\n";
$num_displacements = 0;
for ($m=0; $m<$num_displacecars; $m++) {
    if($singleflag){
        if($defaultflag){
            $filename = "DISPLACECAR";
        } else {
            $filename = $ARGV[0];
        }
    } else {
        $filename = $ARGV[$m+1];
    }
    ($displacecar->[$m],$total_atoms) = read_othercar($filename);
    for ($i=0; $i<$total_atoms; $i++) {
        for ($j=0; $j<3; $j++) {
            if ($displacecar->[$m][$i][$j] != 0) {
                $num_displacements++;
            }
        }
    }
}

#print "Num Displacements: $num_displacements\n";
#print "Total atoms: $total_atoms\n";

#______________________________________________________________________
# check which forces correspond to displacements, build matrix
#----------------------------------------------------------------------
print "Building Hessian matrix\n";
$current_displacement = 0;
for ($m=0; $m<$num_displacecars; $m++) {
    for ($i=0; $i<$total_atoms; $i++) {
        for ($j=0; $j<3; $j++) {
            if ($displacecar->[$m][$i][$j] != 0) {
                $this_displacement = 0;
                for ($n=0; $n<$num_displacecars; $n++) {
                    for ($k=0; $k<$total_atoms; $k++) {
                        for ($l=0; $l<3; $l++) {
                            if ($displacecar->[$n][$k][$l] != 0) {
                                $matrix->[$current_displacement][$this_displacement]
                                       = $forces->[$current_displacement][$k][$l]/$displacecar->[$m][$i][$j];
                                $matrix->[$current_displacement][$this_displacement] /= (sqrt($mass[$k]*$mass[$i]));
                                $this_displacement++;
                            }
                        }
                    }
                }
                $current_displacement++;
            }
        }
    }
}

#print "writing matrix\n";
#for ($i=0;$i<@{$matrix};$i++) {
#    for ($j=0;$j<@{$matrix->[$i]};$j++) {
#        printf "%12.8f", $matrix->[$i][$j]."  ";
#    }
#    print "\n";
#}

# Write the mass-scaled Hessian matrix
# Temporary fix:
$num_dof = $num_displacements;

$B = Math::MatrixReal->new($num_dof,$num_dof);
open FRQ , ">freq.mat";
for($i=0; $i<$num_dof; $i++) {
    $B->assign($i+1,$i+1,$matrix->[$i][$i]);
    for($j=$i+1; $j<$num_dof; $j++) {
        $B->assign($i+1,$j+1,0.5*($matrix->[$i][$j]+$matrix->[$j][$i]));
        $B->assign($j+1,$i+1,0.5*($matrix->[$i][$j]+$matrix->[$j][$i]));
    }
    for ($j=0; $j<@{$matrix->[$i]}; $j++) {
        printf FRQ "%16.8f", $B->element($i+1,$j+1)."  ";
    }
    print FRQ "\n";
}
close FRQ;

# Solve the eigensystem
print "Solving the eigensystem\n";
($l, $V) = $B->sym_diagonalize();

# Sort the eigenvalues. Use selection sort. $t will be the order we want to print out
for($i=0; $i<$num_dof; $i++) {
    $w[$i] = $l->element($i+1,1);
    $t[$i] = $i; 
}
for($i=0 ; $i<$num_dof ; $i++) {
    $iptr = $i ;
    for($j=$i+1; $j<$num_dof; $j++) {
        if($w[$j] < $w[$iptr]){ $iptr = $j; }
    }
    if($i != $iptr) {
        $tmp = $w[$i];
        $w[$i] = $w[$iptr];
        $w[$iptr] = $tmp;
        $tmp = $t[$i];
        $t[$i] = $t[$iptr];
        $t[$iptr] = $tmp;
        # Sort the eigenvectors as well
        $V->swap_col($i+1,$iptr+1);
    } 
}

# Write out the eigenvalues, frequencies in cm^{-1} to the STDOUT and 
# omega^{2} to 'eigs.dat'
print "Writing the eigenvalues\n";
open EIG , ">eigs.dat";
open FREQ , ">freq.dat";
print "\n";
for($i=0; $i<$num_dof; $i++) {
    if($w[$i] > 0) {
        printf "%12.6f cm^{-1} ... 0 \n",sqrt($w[$i]) * 521.47;
        printf FREQ "%12.6f cm^{-1} ... 0 \n",sqrt($w[$i]) * 521.47;
    }else {
        printf "%12.6f cm^{-1} ... 1 \n",sqrt(abs($w[$i])) * 521.47;
        printf FREQ "%12.6f cm^{-1} ... 1 \n",sqrt(abs($w[$i])) * 521.47;
    }
    printf EIG "%25.15g \n",$w[$i];
}
print "\n";
close EIG;
close FREQ;

# Write out the modes
open MDE , ">modes.dat";
for($i=0; $i<$num_dof; $i++) {
    for ($j=0; $j<@{$matrix->[$i]}; $j++) {
        printf MDE "%16.8f", $V->element($i+1,$j+1)."  ";
    }
    print MDE "\n";
}
close MDE;

