#!/usr/bin/env perl
#;-*- Perl -*-

# This program extracts the dynamical matrix information (the forces vs displacement)
# from the OUTCAR files.  

use FindBin qw($Bin);
use lib "$Bin";
use Vasp;

@args = @ARGV;
@args >= 4 || die "usage: dymreorder.pl <number of reference DISPLACECARS> <list of DISPLACECARs> <number of new DISPLACECARs> <list of DISPLACECARs> <new matrix>\n";

$numdispcar1 = $args[0];
$numdispcar2 = $args[$args[0]+1];

#----------------------------------------------------------------------
# read 1st set of displacecars, count displacements
#----------------------------------------------------------------------
$num_displacements = 0;
for ($m=0;$m<$numdispcar1;$m++) {
    $filename=$args[$m+1];
    ($displacecar1->[$m],$total_atoms) = read_othercar($filename);
    for ($i=0; $i<$total_atoms; $i++) {
        for ($j=0;$j<3;$j++) {
            if ($displacecar1->[$m][$i][$j] != 0) {
                $num_displacements1++;
            }
        }
    }
}  

#----------------------------------------------------------------------
# read 2st set of displacecars, count displacements
#----------------------------------------------------------------------
$num_displacements = 0;
for ($m=0; $m<$numdispcar2; $m++) {
    $filename = $args[$m+$args[0]+2];
    ($displacecar2->[$m],$total_atoms) = read_othercar($filename);
    for ($i=0; $i<$total_atoms; $i++) {
        for ($j=0; $j<3; $j++) {
            if ($displacecar2->[$m][$i][$j] != 0) {
                $num_displacements2++;
            }
        }
    }
}  

($num_displacements1==$num_displacements2) || 
  die "DISPLACECARs do not have same number of displacements!\n";

#----------------------------------------------------------------------
# read matrix
#----------------------------------------------------------------------
#print $args[@args-1]."\n";
$j = 0;
open (IN,$args[@args-1]) || die ("Can't open file $!");
while (<IN>) {
    $_ =~ s/^\s+//;
    @temp = split(/\s+/,$_);
    $kk = @temp;
    for ($k=0; $k<@temp; $k++) {
        $matrix->[$j][$k] = $temp[$k];
    }
    $j++;
}
close(IN);

#print "$j $num_displacements1\n";
($j==$num_displacements1) || die "Matrix does not have same number of displacements!\n";

#----------------------------------------------------------------------
# reorder matrix
#----------------------------------------------------------------------

#______________________________________________________________________
# check which forces correspond to displacements, build matrix
#----------------------------------------------------------------------
$index_a1 = 0;
for ($m=0; $m<$numdispcar2; $m++) {
    for ($i=0; $i<$total_atoms; $i++) {
        for ($j=0; $j<3; $j++) {
            if ($displacecar2->[$m][$i][$j] != 0) {
                $index_b1 = 0;
                for ($n=0; $n<$numdispcar1; $n++) {
                    for ($k=0; $k<$total_atoms; $k++) {
                        for ($l=0; $l<3; $l++) {
                            if ($displacecar1->[$n][$k][$l] != 0) {
                                if ($k==$i && $l==$j) {
                                    $index_a2 = 0;
                                    for ($p=0; $p<$numdispcar2; $p++) {
                                        for ($q=0; $q<$total_atoms; $q++) {
                                            for ($r=0; $r<3; $r++) {
                                                if ($displacecar2->[$p][$q][$r] != 0) {
                                                    $index_b2 = 0;
                                                    for ($s=0; $s<$numdispcar1; $s++) {
                                                        for ($t=0; $t<$total_atoms; $t++) {
                                                            for ($v=0; $v<3; $v++) {
                                                                if ($displacecar1->[$s][$t][$v] != 0) {
                                                                    if ($t==$q && $v==$r) {
                                                                        $newmatrix->[$index_b1][$index_b2] = $matrix->[$index_a1][$index_a2];
                                                                    }
                                                                    $index_b2++;
                                                                }
                                                            }
                                                        }
                                                    }
                                                    $index_a2++;
                                                }
                                            }
                                        }
                                    }
                                }
                                $index_b1++;
                            }
                        }
                    }
                }
                $index_a1++;
            }
        }
    }
}

for ($i=0; $i<@{$newmatrix}; $i++) {
    for ($j=0; $j<@{$newmatrix->[$i]}; $j++) {
        printf "%9.4f", $newmatrix->[$i][$j]."  ";
    }
    print "\n";
}


