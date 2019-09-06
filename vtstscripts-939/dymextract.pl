#!/usr/bin/env perl
#;-*- Perl -*-

use FindBin qw($Bin);
use lib "$Bin";
use Vasp;

# This program extracts the dynamical matrix information (the forces vs displacement)
# GH: I'm actually not sure why we have this script; consider replacing or removing it

@args = @ARGV;
@args == 3 || die "usage: extract_matrix.pl <DISPLACECAR 1> <DISPLACECAR 2> <matrix for 1>\n";

#----------------------------------------------------------------------
# read displacecars, count displacements
#----------------------------------------------------------------------
for ($m=0; $m<2; $m++) {
    $filename = $args[$m];
    ($displacecar->[$m],$total_atoms)
    = read_othercar($filename);
}  

#----------------------------------------------------------------------
# read matrix
#----------------------------------------------------------------------
$i = 0;
open (IN,$args[2]);
while (<IN>) {
    $line = $_;
    $line =~ s/^\s+//;
    @line = split(/\s+/,$line);
    for ($j=0; $j<@line; $j++) {
        $matrix->[$i][$j] = $line[$j];
    }
    $i++;
}
close(IN);

#$dof=@{$matrix};
#print $dof."\n";

#______________________________________________________________________
# check which forces correspond to displacements, build matrix
#----------------------------------------------------------------------
$current_displacement = 0;
for ($i=0; $i<$total_atoms; $i++) {
    for ($j=0; $j<3; $j++) {
        if ($displacecar->[0][$i][$j] != 0) {
            if ($displacecar->[1][$i][$j] == 0) {
                for ($k=0;$k<@{$matrix};$k++) {
                    $matrix->[$current_displacement][$k] = 9999;
                    $matrix->[$k][$current_displacement] = 9999;
                }
                # print "Line $current_displacement erased\n";
            }
            $current_displacement++;
        }
   }
}

for ($i=0; $i<@{$matrix}; $i++) {
    $return = 0;
    for ($j=0; $j<@{$matrix->[$i]}; $j++) {
        if ($matrix->[$i][$j] != 9999) {
            printf "%9.4f", $matrix->[$i][$j]."  ";
            $return = 1;
        }
    }
    if ($return == 1) {
        print "\n";
    }
}

