#!/usr/bin/env perl
#;-*- Perl -*-

# This program takes two DISPLACECAR files and creates a third which contains the 
# degrees of freedom that aren't common to both.

use FindBin qw($Bin);
use lib "$Bin";
use Vasp;

# variables needed by the script

@args = @ARGV;
@args == 2 || die "usage: dymcmpdisp.pl <DISPLACECAR 1 filename> <DISPLACECAR 2 filename>";

$displacefile1 = $args[0];
$displacefile2 = $args[1];

($coordinates1,$total_atoms)
  = read_othercar($displacefile1);

($coordinates2,$total_atoms)
  = read_othercar($displacefile2);

for ($i=0; $i<$total_atoms; $i++) {
    for ($j=0; $j<3; $j++) {
        if ($coordinates1->[$i][$j] != $coordinates2->[$i][$j]) {
            $new_displace = $coordinates1->[$i][$j];
            if ($new_displace == 0) {
                $new_displace = $coordinates2->[$i][$j];
            }
        } else {
            $new_displace = 0;
        }
        print $new_displace." ";
    }
    print ($i+1);
    print "\n";
}
