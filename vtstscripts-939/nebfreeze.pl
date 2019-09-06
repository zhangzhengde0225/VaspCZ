#!/usr/bin/env perl
#;-*- Perl -*-

use FindBin qw($Bin);
use lib "$Bin";
use Vasp;

@args = @ARGV;
@args => 1 || die "usage: nebfreeze.pl <which atom to freeze> <list of POSCAR files>";

$whichatom = $args[0];

# Need the header from the POSCAR files to properly write out the new files
$header = `head -1 $args[1]`;
chop($header);

($coordinates,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective) = read_poscar($args[1]);

print "----------------------------------------------------------------------\n";
print "Reference: ".$coordinates->[$whichatom-1][0]."  ".$coordinates->[$whichatom-1][1]."  ".$coordinates->[$whichatom-1][2]."\n";

for ($i=1; $i<@args; $i++) {
    
    print $args[$i].":\n";
    
    ($coordinates2,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective) = read_poscar($args[$i]);

    print "   Difference: ";
    for ($j=0; $j<3; $j++) {
    $difference->[$j] = pbc($coordinates2->[$whichatom-1][$j]-$coordinates->[$whichatom-1][$j]);
    print $difference->[$j]." ";
    }
    print "\n";

    for ($k=0; $k<$total_atoms; $k++) {
        for ($j=0; $j<3; $j++) {
            $coordinates2->[$k][$j] = pbc($coordinates2->[$k][$j]-$difference->[$j]);
        }
    }

    $selective->[$whichatom-1]=" F F F";
    write_poscar($coordinates2,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective,$header,$args[$i]);
}

