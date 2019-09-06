#!/usr/bin/env perl
#;-*- Perl -*-

use FindBin qw($Bin);
use lib "$Bin";
use Vasp;

$blackhole = "/dev/null";

# Count the number of images present 
$im = 1;
$dir = set_dir($im);
while(-e $dir) {
    $im++;
    $dir = set_dir($im);
}
$nim = --$im;
$nn = 2*$nim-1;
print "There are $nim images -> will be $nn \n";

# Save the old images and generate new ones
for($im = 1 ; $im <= $nim ; $im++) {
    $dir = set_dir($im);
`   mv $dir $dir"_prev"`;
}
for($im = 1 ; $im <= $nn ; $im++) {
    $dir = set_dir($im);
    mkdir $dir;
}

# Move the original images to their new homes
for($im = 1 ; $im <= $nim ; $im++) {
    $i = 2*$im-1;
    $new = set_dir($i);
    $old = set_dir($im);
    $old = $old."_prev";
    `cp $old/POSCAR $new`;
    `cp $old/MODECAR $new`;
}

# Interpolate the new images
for($im = 2; $im < $nn; $im = $im+2) {
    $dir = set_dir($im);
    $tail = set_dir($im-1);
    $head = set_dir($im+1);
    `$bin/posinterp.pl $tail/POSCAR $head/POSCAR 0.5`;
    `mv POSCAR.out $dir/POSCAR`;
}

# Make a new MODECAR file, where we interpolate between the old images
($modecar,$nl) = read_othercar("01/MODECAR");
$natoms = $nl/(2*$nim);
undef $new;
# -> First image
for($i = 0; $i < $natoms; $i++) {
    for($j = 0; $j < 3; $j++) {
        $new->[$i][$j] = $modecar->[$i][$j];
    }
}
# -> Middle images
for($im = 2; $im < $nn; $im++) {
    $a = ($im-1)*$natoms;
    if(!($im%2)) {
        $b = ($im/2-1)*$natoms;
        $c = ($im/2)*$natoms;
        for($i = 0; $i < $natoms; $i++) {
            for($j = 0; $j < 3; $j++) {
                $new->[$a+$i][$j] = 0.5*($modecar->[$b+$i][$j] + $modecar->[$c+$i][$j]);
            }
        }
    } else {
        $b = (($im+1)/2-1)*$natoms;
        for($i = 0; $i < $natoms; $i++) {
            for($j = 0; $j < 3; $j++) {
                $new->[$a+$i][$j] = $modecar->[$b+$i][$j];
            }
        }
    }
}
# -> Last image
$a = ($nn-1)*$natoms;
$b = (($nn+1)/2-1)*$natoms;
for($i = 0 ; $i < $natoms ; $i++) {
    for($j = 0 ; $j < 3 ; $j++) {
        $new->[$a+$i][$j] = $modecar->[$b+$i][$j];
    }
}

# Re-normalize the mode, multiply with sqrt(2) since this is only half the elements 
$a = $nn*$natoms;
$norm = sqrt(2)*magnitude($new,$a);
for($i = 0; $i < $a; $i++) {
    for($j = 0; $j < 3; $j++) {
        $new->[$i][$j] /= $norm;
    }
}

# Right out the whole mode
if(-e "NEWMODECAR") {
    unlink "NEWMODECAR";
}
open NEW , ">NEWMODECAR";
for($i = 0; $i < $nn*$natoms; $i++) {
    printf NEW "%20.10e %20.10e %20.10e \n", $new->[$i][0],$new->[$i][1],$new->[$i][2];
}
$a = $nn*$natoms;
for($im = 0; $im < $nn; $im++) {
    $jm = $nn - $im;
    $a = ($jm-1)*$natoms;
    for($i = 0; $i < $natoms; $i++) {
      printf NEW "%20.10e %20.10e %20.10e \n", $new->[$a+$i][0],$new->[$a+$i][1],$new->[$a+$i][2];
    }
}
close NEW;

for($im = 1; $im <= $nn; $im++) {
    $dir = set_dir($im);
    `cp NEWMODECAR $dir/MODECAR`;
}

# -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.- #

sub set_dir {
    my $num = shift;
    my $dir;
    if($num < 10) {
        $dir = "0$num";
    } elsif($num < 100) {
        $dir = "$num";
    }
    return($dir);
}
