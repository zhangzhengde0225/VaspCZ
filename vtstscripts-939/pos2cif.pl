#! /usr/bin/env perl

use FindBin qw($Bin);
use lib "$Bin";
use Vasp;
use Math::Trig;
$fact = 180/pi;

@args = @ARGV;
@args >= 1 || die "usage: pos2cif.pl <POSCAR>\n";
$inputfilename = $args[0];
$filetype = "vasp";

$inputfile = "";
open (IN,$inputfilename);
while (<IN>) {
    $_ =~ s/^\s+//g;
    $inputfile .= $_;
}
close (IN);

@inputfile = split(/\n/,$inputfile);
$header = $inputfile[0];

if($filetype eq "vasp") {
    print "\nConverting vasp to cif\n";
    @header = split(/\s+/,$header);

    ($coordinates,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective,$description)
    = read_poscar($inputfilename);
    set_bc($coordinates,$total_atoms);

    if(@args >= 2) {
        $outputfilename = $args[1];
    } else {
        $outputfilename = $inputfilename.".cif";
    }

    open (OUT,">$outputfilename");
    $description =~ s/^\s+//;
    @atom_types = split(/\s+/,$description);
    $description =~ s/ //g;
    print OUT "data_".$description."\n";
    print OUT "_audit_creation_method           'pos2cif.pl'"."\n";

    for($j=0; $j<3; $j++) {
        $vector1->[0][$j] = $basis->[$j][0];
        $vector2->[0][$j] = $basis->[$j][1];
        $vector3->[0][$j] = $basis->[$j][2];
    }
    $mag1 = magnitude($vector1,1);
    $mag2 = magnitude($vector2,1);
    $mag3 = magnitude($vector3,1);

    $angle1 = acos(dot_product($vector2,$vector3,1)/$mag2/$mag3)*$fact;
    $angle2 = acos(dot_product($vector1,$vector3,1)/$mag1/$mag3)*$fact;
    $angle3 = acos(dot_product($vector1,$vector2,1)/$mag1/$mag2)*$fact;

print OUT ("_cell_length_a               ". $mag1."\n"."_cell_length_b               ". $mag2."\n"."_cell_length_c               ". $mag3."\n");
print OUT ("_cell_angle_alpha            ". $angle1."\n"."_cell_angle_beta             ". $angle2."\n"."_cell_angle_gamma            ". $angle3);
print OUT q{
_symmetry_space_group_H-M        'P1'
_symmetry_Int_Tables_number      '1'
_symmetry_cell_setting           'triclinic'
loop_
_symmetry_equiv_pos_as_xyz
x,y,z
};
print OUT q{
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_occupancy
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_U_iso_or_equiv
};

    $temp = @{$num_atoms};
    $count = 0;
    for ($i=0; $i<@atom_types; $i++) {
        for($j=0; $j<$num_atoms->[$i]; $j++) {
            print OUT $atom_types[$i],$j+1;
            print OUT "  ".$atom_types[$i]."  "."1.0000"."  ";
            for ($k = 0; $k < 3; $k++) {
                printf OUT "%5.5f", $coordinates->[$j+$count][$k];
                print OUT "  ";
            }
            print OUT "0.0000"."\n";
        }
        $count += $num_atoms->[$i];
    }
}
