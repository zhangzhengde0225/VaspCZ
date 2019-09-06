#!/usr/bin/env perl
#;-*- Perl -*-
################################################################################
$var[0] = 0;
$var[1] = 0;
$var[2] = 0;
@args = @ARGV;
$inputfilename = $args[0];
if (@args==4) {
    $var[0] = $args[1];
    $var[1] = $args[2];
    $var[2] = $args[3];
}
print "shift: @var \n";
if ($inputfilename eq "") {
    $inputfilename = "POSCAR";
}
$inputfile = "";
open (IN,$inputfilename);
while (<IN>) {
    $_ =~ s/^\s+//g;
    $inputfile .= $_;
}
close (IN);
@inputfile = split(/\n/,$inputfile);
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
$descript = $inputfile[0];
$scale = $inputfile[1];
$basis[0] = $inputfile[2];
$basis[1] = $inputfile[3];
$basis[2] = $inputfile[4];
$natoms = $inputfile[5];
@natoms = split(/\s+/,$natoms);
$ntypes = @natoms;
$totatoms = 0;
for ($i=0; $i<@natoms; $i++) {
    $totatoms += $natoms[$i];
}
$type = $inputfile[6];
$s = 7;
if (substr($type,0,1) eq "S") {
    $s = 8;
    $type = $inputfile[7];
}
@coords = @inputfile[$s..$s+$totatoms-1];

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
################################################################################
#   print POSCAR
################################################################################

open(VASP,">POSCAR.out") || die "error writing output\n";
print VASP "$descript\n";
printf VASP " %.17f\n", $scale;
for ($i=0; $i<3; $i++) {
    $_ = $basis[$i];
    chop();
    @temp = split();
    printf VASP " %22.16f%22.16f%22.16f\n",@temp;
}
print VASP "@natoms\n";
print VASP "Selective dynamics\n";
print VASP "Direct\n";
for ($i=0; $i<$totatoms; $i++) {
    $_ = $coords[$i];
    chop();
    (@temp) = split();
    @tum = @temp[0..2];
    for ($nie=0; $nie<3; $nie++) {
        $tum[$nie] = $tum[$nie] + $var[$nie];
        if ($tum[$nie]>1) {
            $tum[$nie]--;
        }
        if ($tum[$nie]<0) {
            $tum[$nie]++;
        }
    }
    printf VASP "%20.16f%20.16f%20.16f   T   T   T\n", @tum;
}
close(VASP);

