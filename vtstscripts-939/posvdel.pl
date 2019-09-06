#!/usr/bin/env perl
#;-*- Perl -*-

@ARGV <= 1 || die "usage: posvstrip.pl <POSCAR>\n";

$poscar = "POSCAR";
if(@ARGV>0) { $poscar = $ARGV[0]; }
$posout = $poscar."_nov";

open(IN,"$poscar") || die "cannot open $poscar for input \n";
open(OUT,">$posout") || die "cannot open $posout for output \n";

$header = "";
$image = 0;

while(1){

    $image++;

    # header of 6 lines (vasp5 format only)
    if($header eq "") {
        for($i=0; $i<6; $i++) {
            $line = <IN>;
            print OUT $line;
        }
    } else {
        print OUT $header;
    }

    # add up the number of atoms 
    $line = <IN>;
    print OUT $line;
    chomp($line);
    $line =~ s/^\s+//g;
    @natomsa = split(/\s+/,$line);
    $nelements = @natomsa;
    $natoms = 0;
    for($i=0; $i<$nelements; $i++) {
        $natoms += $natomsa[$i]; }
    $line = <IN>;
    print OUT $line;

    # print out all the atom positions
    for($i=0; $i<$natoms; $i++) {
        $line = <IN>;
        print OUT $line;
    }

    # try to detect if there is velocity data or another POSCAR file
    $header = "";
    if(!($line1 = <IN>)) { # quit if this is the end of the file
        print "image $image velocity data: none\n";
        last; 
    }
    $header .= $line1;
    chomp($line1);
    $line1 =~ s/^\s+//g;
    $line2 = <IN>;
    $header .= $line2;
    chomp($line2);
    $line2 =~ s/^\s+//g;
    @line2a = split(/\s+/,$line2);
    if (($line1 eq "") && (scalar(@line2a) == 3)) {
        # detected velocities
        print "image $image velocity data: deleting\n";
        for($i=1; $i<$natoms; $i++) { $line=<IN>; }
        $header = "";
    } else {
        # no velocity data
        print "image $image velocity data: none\n";
        for($i=2; $i<6; $i++) {
            $header .= <IN>;
        }
    }
}

system("mv $posout $poscar");
