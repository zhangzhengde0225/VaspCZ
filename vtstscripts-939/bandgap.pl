#!/usr/bin/env perl
#;-*- Perl -*-

print "\n";

if(@ARGV>0) {
    $arg1=@ARGV[0];
    if ($arg1 =~ m/(cut=)(\d*)(\.*)(\d+)/ ) {
        $cutoff = "$2$3$4";
    }else{
        $cutoff = 0.01;
        print "Unrecognized 'cut=' value.  Cutoff set to 0.01\n";
        print  " \n";
    }
} else {
    $cutoff = 0.01;
    print "Usage 'bandgap.pl cut=__' where __ is the maximum occupancy to consider a band unoccupied. The default is 0.1\n";
    print " \n";
}
    
if ($cutoff == 0.0) {
    print "All bands with nonzero occupancy will be considered occupied\n";
    print " \n";
    $cutoff = 0.000000001;
    # cutoff changed to small value so that 0.0 is less than $cutoff
}else{
    print "Bands with occupancy less than $cutoff will be considered unoccupied\n";
    print " \n";
}

open (OUTCAR, "OUTCAR") || die "Cannot open OUTCAR!";

until (eof OUTCAR) {
    $line = <OUTCAR>;
    if ($line =~ m/ISPIN/) {
        $line =~ m/(\d{1})/;
        $ispin = $1;
    }elsif ($line =~ m/NBANDS/) {
        $line =~ m/(.)(NBANDS=)(\s*)(\d+)(\s*)/;
        $nbands = $4;
        $line =~ m/(NKPTS =)(\s*)(\d+)(\s*)/;
        $kpoints = $3;
    }elsif ($line =~ m/(Iteration)(\s*)(\d+)/) {
        $iteration = $3;
    }
}    

if ($ispin == 2) {
    print "ISPIN = ",$ispin,": Spin-polarized calculation\n";
} elsif ($ispin == 1) {
    print "ISPIN = ",$ispin,": Non-spin calculation\n";
} else {
    print "$ispin\n";
    die "ISPIN incorrect!\n";
}
print "\n";
print "There are $nbands bands and $kpoints k-points\n";
print "The band gap is calculated for the following ionic step: $iteration\n";
print "\n";

if ($ispin == 1 && $cutoff > 2.0) {
    die "WARNING: Bands hold only 2 electrons.  Set cutoff less or equal to than 2.0\n";
} elsif ($ispin == 2 && $cutoff > 1.0) {
    die "WARNING: Bands hold only 1 electron.  Set cutoff less than or equal to 1.0\n"; 
}

open (OUTCAR, "OUTCAR") || die "Cannot open OUTCAR!";

$line = " ";

if ($ispin == 2) {
    
    # Get to the appropriate place in the file
    while ($line !~ m/(Iteration)(\s*)($iteration)/) {
        $line = <OUTCAR>;
    }
    while ($line !~ m/(spin component 1)/) {
        $line = <OUTCAR>;
    }
    <OUTCAR>;
    
    # initialize some of the variables
    $i = 1; # counter
    $spinupHOMO = -100000.;
    $spinupLUMO = 100000.;
    
    # look at each k-point in the spin up portion
    while ($i <= $kpoints) {
        $line = <OUTCAR>;
        
        # Get the positions of kx, ky, and kz
        $line =~ m/(k-point)(\s*)(\d+)(\s*)(:)(\s+)(-*)(\d+)(\.)(\d+)(\s+)(-*)(\d+)(\.)(\d+)(\s+)(-*)(\d+)(\.)(\d+)/;
        $kx = "$7$8$9$10";
        $ky = "$12$13$14$15";
        $kz = "$17$18$19$20";

        # Read over all of the bands until the homo and lumo have been found
        <OUTCAR>;
        undef $lumo;
        $j = 1;
        while ($j <= $nbands) {
            $line = <OUTCAR>;
            $line =~ m/(\d+)(\s+)(-*)(\d+)(\.)(\d{4})(\s+)(-*)(\d+)(\.)(\d{3})/;
            $bandenergy = "$3$4$5$6";
            $occupation = "$8$9$10$11";
            
            # find the HOMO
            if ($occupation < 0.0) {
                die "negative occupation found in band $j of k-point $i\n"
                }elsif ($occupation >= $cutoff) {
                    $homo = $bandenergy;
                }elsif (defined $lumo) {
                        # if the lumo has been defined, keep reading until the end of this k-point
                    }else{
                    $lumo = $bandenergy;    
            }
            $j++;
        }
        
        # Now that the homo and lumo are found, find if they're the highest over all k
        if ($homo > $spinupHOMO) {
            $spinupHOMO = $homo;
            $spinupHOMOkx = $kx;
            $spinupHOMOky = $ky;
            $spinupHOMOkz = $kz;
        }
        if ($lumo < $spinupLUMO) {
            $spinupLUMO = $lumo;
            $spinupLUMOkx = $kx;
            $spinupLUMOky = $ky;
            $spinupLUMOkz = $kz;
        }
        <OUTCAR>;
        $i++;
    }
    
    # Do the same thing for the spin down component
    while ($line !~ m/(spin component 2)/) {
        $line = <OUTCAR>;
    }
    <OUTCAR>;
    
    # initialize some of the variables
    $i = 1; # counter
    $spindownHOMO = -100000.;
    $spindownLUMO = 100000.;
    
    # look at each k-point in the spin up portion
    while ($i <= $kpoints) {
        $line = <OUTCAR>;
        
        # Get the positions of kx, ky, and kz
        $line =~ m/(k-point)(\s*)(\d+)(\s*)(:)(\s+)(-*)(\d+)(\.)(\d+)(\s+)(-*)(\d+)(\.)(\d+)(\s+)(-*)(\d+)(\.)(\d+)/;
        $kx = "$7$8$9$10";
        $ky = "$12$13$14$15";
        $kz = "$17$18$19$20";

        # Read over all of the bands until the homo and lumo have been found
        <OUTCAR>;
        undef $lumo;
        $j = 1;
        while ($j <= $nbands) {
            $line = <OUTCAR>;
            $line =~ m/(\d+)(\s+)(-*)(\d+)(\.)(\d{4})(\s+)(-*)(\d+)(\.)(\d{3})/;
            $bandenergy = "$3$4$5$6";
            $occupation = "$8$9$10$11";
            
            # find the HOMO
            if ($occupation < 0.0) {
                die "negative occupation found in band $j of k-point $i\n"
                }elsif ($occupation >= $cutoff) {
                    $homo = $bandenergy;
                }elsif (defined $lumo) {
                        # if the lumo has been defined, keep reading until the end of this k-point
                    }else{
                    $lumo = $bandenergy;    
            }
            $j++;
        }
        
        # Now that the homo and lumo are found, find if they're the highest over all k
        if ($homo > $spindownHOMO) {
            $spindownHOMO = $homo;
            $spindownHOMOkx = $kx;
            $spindownHOMOky = $ky;
            $spindownHOMOkz = $kz;
        }
        if ($lumo < $spindownLUMO) {
            $spindownLUMO = $lumo;
            $spindownLUMOkx = $kx;
            $spindownLUMOky = $ky;
            $spindownLUMOkz = $kz;
        }
        <OUTCAR>;
        $i++;
    }
    if ($spinupHOMO == $spindownHOMO) {
        print "The top of the valence band is $spinupHOMO eV and it occurs at\n";
        print "the k-point of $spinupHOMOkx $spinupHOMOky $spinupHOMOkz\n";
        print "\n";
        print "The bottom of the conduction band is $spinupLUMO eV and it occurs at\n";
        print "the k-point of $spinupLUMOkx $spinupLUMOky $spinupLUMOkz\n";
        print "\n";
        $gap = $spinupLUMO - $spinupHOMO;
        if ($spinupHOMOkx == $spinupLUMOkx && $spinupHOMOky == $spinupLUMOky && $spinupHOMOkz == $spinupHOMOkz) {
            $direction = "a direct";
        }else{
            $direction = "an indirect";
        }
    }else{
        print "The top of the conduction band differs for majority- and minority-spin states\n";
        print " \n";
        if ($spinupHOMO > $spindownHOMO) {
            $top = $spinupHOMO;
            $topkx = $spinupHOMOkx;
            $topkz = $spinupHOMOkz;
            $topky = $spinupHOMOky;
            $topdir = "majority";
        }else{
            $top = $spindownHOMO;
            $topkx = $spindownHOMOkx;
            $topky = $spindownHOMOky;
            $topkz = $spindownHOMOkz;    
            $topdir = "minority";
        }
        print "The top of the valence band is a $topdir spin state with energy $top eV and it occurs at\n";
        print "the k-point $topkx, $topky, $topkz\n";
        print " \n";
        $gap = $spinupLUMO - $top;
        print "The bottom of the conduction band is $spinupLUMO eV and it occurs at\n";
        print "the k-point of $spinupLUMOkx $spinupLUMOky $spinupLUMOkz in reciprocal coordinates\n";
        print "\n";
        if ($topkx == $spinupLUMOkx && $topky == $spinupLUMOky && $topkz == $spinupHOMOkz) {
            $direction = "a direct";
        }else{
            $direction = "an indirect";
        }
    }

}else{

    # Get to the appropriate place in the file
    while ($line !~ m/(Iteration)(\s*)($iteration)/) {
            $line = <OUTCAR>;
    }
    while ($line !~ m/(E-fermi)/) {
        $line = <OUTCAR>;
    }
    <OUTCAR>;
    <OUTCAR>;
    
    # initialize some of the variables
    $i = 1; # counter
    $spinupHOMO = -100000.;
    $spinupLUMO = 100000.;
    
    # look at each k-point
    while ($i <= $kpoints) {
        $line = <OUTCAR>;
        
        # Get the positions of kx, ky, and kz
        $line =~ m/(k-point)(\s*)(\d+)(\s*)(:)(\s+)(-*)(\d+)(\.)(\d+)(\s+)(-*)(\d+)(\.)(\d+)(\s+)(-*)(\d+)(\.)(\d+)/;
        $kx = "$7$8$9$10";
        $ky = "$12$13$14$15";
        $kz = "$17$18$19$20";
        

        # Read over all of the bands until the homo and lumo have been found
        <OUTCAR>;
        undef $lumo;
        $j = 1;
        while ($j <= $nbands) {
            $line = <OUTCAR>;
            $line =~ m/(\d+)(\s+)(-*)(\d+)(\.)(\d{4})(\s+)(-*)(\d+)(\.)(\d{3})/;
            $bandenergy = "$3$4$5$6";
            $occupation = "$8$9$10$11";
            
            # find the HOMO
            if ($occupation < 0.0) {
                die "negative occupation found in band $j of k-point $i\n"
                }elsif ($occupation >= $cutoff) {
                    $homo = $bandenergy;
                }elsif (defined $lumo) {
                        # if the lumo has been defined, keep reading until the end of this k-point
                    }else{
                    $lumo = $bandenergy;    
            }
            $j++;
        }
        
        # Now that the homo and lumo are found, find if they're the highest over all k
        if ($homo > $spinupHOMO) {
            $spinupHOMO = $homo;
            $spinupHOMOkx = $kx;
            $spinupHOMOky = $ky;
            $spinupHOMOkz = $kz;
        }
        if ($lumo < $spinupLUMO) {
            $spinupLUMO = $lumo;
            $spinupLUMOkx = $kx;
            $spinupLUMOky = $ky;
            $spinupLUMOkz = $kz;
        }
        <OUTCAR>;
        $i++;
    }
    print "The top of the valence band is $spinupHOMO eV and it occurs at\n";
    print "the k-point of $spinupHOMOkx $spinupHOMOky $spinupHOMOkz\n";
    print "\n";
    print "The bottom of the conduction band is $spinupLUMO eV and it occurs at\n";
    print "the k-point of $spinupLUMOkx $spinupLUMOky $spinupLUMOkz in reciprocal coordinates\n";
    print "\n";
    $gap = $spinupLUMO - $spinupHOMO;
    if ($spinupHOMOkx == $spinupLUMOkx && $spinupHOMOky == $spinupLUMOky && $spinupHOMOkz == $spinupHOMOkz) {
        $direction = "a direct";
    }else{
        $direction = "an indirect";
    }
}

$gap =~ m/(\d+)(\.)(\d{3})/;
$gap = "$1$2$3";

print "This produces $direction band gap of $gap eV\n";
print "\n";
