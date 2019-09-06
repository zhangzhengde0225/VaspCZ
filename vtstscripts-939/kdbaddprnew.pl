#!/usr/bin/env perl

#=========================================================================================================

    # kdbaddpr.pl
    # Adds a directory to the 'database' directory tree given a pr000N process directory from an 
    # akmc.pl run.

#=========================================================================================================

#    use Carp 'verbose';
#    $SIG{"INT" } = sub { Carp::confess( @_ ) };

    use strict;
    use FindBin qw($Bin);
    use lib "$Bin";
    use Vasp;
    use kdbutil;
    use File::Copy;
    use Storable qw(dclone);
    use Math::Trig;    
    use Getopt::Long;
    
#=========================================================================================================
# Globals

    # Angstroms an atom must move in order to be included as a mobile atom.
    my $cutoff = 0.7; 
    my $duplicateCutoff = 0.5;
    my $neighborCutoff = 0.3;
    my %radii = %covalentRadii;
    my $twoPi = 6.2831853071795862;


    use constant ANGSTROM_FUDGE          => 0.5;

#=========================================================================================================

    my $debug = '';
    GetOptions("debug"=>\$debug);

    # Store the input filenames.
    my $min1_file_name = $ARGV[0];
    my $min2_file_name = $ARGV[1];
    my $saddle_file_name = $ARGV[2];
    my $mode_file_name = $ARGV[3];

    # Load the POSCAR files for the 2 minima and the saddle.
    my @min1 = read_poscar($min1_file_name);
    my @saddle = read_poscar($saddle_file_name);
    my @min2 = read_poscar($min2_file_name);

    # Load the database location from the environment variable.
    my $dbase = kdbHome();

    # Determine the target directory based on the description.
    my $desc = join("_", sort(split(/ /, @saddle[7])));

    # The total number of atoms.
    my $numAtoms = @saddle[POSCAR_TOTAL_ATOMS];
    
    # ALGO Calculate the distance from each minimum to the saddle and betwixt each minimum, for 
    # ALGO each atom.
    my $min1Min2Diff = pbc_difference($min1[0], $min2[0], $saddle[4]);
    my $min1Min2Car = dirkar($min1Min2Diff, $saddle[1], $saddle[2], $saddle[4]);
    my (@min1SaddleDist, @min2SaddleDist, @min1Min2Dist);
    for(my $i = 0; $i < $numAtoms; $i++)
    {
        $min1Min2Dist[$i] =   sqrt($min1Min2Car->[$i][0] ** 2 + 
                                   $min1Min2Car->[$i][1] ** 2 +
                                   $min1Min2Car->[$i][2] ** 2); 
    }

    # ALGO Make a list of the atoms that move more than $cutoff.
    my %mobileAtoms = ();
    my @mobile = ();
    my $mobileAtomCount = 0;
    my $maxMove = 0.0;
    my $maxMoveAtom = 0;
    for(my $i = 0; $i < $numAtoms; $i++)
    {
        if($min1Min2Dist[$i] >= $cutoff)
        {
            $mobileAtoms{$i} = 1;
            push @mobile, $i;
        }
        if($mobileAtoms{$i} == 1)
        {
            $mobileAtomCount++;
        }
        if($min1Min2Dist[$i] > $maxMove)
        {
            $maxMove = $min1Min2Dist[$i];
            $maxMoveAtom = $i;
        }
    }
    
    # ALGO If no atoms met the cutoff criteria, choose the one that moved the most.
    if($mobileAtomCount == 0)
    {
        $mobileAtoms{$maxMoveAtom} = 1;
        push @mobile, $maxMoveAtom;
        $mobileAtomCount++;
    }
    
    # ALGO Make a list of nearest neighbors to the mobile atoms.
    my @min1Neighbors = ();
    my @min2Neighbors = ();
    my @saddleNeighbors = ();
    # Loop over all atoms.
    for(my $i = 0; $i < $saddle[4]; $i++)
    {
        # If this is a mobile atom, try to find the nearest neightbors.
        if($mobileAtoms{$i})
        {
            # Get the ith atom element.
            my $itype = nthAtomElement($i, \@saddle);
            
            # Loop over all atoms...
            for(my $j = 0; $j < $saddle[4]; $j++)
            {
                # If $j is not a mobile atom, check the distance.
                if(!$mobileAtoms{$j})
                {
                    my $diff = pbc_difference([$min1[0][$i]], [$min1[0][$j]], 1);
                    my $car = dirkar($diff, $saddle[1], $saddle[2], 1);
                    my $min1Dist = sqrt($car->[0][0] ** 2 + $car->[0][1] ** 2 + $car->[0][2] ** 2);
                    $diff = pbc_difference([$min2[0][$i]], [$min2[0][$j]], 1);
                    $car = dirkar($diff, $saddle[1], $saddle[2], 1);
                    my $min2Dist = sqrt($car->[0][0] ** 2 + $car->[0][1] ** 2 + $car->[0][2] ** 2);
                    $diff = pbc_difference([$saddle[0][$i]], [$saddle[0][$j]], 1);
                    $car = dirkar($diff, $saddle[1], $saddle[2], 1);
                    my $saddleDist = sqrt($car->[0][0] ** 2 + $car->[0][1] ** 2 + $car->[0][2] ** 2);
                    # Get the $j element type.
                    my $jtype = nthAtomElement($j, \@saddle);
                    my $maxDist = $radii{$jtype} + $radii{$itype};
                    $maxDist += $maxDist * $neighborCutoff;
                    if($min1Dist < $maxDist || $min2Dist < $maxDist || $saddleDist < $maxDist)
                    {
                        if($min1Dist < $min2Dist && $min1Dist < $saddleDist)
                        {
                            push @min1Neighbors, $j;
                        }
                        elsif($min2Dist < $min1Dist && $min2Dist < $saddleDist)
                        {
                            push @min2Neighbors, $j;
                        }
                        elsif($saddleDist < $min2Dist && $saddleDist < $min1Dist)
                        {
                            push @saddleNeighbors, $j;
                        }
                    }
                }
            }
        }    
    }

    # Create a directory for the process we are about to add.
    if(! -d "$dbase/$desc")
    {
        mkdir("$dbase/$desc");
    }
    
    # ALGO Remove periodic boundary conditions by clumping all the relevant atoms together.
    my @selected = ();
    push @selected, @mobile;
    push @selected, @min1Neighbors;
    push @selected, @min2Neighbors;
    push @selected, @saddleNeighbors;
    # Uniqueify selected.
    my @unique = ();
    my %seen = ();
    foreach my $item(@selected)
    {
        if(!$seen{$item})
        {
            $seen{$item} = 1;
            push @unique, $item;
        }
    }
    @selected = @unique;
    # save a copy of selected for the mode info later...
    my @selectedOriginal = @{dclone(\@selected)};
    
    # Drop unselected atoms.
    my($min1, $mapping) = carReduce(\@min1, \@selected);
    my($min2, $foo1) = carReduce(\@min2, \@selected);
    my($saddle, $foo2) = carReduce(\@saddle, \@selected);
    @min1 = @$min1;
    @min2 = @$min2;
    @saddle = @$saddle;
    for(my $i = 0; $i < @mobile; $i++)
    {
        $mobile[$i] = $mapping->{$mobile[$i]};
    }
    $maxMoveAtom = $mapping->{$maxMoveAtom};
    
    
    my @min1clumped = @{dclone(\@min1)};
    my %seen = ();
    my %pbcd = (0 => 1);
    my $atom = 0;
    
    while(scalar keys %seen < $min1clumped[POSCAR_TOTAL_ATOMS])
    {
        for(my $i = 0; $i < $min1clumped[POSCAR_TOTAL_ATOMS]; $i += 1)
        {
            if($i != $atom && !$pbcd{$i})
            {
                my $diff = pbc_difference([$min1[0][$i]], [$min1[0][$atom]], 1);
                my $car = dirkar($diff, $min1clumped[1], $min1clumped[2], 1);
                my $dist = sqrt($car->[0][0] * $car->[0][0] + $car->[0][1] * $car->[0][1] + $car->[0][2] * $car->[0][2]);
                if($dist < 3.3)
                {
                    $diff = pbc_difference([$min1clumped[0][$i]], [$min1clumped[0][$atom]], 1);
                    $min1clumped[0][$i][0] = $min1clumped[0][$atom][0] + $diff->[0][0];
                    $min1clumped[0][$i][1] = $min1clumped[0][$atom][1] + $diff->[0][1];
                    $min1clumped[0][$i][2] = $min1clumped[0][$atom][2] + $diff->[0][2];
                    $pbcd{$i} = 1;
                }
            }
        }
        $seen{$atom} = 1;
        foreach my $item(keys %pbcd)
        {
            if(!$seen{$item})
            {
                $atom = $item;
                last;
            }
        }
    }
    

    my $diff1s = pbc_difference($saddle[0], $min1[0], $min1[POSCAR_TOTAL_ATOMS]);
    my $diff12 = pbc_difference($min2[0], $min1[0], $min1[POSCAR_TOTAL_ATOMS]);
    my @min2clumped = @{dclone(\@min1clumped)};
    my @saddleClumped = @{dclone(\@min1clumped)};
    for(my $i = 0; $i < $min1clumped[POSCAR_TOTAL_ATOMS]; $i++)
    {
        $min2clumped[POSCAR_COORDINATES][$i][0] += $diff12->[$i][0];
        $min2clumped[POSCAR_COORDINATES][$i][1] += $diff12->[$i][1];
        $min2clumped[POSCAR_COORDINATES][$i][2] += $diff12->[$i][2];
        $saddleClumped[POSCAR_COORDINATES][$i][0] += $diff1s->[$i][0];
        $saddleClumped[POSCAR_COORDINATES][$i][1] += $diff1s->[$i][1];
        $saddleClumped[POSCAR_COORDINATES][$i][2] += $diff1s->[$i][2];
    }
    my @min1clumpedCOM = centerOfMass(\@min1clumped);
    for(my $i = 0; $i < $min1clumped[POSCAR_TOTAL_ATOMS]; $i++)
    {
        $min1clumped[POSCAR_COORDINATES][$i][0] -= $min1clumpedCOM[0];
        $min1clumped[POSCAR_COORDINATES][$i][1] -= $min1clumpedCOM[1];
        $min1clumped[POSCAR_COORDINATES][$i][2] -= $min1clumpedCOM[2];
        $min2clumped[POSCAR_COORDINATES][$i][0] -= $min1clumpedCOM[0];
        $min2clumped[POSCAR_COORDINATES][$i][1] -= $min1clumpedCOM[1];
        $min2clumped[POSCAR_COORDINATES][$i][2] -= $min1clumpedCOM[2];
        $saddleClumped[POSCAR_COORDINATES][$i][0] -= $min1clumpedCOM[0];
        $saddleClumped[POSCAR_COORDINATES][$i][1] -= $min1clumpedCOM[1];
        $saddleClumped[POSCAR_COORDINATES][$i][2] -= $min1clumpedCOM[2];
    }
    # Convert everything to Kartesian coordinates. (db already in Kartesians)
    @min1 = DirKar(\@min1clumped);
    @min2 = DirKar(\@min2clumped);
    @saddle = DirKar(\@saddleClumped);
    

    if($debug)
    {
        writeXYZ(\@saddle, "stripped.xyz");
    }

    # get the saddle center of mass (coords, actually).
    my @saddleCOM = centerOfMass(\@saddle);
    
    # Move the QC saddle and minimums so that the saddle COM is at the origin.
    for(my $a = 0; $a < $saddle[POSCAR_TOTAL_ATOMS]; $a++)
    {
        $saddle[POSCAR_COORDINATES][$a][0] -= $saddleCOM[0];
        $saddle[POSCAR_COORDINATES][$a][1] -= $saddleCOM[1];
        $saddle[POSCAR_COORDINATES][$a][2] -= $saddleCOM[2];
        $min1[POSCAR_COORDINATES][$a][0] -= $saddleCOM[0];
        $min1[POSCAR_COORDINATES][$a][1] -= $saddleCOM[1];
        $min1[POSCAR_COORDINATES][$a][2] -= $saddleCOM[2];
        $min2[POSCAR_COORDINATES][$a][0] -= $saddleCOM[0];
        $min2[POSCAR_COORDINATES][$a][1] -= $saddleCOM[1];
        $min2[POSCAR_COORDINATES][$a][2] -= $saddleCOM[2];
    }
    
    # Find the proper directory number for this process and check for dupes.
    my $diri = 0;
    
    my @saddleElements = CarElementTable(\@saddle);

    my @saddleDistanceTable = atomicDistanceTableKar(\@saddle);
    
    # ALGO Loop over each kdb entry that satisfies the incoming description.
    my $worstDist = 1000.0;
    KDB: while(-d "$dbase/$desc/$diri")
    {
        # Load the database saddle structure.
        my @dbSaddle = loadXYZ2CAR("$dbase/$desc/$diri/saddle.xyz");
        
        # Make sure the number of atoms in each structure is the same.
        if($dbSaddle[POSCAR_TOTAL_ATOMS] != $saddle[POSCAR_TOTAL_ATOMS])
        {
            $diri++;
            next KDB;
        }

        # Make sure the kind and number of each element is equivalent. Assumes the number
        # of each type of atom is the same for the db and incoming cars.

        # Get an element count hash for saddle.
        my %saddleHash = carElementCountHash(\@saddle);
        my %dbSaddleHash = carElementCountHash(\@dbSaddle);
        
        # Make sure each hash key in saddle is copied in dbSaddle.
        while((my $k, my $v) = each %saddleHash)
        {
            if(!$dbSaddleHash{$k})
            {
                $diri++;
                next KDB;
            }
            elsif($dbSaddleHash{$k} != $v)
            {
                $diri++;
                next KDB;
            }
        }

        # Make sure each hash key in dbSaddle is copied in saddle.
        while((my $k, my $v) = each %dbSaddleHash)
        {
            if(!$saddleHash{$k})
            {
                $diri++;
                next KDB;
            }
            elsif($saddleHash{$k} != $v)
            {
                $diri++;
                next KDB;
            }
        }

        my @dbElements = CarElementTable(\@dbSaddle);
        my @elementMatches = ();
        for(my $i = 0; $i < $saddle[POSCAR_TOTAL_ATOMS]; $i++)
        {
            for(my $j = 0; $j < $saddle[POSCAR_TOTAL_ATOMS]; $j++)
            {
                $elementMatches[$i][$j] = matchDescription($saddleElements[$i], $dbElements[$j]);
            }
        }
        
        my @dbSaddleDistanceTable = atomicDistanceTableKar(\@dbSaddle);        
        
        my $mappings = [];
        for(my $i = 0; $i < $saddle[POSCAR_TOTAL_ATOMS]; $i++)
        {
            push(@$mappings, [{0=>$i}, {$i=>0}]);
        }
        
        dieIfMapped($mappings, \@dbSaddleDistanceTable, \@saddleDistanceTable);
        
        $diri++;
    } # End loop over kdbs.

    sub dieIfMapped
    {
        my ($mappings, $dbDist, $carDist) = @_;    
        my $nAtoms = scalar(@$dbDist);
        my $newMapping = [];
        MAPPING: foreach my $mapping(@$mappings)
        {
            my $map = $mapping->[0];
            my $rmap = $mapping->[1];
            my $dba = 0;
            while(exists $map->{$dba})
            {
                $dba++;
            }
            if($dba >= $nAtoms)
            {
                print "duplicate of $dbase/$desc/$diri\n";
                exit;
            }
            CA: for(my $ca = 0; $ca < $nAtoms; $ca++)
            {
                if($rmap->{$ca})
                {
                    next CA;    
                }
                DBM: foreach my $dbm(keys %$map)
                {
                    if($dbm == $dba)
                    {
                        next DBM;
                    }
                    if(abs($dbDist->[$dba][$dbm] - $carDist->[$ca][$map->{$dbm}]) > ANGSTROM_FUDGE)
                    {   
                        next CA;
                    }
                }
                my %newMap = %$map;
                my %newRmap = %$rmap;
                $newMap{$dba} = $ca;
                $newRmap{$ca} = $dba;
                dieIfMapped([[\%newMap, \%newRmap]], $dbDist, $carDist);
            }
        }
    }

    # If everything is okay, add the process to the database.
    # Create the directory.
    my $home = "$dbase/$desc/$diri";
    mkdir("$home", 0777);
    
    # Save the original files.
    copy($saddle_file_name, "$home/.saddle.car");
    copy($min1_file_name, "$home/.min1.car");
    copy($min2_file_name, "$home/.min2.car");
    copy($mode_file_name, "$home/.mode.car");
    
    # -- Save xyz files ---
    writeXYZ(\@min1, "$home/min1.xyz");
    writeXYZ(\@min2, "$home/min2.xyz");
    writeXYZ(\@saddle, "$home/saddle.xyz");
    
    # -- Save the MODECAR info. ---
    # Load the MODECAR and the kdb ouput file
    my @mode = read_othercar($mode_file_name);
    open (MODE, ">$home/mode");
    
    foreach my $selection (sort { $a <=> $b } @selectedOriginal)
    {
        print MODE sprintf("% .16e   % .16e   % .16e\n", ($mode[0][$selection][0], 
                                                          $mode[0][$selection][1], 
                                                          $mode[0][$selection][2]));
    }
    close MODE;

#    # Touch the mirror file if the mirror flag is set.
#    if($mirror_flag == 1)
#    {
#        system("touch $dbase/$desc/$diri/mirror");
#    }
    
    # Save the original process directory.
    open (INFO, ">$home/.info");
    print INFO "kdb directory: $home\n";
    close INFO;
    open (MOBILE, ">$home/mobile");
    for(my $mob = 0; $mob < @mobile; $mob++)
    {
        print MOBILE "$mobile[$mob]\n";
    }
    close MOBILE;

    # Finished with this process.
    printf("good.\n");


#===================================================================================================
# Subroutine carReduce($car, $atoms)
#
# Takes a vasp car structure and gets rid of any atoms not listed in $selected.  Returns the new, 
# reduced vasp car structure and a mapping of $selected to the new atoms.
#---------------------------------------------------------------------------------------------------
    sub carReduce
    {
        my ($old, $selected) = @_;
        my %atoms = map {$_ => 1} @$selected;
        my $mapping = {};
        my $new = dclone($old);
        $new->[POSCAR_COORDINATES] = ();
        $new->[POSCAR_NUM_ATOMS] = ();
        $new->[POSCAR_TOTAL_ATOMS] = scalar(@$selected);
        $new->[POSCAR_SELECTIVE] = [];
        $new->[POSCAR_DESCRIPTION] = "";
        my @splitDesc = split(" ", $old->[POSCAR_DESCRIPTION]);
        my $j = 0;
        my $comp = 0;
        my $compCount = 0;
        for(my $i = 0; $i < $old->[POSCAR_TOTAL_ATOMS]; $i++)
        {
            if($atoms{$i})
            {
                $new->[POSCAR_COORDINATES][$j] = ();
                $new->[POSCAR_COORDINATES][$j][0] = $old->[POSCAR_COORDINATES][$i][0];
                $new->[POSCAR_COORDINATES][$j][1] = $old->[POSCAR_COORDINATES][$i][1];
                $new->[POSCAR_COORDINATES][$j][2] = $old->[POSCAR_COORDINATES][$i][2];
                $new->[POSCAR_NUM_ATOMS][$comp] += 1;
                $new->[POSCAR_SELECTIVE]->[$j] = $old->[POSCAR_SELECTIVE]->[$i];
                if(index($new->[POSCAR_DESCRIPTION], $splitDesc[$comp]) == -1)
                {
                    $new->[POSCAR_DESCRIPTION] .= $splitDesc[$comp] . " ";
                }
                $mapping->{$i} = $j;
                $j++;
            }
            $compCount++;
            if($old->[POSCAR_NUM_ATOMS] == $compCount)
            {
                $comp++;
                $compCount = 0;
            }
        }
        return $new, $mapping;
    }














#    # See if the two processes are mirror images, and if so, flag it to use only one process in the query.
#    my $mirror_flag = 0;
#    my $mag1 = magnitude($min1SaddleDiff, $min1[4]);
#    my $mag2 = magnitude($min2SaddleDiff, $min2[4]);
#    if(abs($mag1 - $mag2) < 0.1) #TODO: Parameterize this cutoff.
#    {
#        my $unit = unit($min1Min2Diff, $min2[4]);
#        my $dot = dot_product($min2SaddleDiff, $unit, $min2[4]);
#        my $newvec = vsum(vmult($unit, $dot * -2.0, $min2[4]), $min2SaddleDiff, $min2[4]);
#        my $diff = magnitude(pbc_difference($newvec, $min1SaddleDiff, $min1[4]), $min1[4]);
#        if($diff < 0.01)
#        {
#            $mirror_flag = 1;
#        }
#    }
