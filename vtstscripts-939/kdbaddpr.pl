#!/usr/bin/env perl

#=========================================================================================================

    # kdbaddpr.pl
    # Adds a directory to the 'database' directory tree given a pr000N process directory from an 
    # akmc.pl run.

#=========================================================================================================

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
    my $springK = 10.0;
    my $maxIter = 1000;

    # The indices of the read_poscar data format.
    use constant POSCAR_COORDINATES      => 0;
    use constant POSCAR_BASIS            => 1;
    use constant POSCAR_LATTICE          => 2;
    use constant POSCAR_NUM_ATOMS        => 3;
    use constant POSCAR_TOTAL_ATOMS      => 4;
    use constant POSCAR_SELECTIVE_FLAG   => 5;
    use constant POSCAR_SELECTIVE        => 6;
    use constant POSCAR_DESCRIPTION      => 7;
    use constant XYZ_POSCAR_SWITCH       => 8;

    use constant ANGSTROM_FUDGE          => 0.5;

#=========================================================================================================

    my $debug = '';
    GetOptions("debug"=>\$debug);

    # Failure flag.
    my $goodToGo = 1;

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
    my $numAtoms = @saddle[4];
    
    # ALGO Calculate the distance from each minimum to the saddle and betwixt each minimum, for 
    # ALGO each atom.
    my $min1SaddleDiff = pbc_difference($min1[0], $saddle[0], $saddle[4]);
    my $min1SaddleCar = dirkar($min1SaddleDiff, $saddle[1], $saddle[2], $saddle[4]);
    my $min2SaddleDiff = pbc_difference($min2[0], $saddle[0], $saddle[4]);
    my $min2SaddleCar = dirkar($min2SaddleDiff, $saddle[1], $saddle[2], $saddle[4]);
    my $min1Min2Diff = pbc_difference($min1[0], $min2[0], $saddle[4]);
    my $min1Min2Car = dirkar($min1Min2Diff, $saddle[1], $saddle[2], $saddle[4]);
    my (@min1SaddleDist, @min2SaddleDist, @min1Min2Dist);
    for(my $i = 0; $i < $numAtoms; $i++)
    {
        $min1SaddleDist[$i] = sqrt($min1SaddleCar->[$i][0] ** 2 + 
                                   $min1SaddleCar->[$i][1] ** 2 +
                                   $min1SaddleCar->[$i][2] ** 2); 
        $min2SaddleDist[$i] = sqrt($min2SaddleCar->[$i][0] ** 2 + 
                                   $min2SaddleCar->[$i][1] ** 2 +
                                   $min2SaddleCar->[$i][2] ** 2); 
        $min1Min2Dist[$i] =   sqrt($min1Min2Car->[$i][0] ** 2 + 
                                   $min1Min2Car->[$i][1] ** 2 +
                                   $min1Min2Car->[$i][2] ** 2); 
    }

    # See if the two processes are mirror images, and if so, flag it to use only one process in the query.
    my $mirror_flag = 0;
    my $mag1 = magnitude($min1SaddleDiff, $min1[4]);
    my $mag2 = magnitude($min2SaddleDiff, $min2[4]);
    if(abs($mag1 - $mag2) < 0.1) #TODO: Parameterize this cutoff.
    {
        my $unit = unit($min1Min2Diff, $min2[4]);
        my $dot = dot_product($min2SaddleDiff, $unit, $min2[4]);
        my $newvec = vsum(vmult($unit, $dot * -2.0, $min2[4]), $min2SaddleDiff, $min2[4]);
        my $diff = magnitude(pbc_difference($newvec, $min1SaddleDiff, $min1[4]), $min1[4]);
        if($diff < 0.01)
        {
            $mirror_flag = 1;
        }
    }
    

    # ALGO Make a list of the atoms that move more than $cutoff.
    my %mobileAtoms = ();
    my @mobile = ();
    my $mobileAtomCount = 0;
    my $maxMove = 0.0;
    my $maxMoveAtom = 0;
    for(my $i = 0; $i < $numAtoms; $i++)
    {
#        if($min1SaddleDist[$i] >= $cutoff)
#        {
#            $mobileAtoms{$i} = 1;
#            push @mobile, $i;
#        }
#        elsif($min2SaddleDist[$i] >= $cutoff)
#        {
#            $mobileAtoms{$i} = 1;
#            push @mobile, $i;
#        }
#        elsif($min1Min2Dist[$i] >= $cutoff)
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
                    my $min1Dist = sqrt($car->[0][0] ** 2 + $car->[0][1] ** 2 + 
                                        $car->[0][2] ** 2);
                    $diff = pbc_difference([$min2[0][$i]], [$min2[0][$j]], 1);
                    $car = dirkar($diff, $saddle[1], $saddle[2], 1);
                    my $min2Dist = sqrt($car->[0][0] ** 2 + $car->[0][1] ** 2 + 
                                        $car->[0][2] ** 2);
                    $diff = pbc_difference([$saddle[0][$i]], [$saddle[0][$j]], 1);
                    $car = dirkar($diff, $saddle[1], $saddle[2], 1);
                    my $saddleDist = sqrt($car->[0][0] ** 2 + $car->[0][1] ** 2 + 
                                          $car->[0][2] ** 2);
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
    my $i = 0;
    while($i < $min1[4])
    {
        my %is_selected;
        for(@selected) {$is_selected{$_} = 1;}
        if(!$is_selected{$i})
        {
            @min1 = carDropAtom($i, \@min1);
            @min2 = carDropAtom($i, \@min2);
            @saddle = carDropAtom($i, \@saddle);
            for(my $j = 0; $j < @selected; $j++)
            {
                if($selected[$j] >= $i)
                {
                    $selected[$j]--;
                }
            }
            for(my $j = 0; $j < @mobile; $j++)
            {
                if($mobile[$j] >= $i)
                {
                    $mobile[$j]--;
                }
            }
            if($maxMoveAtom >= $i)
            {
                $maxMoveAtom--;
            }
        }
        else
        {
            $i++;
        }
    }            

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
                if($dist < 3.0)
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
    
    if($debug)
    {
        my @temp = centerOfMass(\@saddle);
        print "    Incoming saddle center of mass: ($temp[0], $temp[1], $temp[2])\n";
    }
    
    # Find the proper directory number for this process and check for dupes.
    my $diri = 0;
    my $passedFilter = 1;
    
    # ALGO Find the two atoms farthest from the origin.
    # Find the farthest.
    my $saddleFarthestAtom = -1;
    my $saddleFarthestDist = 0.0;
    for(my $a = 0; $a < @saddle[POSCAR_TOTAL_ATOMS]; $a++)
    {
        my @dist = ($saddle[POSCAR_COORDINATES][$a][0],
                    $saddle[POSCAR_COORDINATES][$a][1],
                    $saddle[POSCAR_COORDINATES][$a][2]);
        my $dist = sqrt($dist[0] * $dist[0] + $dist[1] * $dist[1] + $dist[2] * $dist[2]);
        if($dist > $saddleFarthestDist)
        {
            $saddleFarthestDist = $dist;
            $saddleFarthestAtom = $a;
        }                    
    }

    
    # Find the second farthest.
    my $saddleSecondFarthestAtom = -1;
    my $saddleSecondFarthestDist = 0.0;
    for(my $a = 0; $a < @saddle[POSCAR_TOTAL_ATOMS]; $a++)
    {
        my @dist = ($saddle[POSCAR_COORDINATES][$a][0],
                    $saddle[POSCAR_COORDINATES][$a][1],
                    $saddle[POSCAR_COORDINATES][$a][2]);
        my $dist = sqrt($dist[0] * $dist[0] + $dist[1] * $dist[1] + $dist[2] * $dist[2]) ;
        if($dist > $saddleSecondFarthestDist && $a != $saddleFarthestAtom)
        {
            $saddleSecondFarthestDist = $dist;
            $saddleSecondFarthestAtom = $a;
        }                    
    }
    

    my @dist = ($saddle[POSCAR_COORDINATES][$saddleFarthestAtom][0] - $saddle[POSCAR_COORDINATES][$saddleSecondFarthestAtom][0], 
                $saddle[POSCAR_COORDINATES][$saddleFarthestAtom][1] - $saddle[POSCAR_COORDINATES][$saddleSecondFarthestAtom][1], 
                $saddle[POSCAR_COORDINATES][$saddleFarthestAtom][2] - $saddle[POSCAR_COORDINATES][$saddleSecondFarthestAtom][2]);
    my $farthestSecondFarthestDist = sqrt($dist[0] * $dist[0] + $dist[1] * $dist[1] + $dist[2] * $dist[2]);
    
    my @saddleElements = CarElementTable(\@saddle);
    
    my @farthestSaddleUnit = ($saddle[POSCAR_COORDINATES][$saddleFarthestAtom][0],
                              $saddle[POSCAR_COORDINATES][$saddleFarthestAtom][1],
                              $saddle[POSCAR_COORDINATES][$saddleFarthestAtom][2]);
    @farthestSaddleUnit = unitVector(\@farthestSaddleUnit);
    my @secondFarthestSaddleUnit = ($saddle[POSCAR_COORDINATES][$saddleSecondFarthestAtom][0],
                                   $saddle[POSCAR_COORDINATES][$saddleSecondFarthestAtom][1],
                                   $saddle[POSCAR_COORDINATES][$saddleSecondFarthestAtom][2]);
    @secondFarthestSaddleUnit = unitVector(\@secondFarthestSaddleUnit);
    
    
    # ALGO Loop over each kdb entry that satisfies the incoming description.
    my $worstDist = 1000.0;
    while(-d "$dbase/$desc/$diri" && $goodToGo == 1)
    {
        # Load the database structures.
        my @dbMin1 = loadXYZ2CAR("$dbase/$desc/$diri/min1.xyz");
        my @dbMin2 = loadXYZ2CAR("$dbase/$desc/$diri/min2.xyz");
        my @dbSaddle = loadXYZ2CAR("$dbase/$desc/$diri/saddle.xyz");

        
        # Make sure the kind and number of each element is equivalent. Assumes the number
        # of each type of atom is the same for min1 min2 and saddle for the db and incoming
        # cars.
        $passedFilter = 1;
        if($passedFilter == 1)
        {
            # Get an element count hash for each min1.
            my %min1Hash = carElementCountHash(\@min1);
            my %dbMin1Hash = carElementCountHash(\@dbMin1);
            # Make sure each hash key in min1 is copied in dbMin1.
            while((my $k, my $v) = each %min1Hash)
            {
                if(!$dbMin1Hash{$k})
                {
                    $passedFilter = 0;
                }
                elsif($dbMin1Hash{$k} != $v)
                {
                    $passedFilter = 0;
                }
            }
            # Make sure each hash key in dbMin1 is copied in min1.
            while((my $k, my $v) = each %dbMin1Hash)
            {
                if(!$min1Hash{$k})
                {
                    $passedFilter = 0;
                }
                elsif($min1Hash{$k} != $v)
                {
                    $passedFilter = 0;
                }
            }
        }                

        if($debug)
        {
            if(!$passedFilter)
            {
                print "    filter: number and type of elements failed. ($diri)\n";
            }
            else
            {
                print "    filter: number and type of elements passed. ($diri)\n";
            }
        }
        

        # ALGO If all prior filters have passed, check to see if the configurations can be
        # ALGO matched through rotation.
        if($passedFilter == 1)
        {

            # Get the center of masses for the database structures.
            my @dbSaddleCOM = centerOfMass(\@dbSaddle);

            # Move the db center of mass to the origin.
            for(my $a = 0; $a < $saddle[POSCAR_TOTAL_ATOMS]; $a++)
            {
                $dbSaddle[POSCAR_COORDINATES][$a][0] -= $dbSaddleCOM[0];
                $dbSaddle[POSCAR_COORDINATES][$a][1] -= $dbSaddleCOM[1];
                $dbSaddle[POSCAR_COORDINATES][$a][2] -= $dbSaddleCOM[2];
            }

            if($debug)
            {
                my @temp = centerOfMass(\@dbSaddle);
                print "    DB new center of mass:          ($temp[0], $temp[1], $temp[2])\n";
            }

            # Create a list of dbSaddle atoms outside $saddleSecondFarthestDist - 0.25 of the origin.
            my @dbPairAtoms = ();
            for(my $i = 0; $i < $dbSaddle[POSCAR_TOTAL_ATOMS]; $i++)
            {
                my @dist = ($dbSaddle[POSCAR_COORDINATES][$i][0],
                            $dbSaddle[POSCAR_COORDINATES][$i][1],
                            $dbSaddle[POSCAR_COORDINATES][$i][2]);
                my $dist = sqrt($dist[0] * $dist[0] + $dist[1] * $dist[1] + $dist[2] * $dist[2]);
                if($dist > $saddleSecondFarthestDist - 0.25)
                {
                    push(@dbPairAtoms, $i)
                }
            }

            if($debug)
            {
                print "    pair atoms: ";
                for(my $q = 0; $q < @dbPairAtoms; $q++)
                {                 
                    print "$dbPairAtoms[$q] ";
                }
                print "\n";
            }

            # Get a list of viable pairs of atoms.
            my @dbAtomPairs = ();
            for(my $i = 0; $i < @dbPairAtoms; $i++)
            {
                my @dist = ($dbSaddle[POSCAR_COORDINATES][$dbPairAtoms[$i]][0],
                            $dbSaddle[POSCAR_COORDINATES][$dbPairAtoms[$i]][1],
                            $dbSaddle[POSCAR_COORDINATES][$dbPairAtoms[$i]][2]);
                my $dist = sqrt($dist[0] * $dist[0] + $dist[1] * $dist[1] + $dist[2] * $dist[2]);
                if(abs($dist - $saddleFarthestDist) < ANGSTROM_FUDGE)
                {
                    for(my $j = 0; $j < @dbPairAtoms; $j++)
                    {
                        if($i != $j)
                        {
                            my @dist = ($dbSaddle[POSCAR_COORDINATES][$dbPairAtoms[$j]][0],
                                        $dbSaddle[POSCAR_COORDINATES][$dbPairAtoms[$j]][1],
                                        $dbSaddle[POSCAR_COORDINATES][$dbPairAtoms[$j]][2]);
                            my $dist = sqrt($dist[0] * $dist[0] + $dist[1] * $dist[1] + $dist[2] * $dist[2]);
                            if(abs($dist - $saddleSecondFarthestDist) < ANGSTROM_FUDGE)
                            {
                                my @dist = ($dbSaddle[POSCAR_COORDINATES][$dbPairAtoms[$i]][0] - $dbSaddle[POSCAR_COORDINATES][$dbPairAtoms[$j]][0], 
                                            $dbSaddle[POSCAR_COORDINATES][$dbPairAtoms[$i]][1] - $dbSaddle[POSCAR_COORDINATES][$dbPairAtoms[$j]][1], 
                                            $dbSaddle[POSCAR_COORDINATES][$dbPairAtoms[$i]][2] - $dbSaddle[POSCAR_COORDINATES][$dbPairAtoms[$j]][2]);
                                my $dist = sqrt($dist[0] * $dist[0] + $dist[1] * $dist[1] + $dist[2] * $dist[2]);
                                if(abs($dist - $farthestSecondFarthestDist) < ANGSTROM_FUDGE)
                                {
                                    push(@dbAtomPairs, [$dbPairAtoms[$i], $dbPairAtoms[$j]]);
                                }
                            }
                        }
                    }   
                }
            }

            my @dbElements = CarElementTable(\@dbSaddle);
            my @minMatches = ();
            for(my $i = 0; $i < $saddle[POSCAR_TOTAL_ATOMS]; $i++)
            {
                for(my $j = 0; $j < $saddle[POSCAR_TOTAL_ATOMS]; $j++)
                {
                    $minMatches[$i][$j] = matchDescription($saddleElements[$i], $dbElements[$j]);
                }
            }
            
            if($debug)
            {
                print "    atom pairs: [$saddleFarthestAtom, $saddleSecondFarthestAtom] -> ";
                for(my $atomPair = 0; $atomPair < @dbAtomPairs; $atomPair++)
                {                 
                    print "($dbAtomPairs[$atomPair]->[0], $dbAtomPairs[$atomPair]->[1]) ";
                }
                print "\n";
            }


            for(my $atomPair = 0; $atomPair < @dbAtomPairs && $goodToGo; $atomPair++)
            {                 
            
                my @localDBSaddle = @{dclone(\@dbSaddle)};
            
                my $dbFarthestAtom = $dbAtomPairs[$atomPair]->[0];
                my $dbSecondFarthestAtom = $dbAtomPairs[$atomPair]->[1];
            
                if($debug)
                {
                    writeXYZ(\@localDBSaddle, "kdbaddpr-movie-$diri-$dbFarthestAtom-$dbSecondFarthestAtom.xyz");
                }

                # -- Rotation 1 --
                # Determine the unit vector from the saddle COM to the farthest DB atom.
                my @farthestDBUnit = ($localDBSaddle[POSCAR_COORDINATES][$dbFarthestAtom][0],
                                      $localDBSaddle[POSCAR_COORDINATES][$dbFarthestAtom][1],
                                      $localDBSaddle[POSCAR_COORDINATES][$dbFarthestAtom][2]);
                @farthestDBUnit = unitVector(\@farthestDBUnit);
                if($debug)
                {
                    print "    farthest DB Unit: ($farthestDBUnit[0], $farthestDBUnit[1], $farthestDBUnit[2])\n";
                    print "    farthest QC Unit: ($farthestSaddleUnit[0], $farthestSaddleUnit[1], $farthestSaddleUnit[2])\n";
                }
                # Calculate the rotation axis and angle needed to rotate farthestDBUnit onto farthestSaddleUnit
                my @crossProduct = crossProduct(\@farthestDBUnit, \@farthestSaddleUnit);
                # If the cross product is the zero vector, then the vectors are the same vector and we don't need to do the first rotation.
                my @axisFirstRotation = ();
                my $thetaFirstRotation = 0;
                if ($crossProduct[0] != 0 or $crossProduct[1] != 0 or $crossProduct[2] != 0)
                {
                    @axisFirstRotation = unitVector(\@crossProduct);
                    my $dp = dotProduct(\@farthestDBUnit, \@farthestSaddleUnit);
                    # WITHOUT THIS CHECK, acos(1) will be called through Math::Complex and slow everything down.
                    if ($dp > 0.99999999999)        
                    {
                        $dp = 0.99999999999;
                    }
                    if ($dp < -0.9999999)
                    {
                        $dp = -0.99999999;
                    }
                    $thetaFirstRotation = acos($dp);
                }
                else
                {
                    if($debug)
                    {
                        print "    fishy first rotation\n";
                    }
                    @axisFirstRotation = (1, 0, 0);
                    $thetaFirstRotation = 0;
                }
                if($debug)
                {
                    print "    first rotation: ($axisFirstRotation[0], $axisFirstRotation[1], $axisFirstRotation[2]) , -$thetaFirstRotation\n";
                }
                
                # Rotate $localDBSaddle so farthestDBUnit aligns with farthestSaddleUnit.
                for(my $a = 0; $a < $localDBSaddle[POSCAR_TOTAL_ATOMS]; $a++)
                {
                    my $x = $localDBSaddle[POSCAR_COORDINATES][$a][0];
                    my $y = $localDBSaddle[POSCAR_COORDINATES][$a][1];
                    my $z = $localDBSaddle[POSCAR_COORDINATES][$a][2];
                    my ($xp, $yp, $zp) = rotatePointAxis($x, $y, $z, $axisFirstRotation[0], $axisFirstRotation[1], $axisFirstRotation[2], -$thetaFirstRotation);
                    $localDBSaddle[POSCAR_COORDINATES][$a][0] = $xp;
                    $localDBSaddle[POSCAR_COORDINATES][$a][1] = $yp;
                    $localDBSaddle[POSCAR_COORDINATES][$a][2] = $zp;
                }

                if($debug)
                {
                    @farthestDBUnit = ($localDBSaddle[POSCAR_COORDINATES][$dbFarthestAtom][0],
                                       $localDBSaddle[POSCAR_COORDINATES][$dbFarthestAtom][1],
                                       $localDBSaddle[POSCAR_COORDINATES][$dbFarthestAtom][2]);
                    @farthestDBUnit = unitVector(\@farthestDBUnit);
                    my $dp = dotProduct(\@farthestDBUnit, \@farthestSaddleUnit);
                    print "    first rotation alignment: $dp\n";
                }

                if($debug)
                {
                    appendXYZ(\@localDBSaddle, "kdbaddpr-movie-$diri-$dbFarthestAtom-$dbSecondFarthestAtom.xyz");
                }
                
                
                # --- Rotation 2 ---
                my @CCb = ($saddle[POSCAR_COORDINATES][$saddleSecondFarthestAtom][0],
                           $saddle[POSCAR_COORDINATES][$saddleSecondFarthestAtom][1],
                           $saddle[POSCAR_COORDINATES][$saddleSecondFarthestAtom][2]);
                my $dp = dotProduct(\@CCb, \@farthestSaddleUnit);
                my @VCb = ($saddle[POSCAR_COORDINATES][$saddleSecondFarthestAtom][0] - ($dp * $farthestSaddleUnit[0]), 
                           $saddle[POSCAR_COORDINATES][$saddleSecondFarthestAtom][1] - ($dp * $farthestSaddleUnit[1]), 
                           $saddle[POSCAR_COORDINATES][$saddleSecondFarthestAtom][2] - ($dp * $farthestSaddleUnit[2]));
                my @CDb = ($localDBSaddle[POSCAR_COORDINATES][$dbSecondFarthestAtom][0],
                           $localDBSaddle[POSCAR_COORDINATES][$dbSecondFarthestAtom][1],
                           $localDBSaddle[POSCAR_COORDINATES][$dbSecondFarthestAtom][2]);
                $dp = dotProduct(\@CDb, \@farthestSaddleUnit);
                my @VDb = ($localDBSaddle[POSCAR_COORDINATES][$dbSecondFarthestAtom][0] - ($dp * $farthestSaddleUnit[0]), 
                           $localDBSaddle[POSCAR_COORDINATES][$dbSecondFarthestAtom][1] - ($dp * $farthestSaddleUnit[1]), 
                           $localDBSaddle[POSCAR_COORDINATES][$dbSecondFarthestAtom][2] - ($dp * $farthestSaddleUnit[2]));
                @VCb = unitVector(\@VCb);
                @VDb = unitVector(\@VDb);
                my @crossProduct = crossProduct(\@VDb, \@VCb);
                # If the cross product is the zero vector, then the vectors are the same vector and we don't need to do the second rotation.
                my @axisSecondRotation = ();
                my $thetaSecondRotation = 0;
                if ($crossProduct[0] != 0 or $crossProduct[1] != 0 or $crossProduct[2] != 0)
                {
                    @axisSecondRotation = unitVector(\@crossProduct);
                    $dp = dotProduct(\@VCb, \@VDb);
                    # WITHOUT THIS CHECK, acos(1) will be called through Math::Complex and slow everything down.
                    if ($dp > 0.99999999999)        
                    {
                        $dp = 0.99999999999;
                    }
                    if ($dp < -0.9999999)
                    {
                        $dp = -0.99999999;
                    }
                    $thetaSecondRotation = acos($dp);
                }
                else
                {
                    if($debug)
                    {
                        print "    fishy second rotation\n";
                    }
                    @axisSecondRotation = (1, 0, 0);
                    $thetaSecondRotation = 0;
                }
                
                if($debug)
                {
                    print "    Second dp: $dp\n";
                    print "    second rotation: ($axisSecondRotation[0], $axisSecondRotation[1], $axisSecondRotation[2]) , -$thetaSecondRotation\n";
                }
                
                # Rotate $localDBSaddle so farthestDBUnit aligns with farthestSaddleUnit.
                for(my $a = 0; $a < $localDBSaddle[POSCAR_TOTAL_ATOMS]; $a++)
                {
                    my $x = $localDBSaddle[POSCAR_COORDINATES][$a][0];
                    my $y = $localDBSaddle[POSCAR_COORDINATES][$a][1];
                    my $z = $localDBSaddle[POSCAR_COORDINATES][$a][2];
                    my ($xp, $yp, $zp) = rotatePointAxis($x, $y, $z, $axisSecondRotation[0], $axisSecondRotation[1], $axisSecondRotation[2], -$thetaSecondRotation);
                    $localDBSaddle[POSCAR_COORDINATES][$a][0] = $xp;
                    $localDBSaddle[POSCAR_COORDINATES][$a][1] = $yp;
                    $localDBSaddle[POSCAR_COORDINATES][$a][2] = $zp;
                }
                if($debug)
                {
                    appendXYZ(\@localDBSaddle, "kdbaddpr-movie-$diri-$dbFarthestAtom-$dbSecondFarthestAtom.xyz");
                }

                my @allDistances = ();
                my @costTable;
                for(my $i = 0; $i < $saddle[POSCAR_TOTAL_ATOMS]; $i++)
                {
                    for(my $j = 0; $j < $saddle[POSCAR_TOTAL_ATOMS]; $j++)
                    {
                        if($minMatches[$j][$i])
                        {
                            my @dist = ($localDBSaddle[POSCAR_COORDINATES][$i][0] - $saddle[POSCAR_COORDINATES][$j][0], 
                                        $localDBSaddle[POSCAR_COORDINATES][$i][1] - $saddle[POSCAR_COORDINATES][$j][1], 
                                        $localDBSaddle[POSCAR_COORDINATES][$i][2] - $saddle[POSCAR_COORDINATES][$j][2]);
                            my $dist = sqrt($dist[0] * $dist[0] + $dist[1] * $dist[1] + $dist[2] * $dist[2]);
                            $costTable[$i][$j] = $dist;
                            push @allDistances, [$i, $j, $dist];
                        }
                    }
                }
                
                # Get the best distance for each db atom.
                my %takenSaddle;
                my %takenDB;
                my @sorted = sort {$a->[2] <=> $b->[2]} @allDistances;
                my $matchCount = 0;
                my $index = 0;
                while($matchCount < $saddle[POSCAR_TOTAL_ATOMS])
                {
                    my $a = $sorted[$index][0];
                    my $b = $sorted[$index][1];
                    if(!exists $takenDB{$a} && !exists $takenSaddle{$b})
                    {
                        $takenDB{$a} = $b;
                        $takenSaddle{$b} = $a;
                        $matchCount++;
                    }
                    $index++;
                }
                
                my $magnitude = 0.777;
                my $dt = 0.1;

                # Minimize the torque.
                my $iterCount = 0;
                my @velocity = (0.0, 0.0, 0.0);
                while($magnitude > 0.001 && $iterCount < $maxIter)
                {
                    $iterCount++;
                    #print "$iterCount, $magnitude\n";
                    my @torque = (0, 0, 0);
                    for(my $a = 0; $a < $localDBSaddle[POSCAR_TOTAL_ATOMS]; $a++)
                    {
                        my($dx, $dy, $dz) = ($localDBSaddle[POSCAR_COORDINATES][$a][0] - $saddle[POSCAR_COORDINATES][$takenDB{$a}][0], 
                                             $localDBSaddle[POSCAR_COORDINATES][$a][1] - $saddle[POSCAR_COORDINATES][$takenDB{$a}][1], 
                                             $localDBSaddle[POSCAR_COORDINATES][$a][2] - $saddle[POSCAR_COORDINATES][$takenDB{$a}][2]);
                        my @F = ($dx * -$springK, $dy * -$springK, $dz * -$springK); 
                        my @T = crossProduct([$localDBSaddle[POSCAR_COORDINATES][$a][0],
                                              $localDBSaddle[POSCAR_COORDINATES][$a][1],
                                              $localDBSaddle[POSCAR_COORDINATES][$a][2]], \@F);
                        $torque[0] += $T[0];
                        $torque[1] += $T[1];
                        $torque[2] += $T[2];
                    }
                    my $dp = dotProduct(\@torque, \@velocity);
                    if($dp < 0.0)
                    {
                        $velocity[0] = $torque[0] * $dt;
                        $velocity[1] = $torque[1] * $dt;
                        $velocity[2] = $torque[2] * $dt;
                        $dt *= 0.9;
                    }
                    else
                    {
                        $velocity[0] += $torque[0] * $dt;
                        $velocity[1] += $torque[1] * $dt;
                        $velocity[2] += $torque[2] * $dt;
                    }
                    $magnitude = sqrt($velocity[0]*$velocity[0] + $velocity[1]*$velocity[1] + $velocity[2]*$velocity[2]);
                    my $theta = $dt * -$magnitude;
                    my $thetaMax = 0.1;
                    if($theta > $thetaMax)
                    {
                        $theta = $thetaMax;
                    }
                    if($theta < -$thetaMax)
                    {
                        $theta = -$thetaMax;
                    }
                    for(my $a = 0; $a < $localDBSaddle[POSCAR_TOTAL_ATOMS]; $a++)
                    {
                        my $x = $localDBSaddle[POSCAR_COORDINATES][$a][0];
                        my $y = $localDBSaddle[POSCAR_COORDINATES][$a][1];
                        my $z = $localDBSaddle[POSCAR_COORDINATES][$a][2];
                        my ($xp, $yp, $zp) = rotatePointAxis($x, $y, $z, $velocity[0], $velocity[1], $velocity[2], $theta);
                        $localDBSaddle[POSCAR_COORDINATES][$a][0] = $xp;
                        $localDBSaddle[POSCAR_COORDINATES][$a][1] = $yp;
                        $localDBSaddle[POSCAR_COORDINATES][$a][2] = $zp;
                    }
                    if($debug)
                    {
                        appendXYZ(\@localDBSaddle, "kdbaddpr-movie-$diri-$dbFarthestAtom-$dbSecondFarthestAtom.xyz");
                    }
                    
                }
                
                if($debug)
                {
                    writeXYZ(\@saddle, "kdbaddpr-movie-final-$diri-$dbFarthestAtom-$dbSecondFarthestAtom.xyz");
                    appendXYZ(\@localDBSaddle, "kdbaddpr-movie-final-$diri-$dbFarthestAtom-$dbSecondFarthestAtom.xyz");
                }

                # Torque is minimized, calculate a final distance for the saddle.
                my $saddleDist = 0.0;
                for(my $a = 0; $a < $localDBSaddle[POSCAR_TOTAL_ATOMS]; $a++)
                {
                    my $x = $localDBSaddle[POSCAR_COORDINATES][$a][0] - $saddle[POSCAR_COORDINATES][$takenDB{$a}][0];
                    my $y = $localDBSaddle[POSCAR_COORDINATES][$a][1] - $saddle[POSCAR_COORDINATES][$takenDB{$a}][1];
                    my $z = $localDBSaddle[POSCAR_COORDINATES][$a][2] - $saddle[POSCAR_COORDINATES][$takenDB{$a}][2];
                    my $tempDist = sqrt($x*$x + $y*$y + $z*$z);
                    if($tempDist > $saddleDist)
                    {
                        $saddleDist = $tempDist;
                    }
                }

                if($debug)
                {
                    print "    score ($dbAtomPairs[$atomPair]->[0], $dbAtomPairs[$atomPair]->[1]): $saddleDist\n";
                }
                
                if($saddleDist < $duplicateCutoff)
                {
                    print "duplicate of $dbase/$desc/$diri. ($saddleDist)\n";
                    $goodToGo = 0;
                }
                
                if($saddleDist < $worstDist)
                {
                    $worstDist = $saddleDist;
                }
                
            } # End loop over atom pairs.
        }
        $diri++;
    }
    # If everything is okay, add the process to the database.
    if($goodToGo == 1)
    {
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

        # Touch the mirror file if the mirror flag is set.
        if($mirror_flag == 1)
        {
            system("touch $dbase/$desc/$diri/mirror");
        }
        
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
        printf("good. (%10f)\n", $worstDist);
    }

















