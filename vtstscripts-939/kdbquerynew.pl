#!/usr/bin/env perl


#===================================================================================================
# Directives
#---------------------------------------------------------------------------------------------------
    use strict;
    use warnings;
    # Add the directory containing the script to the lib path.
    use FindBin qw($Bin);
    use lib "$Bin";


#===================================================================================================
# Includes
#---------------------------------------------------------------------------------------------------
    use IO::Handle;
    use File::Path;
    use File::Copy;
    use Math::Trig;
    use Vasp;
    use kdbutil;
    use Getopt::Long;
    use Storable qw(dclone);
#use Carp 'verbose';
#$SIG{'INT'} = sub { Carp::confess( @_ ) };


#===================================================================================================
# Global Constants
#---------------------------------------------------------------------------------------------------

    my $debug = '';
    my $dedupe = '';
    my $allmob = '';
    my $ScoreCutoff = 0.1;
    my $NeighborFudgePercent = 0.2;
    my $AngstromFudge = 0.1;
    GetOptions("debug"=>  \$debug, 
               "dedupe"=> \$dedupe, 
               "allmob"=> \$allmob, 
               "score=f", \$ScoreCutoff, 
               "af=f",    \$AngstromFudge, 
               "nfp=f",   \$NeighborFudgePercent);
        
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

#===================================================================================================
# Global Variables
#---------------------------------------------------------------------------------------------------
    
    # The maximum number of iterations to minimize the torque.
    my $maxIter = 64;

    # The spring constant used for torque minimization.
    my $springK = 1.0;    
    my $FineCriteria = 0.0001;
    
#===================================================================================================
# Initialization
#---------------------------------------------------------------------------------------------------

    STDOUT->autoflush(1);

    # Get the home directory of the database.
    my $KDBdir = kdbHome();
    
    # $CarFileName is the name of the POSCAR file we are attempting to match against the DB. This
    # is the first argument to kdbquery.pl.
    my $CarFileName = $ARGV[0]; 
    
    # If no POSCAR filename is given, die gracefully with usage instructions.
    if (!$CarFileName)
    {
        die "\nkdbquery.pl usage:\n\n" . 
        "     kdbquery.pl <CAR file of VASP configuration to match>\n\n";
    }
    
    # Attempt to load the POSCAR.
    my $Car = [read_poscar($CarFileName)];
    
    # Get the number of atoms inside the car file.
    my $NumCarAtoms = $Car->[POSCAR_TOTAL_ATOMS];

    # --- Create a local directory (kdbmatches) to spit out the appropriately transformed  --- 
    # --- poscars there.                                                                   ---
    # Remove the directory if it exists, and create it again.
    rmtree('kdbmatches');
    mkpath('kdbmatches');

    # Create the distance table for the target POSCAR.
    print "Calculating Car atomic distance table... ";
    my @CarDistanceTable = atomicDistanceTableDir($Car);
    print "done.\n";
    
    # Create a list of atoms for each car atom sorted by distance from that car atom.
    my %sortedCarAtoms = ();
    my @unsorted = (0..$NumCarAtoms - 1);
    for (my $i = 0; $i < $NumCarAtoms; $i++)
    {
        my @sorted = sort { $CarDistanceTable[$i][$a] <=> $CarDistanceTable[$i][$b] } @unsorted;
        $sortedCarAtoms{$i} = \@sorted;
    }


    # Store the list of elements inside $Desc.
    my $Desc = $Car->[POSCAR_DESCRIPTION];

    my @CarElements = CarElementTable($Car);

    # Create a copy of $Car in cartesian format.
    my @CarKar = DirKar($Car);
    
    # Get the type count for Car.
    my $carTypeCount = neighborTypeCount($Car, \@CarDistanceTable, $NeighborFudgePercent);

#===================================================================================================
# Main algorithm.
#---------------------------------------------------------------------------------------------------

    # Get the list of $Desc matching process subdirectories.
    my @tempDirs = getProcessDirs($KDBdir, $Desc);
    
    # In order to do this one minimum at a time, make an array of the each minimum and saddle.
    # Each process in the db, then, will have two entries.
    my @AllDirs = ();
    foreach my $dir (@tempDirs)
    {
        push @AllDirs, ["$dir", "$dir/min1.xyz", "$dir/saddle.xyz", "$dir/e1"];
        push @AllDirs, ["$dir", "$dir/min2.xyz", "$dir/saddle.xyz", "$dir/e2"];
    }  
    
    # Keep count of the number of matches.
    my $numMatches = 0;

    # Keep track of the matched saddles for duplicate comparisons.
    my @saddleMatches = ();

    # Loop over each directory.
    my $dirIndex = 0;
    foreach my $dirPointer (@AllDirs)
    {

        # Dereference $dirPointer into @dir.
        my @dir = @{$dirPointer};
        
        # Notify the user we are checking the minimum.
        print "\nChecking $dir[1]";

        # Load the data for each minimum. $dir[1] is the path to the db minimum.
        my @dbMinKar = loadXYZ2CAR($dir[1]);
        my @dbMinCar = @{dclone(\@dbMinKar)};
        $dbMinCar[POSCAR_BASIS] = $Car->[POSCAR_BASIS];
        $dbMinCar[POSCAR_LATTICE] = $Car->[POSCAR_LATTICE];
        kardir($dbMinCar[0], $dbMinCar[POSCAR_BASIS], $dbMinCar[POSCAR_LATTICE], $dbMinCar[POSCAR_TOTAL_ATOMS]);
        my $dbNumAtoms = $dbMinCar[POSCAR_TOTAL_ATOMS];
        open(MOBILE, "<$dir[0]/mobile");
        my @mobile = <MOBILE>;
        chomp @mobile;
        close MOBILE;

        my @minElements = CarElementTable(\@dbMinKar);
        my @minMatches = ();
        
        # Create the match table.  This tells us if the atom j of the car we are querying with
        # is a match with atom i of the database minimum, by element. So, if minMatches[$i][$j]
        # is true, then atom i of the database minimum is an elemental match with the jth atom of
        # the Car we are querying with.
        for(my $i = 0; $i < $dbMinKar[POSCAR_TOTAL_ATOMS]; $i++)
        {
            for(my $j = 0; $j < $NumCarAtoms; $j++)
            {
                if($minElements[$i] eq $CarElements[$j])
                {
                    $minMatches[$i][$j] = 1;
                }
                else
                {
                    $minMatches[$i][$j] = 0;
                }
                # Use the following for any-metal matching.
                #$minMatches[$i][$j] = matchDescription($minElements[$i], $CarElements[$j]);
            }
        }

        # Create a distance table for the database minimum.
        my @dbMinDistanceTable = atomicDistanceTableDir(\@dbMinCar);

        # Get the type count for this kdb process minimum.
        my $dbTypeCount = neighborTypeCount(\@dbMinCar, \@dbMinDistanceTable, $NeighborFudgePercent);
        
#        if (!$allmob) # Use the mobile atom with the least number of neighbors.
#        {
#            my $least = 1000;
#            my $one = 0;
#            foreach my $mob(@mobile)
#            {
#                if($dbTypeCount->{$mob}->{'total'} < $least)
#                {
#                    $least = $dbTypeCount->{$mob}->{'total'};
#                    $one = $mob;                
#                }
#            }
#            @mobile = ($one);
#        }
        
        # Determine the process radius for each mobile atom in this kdb process minimum.
        my %mobileRadii = ();
        foreach my $mob(@mobile)
        {
            my $maxDist = 0;
            for (my $a = 0; $a < $dbNumAtoms; $a++)
            {
                if($dbMinDistanceTable[$mob][$a] > $maxDist)
                {
                    $maxDist = $dbMinDistanceTable[$mob][$a];
                    $mobileRadii{$mob} = $maxDist;
                }
            }                
        }
        
        my @allMappings = ();
        my $mobWithLeastMappings = -1;
        my $leastMappings = 100000000;
        foreach my $mob (@mobile)
        {
            if($dbTypeCount->{$mob}->{'total'} < 1)
            {
                next;
            }
            my @mappings = ();
            for my $qa(keys %$carTypeCount)
            {
                if (neighborTypeCountDifference($carTypeCount->{$qa}, $dbTypeCount->{$mob}) < 1)
                {
                    push @mappings, [{$mob => $qa}, {$qa => $mob}];
                }
            }
            my $numMaps = @mappings;
            if($numMaps < $leastMappings)
            {
                $mobWithLeastMappings = $mob;
                $leastMappings = $numMaps;
                @allMappings = @mappings;
            }
            
        }

        my $mob = $mobWithLeastMappings;
        my @newMappings = ();
        for(my $dba = 0; $dba < $dbNumAtoms; $dba++)
        {

            MAPPING: for (my $i = 0; $i < @allMappings; $i++)
            {


                my %maps = %{$allMappings[$i][0]};
                my %rmaps = %{$allMappings[$i][1]};


                if (exists $maps{$dba})
                {
                    my %map = %maps;
                    my %rmap = %rmaps;
                    push @newMappings, [\%map, \%rmap];
                    next MAPPING;
                }

                my $scai = 0; # sortedCarAtom index
                my $qa  = $sortedCarAtoms{$maps{$mob}}->[$scai];
                my $length = @{$sortedCarAtoms{$maps{$mob}}};
                QA: while($scai < $length && $CarDistanceTable[$maps{$mob}][$qa] < $mobileRadii{$mob} + $AngstromFudge)
                {
                    if (exists $rmaps{$qa})
                    {
                        $scai++;
                        $qa = $sortedCarAtoms{$maps{$mob}}->[$scai];
                        next QA;
                    }
                    foreach my $dbm(keys %maps)
                    {
                        if(abs($dbMinDistanceTable[$dbm][$dba] - $CarDistanceTable[$maps{$dbm}][$qa]) > $AngstromFudge)
                        {
                            $scai++;
                            $qa = $sortedCarAtoms{$maps{$mob}}->[$scai];
                            next QA;
                        }
                        if(!$minMatches[$dba][$qa])
                        {
                            $scai++;
                            $qa = $sortedCarAtoms{$maps{$mob}}->[$scai];
                            next QA;
                        }
                    }
                    my %map = %maps;
                    my %rmap = %rmaps;
                    $map{$dba} = $qa;
                    $rmap{$qa} = $dba;
                    push @newMappings, [\%map, \%rmap];
                    $scai++;
                    $qa = $sortedCarAtoms{$maps{$mob}}->[$scai];
                }
            }
            @allMappings = @newMappings;
            @newMappings = ();
        }
#        @allMappings = (@allMappings, @mappings);
    #}


#===============================================================================================
# Main loop over mappings.
#-----------------------------------------------------------------------------------------------
        foreach my $mapping(@allMappings)
        {

            print ".";
            foreach my $mirrorFlag((0, 1))
            {
                my %maps = %{$mapping->[0]};
                my $db_a = $mobile[0];

                # Find the two atoms closest to db_a.
                my $db_b = 0;
                my $db_c = 0;
                # Find the closest:
                my $min = 99999999999;
                for(my $i = 0; $i < $dbNumAtoms; $i++)
                {
                    if($i != $db_a)
                    {
                        my $d = $dbMinDistanceTable[$db_a][$i];
                        if($d < $min)
                        {
                            $db_b = $i;
                            $min = $d;
                        }
                    }
                }
                # Find the second closest:
                my @db_ab = intraAtomicVectorDirKar(\@dbMinCar, $db_a, $db_b); 
                my @db_ab_u = unitVector(\@db_ab);
                $min = 99999999999;
                for(my $i = 0; $i < $dbNumAtoms; $i++)
                {
                    if($i != $db_a && $i != $db_b)
                    {
                        my $d = $dbMinDistanceTable[$db_a][$i];
                        if($d < $min)
                        {
                            my $c = $i;
                            my @db_ac = intraAtomicVectorDirKar(\@dbMinCar, $db_a, $c);        
                            my @db_ac_u = unitVector(\@db_ac);
                            my $dp = abs(dotProduct(\@db_ab_u, \@db_ac_u));
                            # Make sure the three atoms are not in a line.
                            if($dp < 0.9)
                            {
                                $db_c = $i;
                                $min = $d;
                            }
                        }
                    }
                }
                my $q_a = $maps{$db_a};
                my $q_b = $maps{$db_b};
                my $q_c = $maps{$db_c};
	            if($debug)
	            {
	                print "\n            db_a, b, c: $db_a, $db_b, $db_c\n";
	                print "            q_a, b, c:  $q_a, $q_b, $q_c\n";
	            }                        

#===============================================================================================
# Perform initial rough rotation:
#-----------------------------------------------------------------------------------------------

                my @localdbMinKar = @{dclone(\@dbMinKar)};

                # Create mirror if mirrorFlag.
                if($mirrorFlag)
                {
                    for(my $i = 1; $i < $dbNumAtoms; $i++)
                    {
                        $localdbMinKar[POSCAR_COORDINATES][$i][0] += 2.0 * ($localdbMinKar[POSCAR_COORDINATES][0][0] - $localdbMinKar[POSCAR_COORDINATES][$i][0]);
                        $localdbMinKar[POSCAR_COORDINATES][$i][1] += 2.0 * ($localdbMinKar[POSCAR_COORDINATES][0][1] - $localdbMinKar[POSCAR_COORDINATES][$i][1]);
                        $localdbMinKar[POSCAR_COORDINATES][$i][2] += 2.0 * ($localdbMinKar[POSCAR_COORDINATES][0][2] - $localdbMinKar[POSCAR_COORDINATES][$i][2]);
                    }
                }
                
                my @origDbMinKar = @{dclone(\@localdbMinKar)};

                my @localdbMinCar = @{dclone(\@localdbMinKar)};
                $localdbMinCar[POSCAR_BASIS] = $Car->[POSCAR_BASIS];
                $localdbMinCar[POSCAR_LATTICE] = $Car->[POSCAR_LATTICE];
                kardir($localdbMinCar[0], $localdbMinCar[POSCAR_BASIS], $localdbMinCar[POSCAR_LATTICE], $localdbMinCar[POSCAR_TOTAL_ATOMS]);
                
                my @translation = (0, 0, 0);
                
                # Find a translation for the minimum.
	            @translation = interAtomicVectorDirKar(\@localdbMinCar, $Car, $db_a, $q_a, $Car);
	            
	            my $xd = $translation[0] - ($CarKar[POSCAR_COORDINATES][$q_a][0] - $localdbMinKar[POSCAR_COORDINATES][$db_a][0]);
	            my $yd = $translation[1] - ($CarKar[POSCAR_COORDINATES][$q_a][1] - $localdbMinKar[POSCAR_COORDINATES][$db_a][1]);
	            my $zd = $translation[2] - ($CarKar[POSCAR_COORDINATES][$q_a][2] - $localdbMinKar[POSCAR_COORDINATES][$db_a][2]);
    #	        print "\n\nD: $xd $yd $zd\n\n";
	            
	            
	            if($debug)
	            {
	                print "            translation: @translation\n";
	            }                        


                # Shift the new coords by @translation.
                for(my $a = 0; $a < $dbNumAtoms; $a++)
                {
                    $localdbMinKar[POSCAR_COORDINATES][$a][0] += $translation[0];
                    $localdbMinKar[POSCAR_COORDINATES][$a][1] += $translation[1];
                    $localdbMinKar[POSCAR_COORDINATES][$a][2] += $translation[2];
                }
                
                # If $debugging,start the .xyz file.
                if($debug)
                {
                    writeXYZ(\@localdbMinKar, "$dirIndex-$q_b-$q_c.xyz");
                }
                
                # Determine the unit vector from db_a to db_b.
                @db_ab_u = intraAtomicVectorDirKar(\@localdbMinCar, $db_a, $db_b);
                @db_ab_u = unitVector(\@db_ab_u);
                # Determine the unit vector from q_a to q_b.
                my @q_ab_u = intraAtomicVectorDirKar($Car, $q_a, $q_b);
                @q_ab_u = unitVector(\@q_ab_u);

                # Calculate the rotation axis and angle needed to rotate db_ab_u onto q_ab_u
                my @crossProduct = crossProduct(\@db_ab_u, \@q_ab_u);
                
                # If the cross product is the zero vector, then the vectors are the same vector and we don't need to do the first rotation.
                my @axisFirstRotation = ();
                my $thetaFirstRotation = 0;
                if ($crossProduct[0] != 0 or $crossProduct[1] != 0 or $crossProduct[2] != 0)
                {
                    @axisFirstRotation = unitVector(\@crossProduct);
                    my $dp = dotProduct(\@db_ab_u, \@q_ab_u);
                    # WITHOUT THIS CHECK, acos(1) will be called through Math::Complex and slow everything down.
                    if ($dp > 0.99999999999)        
                    {
                        $dp = 0.99999999999;
                    }
                    $thetaFirstRotation = acos($dp);
                }
                else
                {
                    @axisFirstRotation = (1, 0, 0);
                    $thetaFirstRotation = 0;
                }
                
                if($debug)
                {
                    print "            thetaFirstRotation: $thetaFirstRotation\n";
                }

                # Rotate localdbMinKar so db_ab_u aligns with q_ab_u.
                for(my $a = 0; $a < $dbNumAtoms; $a++)
                {
                    if($a != $db_a)
                    {
                        my $x = $localdbMinKar[POSCAR_COORDINATES][$a][0] - $localdbMinKar[POSCAR_COORDINATES][$db_a][0];
                        my $y = $localdbMinKar[POSCAR_COORDINATES][$a][1] - $localdbMinKar[POSCAR_COORDINATES][$db_a][1];
                        my $z = $localdbMinKar[POSCAR_COORDINATES][$a][2] - $localdbMinKar[POSCAR_COORDINATES][$db_a][2];
                        my ($xp, $yp, $zp) = rotatePointAxis($x, $y, $z, $axisFirstRotation[0], $axisFirstRotation[1], $axisFirstRotation[2], -$thetaFirstRotation);
                        $localdbMinKar[POSCAR_COORDINATES][$a][0] = $xp + $localdbMinKar[POSCAR_COORDINATES][$db_a][0];
                        $localdbMinKar[POSCAR_COORDINATES][$a][1] = $yp + $localdbMinKar[POSCAR_COORDINATES][$db_a][1];
                        $localdbMinKar[POSCAR_COORDINATES][$a][2] = $zp + $localdbMinKar[POSCAR_COORDINATES][$db_a][2];
                    }
                }
                
                if($debug)
                {
                    appendXYZ(\@localdbMinKar, "$dirIndex-$q_b-$q_c.xyz");
                }
                # Update dbMinCar.
                @localdbMinCar = @{dclone(\@localdbMinKar)};
                $localdbMinCar[POSCAR_BASIS] = $Car->[POSCAR_BASIS];
                $localdbMinCar[POSCAR_LATTICE] = $Car->[POSCAR_LATTICE];
                kardir($localdbMinCar[0], $localdbMinCar[POSCAR_BASIS], $localdbMinCar[POSCAR_LATTICE], $localdbMinCar[POSCAR_TOTAL_ATOMS]);
                # Calculate the angle for matching q_c and db_c.
                my @q_ab = intraAtomicVectorDirKar($Car, $q_a, $q_b);        
                @q_ab_u = unitVector(\@q_ab);
                my @db_ac = intraAtomicVectorDirKar(\@localdbMinCar, $db_a, $db_c);        
                my $vm = dotProduct(\@db_ac, \@q_ab_u);
                my @db_c_ab = ($localdbMinKar[POSCAR_COORDINATES][$db_c][0] - ($localdbMinKar[POSCAR_COORDINATES][$db_a][0] + $q_ab_u[0] * $vm),
                               $localdbMinKar[POSCAR_COORDINATES][$db_c][1] - ($localdbMinKar[POSCAR_COORDINATES][$db_a][1] + $q_ab_u[1] * $vm),
                               $localdbMinKar[POSCAR_COORDINATES][$db_c][2] - ($localdbMinKar[POSCAR_COORDINATES][$db_a][2] + $q_ab_u[2] * $vm));   
                kardir([\@db_c_ab], $Car->[POSCAR_BASIS], $Car->[POSCAR_LATTICE], 1);
                set_bc([\@db_c_ab], 1);
                dirkar([\@db_c_ab], $Car->[POSCAR_BASIS], $Car->[POSCAR_LATTICE], 1);                
                my @db_c_ab_u = unitVector(\@db_c_ab);
                my @q_c_ab = intraAtomicVectorDirKar($Car, $q_a, $q_c);    
                $vm = dotProduct(\@q_c_ab, \@q_ab_u);
                @q_c_ab = ($CarKar[POSCAR_COORDINATES][$q_c][0] - ($CarKar[POSCAR_COORDINATES][$q_a][0] + $q_ab_u[0] * $vm),
                           $CarKar[POSCAR_COORDINATES][$q_c][1] - ($CarKar[POSCAR_COORDINATES][$q_a][1] + $q_ab_u[1] * $vm),
                           $CarKar[POSCAR_COORDINATES][$q_c][2] - ($CarKar[POSCAR_COORDINATES][$q_a][2] + $q_ab_u[2] * $vm));     
                kardir([\@q_c_ab], $Car->[POSCAR_BASIS], $Car->[POSCAR_LATTICE], 1);
                set_bc([\@q_c_ab], 1);
                dirkar([\@q_c_ab], $Car->[POSCAR_BASIS], $Car->[POSCAR_LATTICE], 1);                
                my @q_c_ab_u = unitVector(\@q_c_ab);
                my $thetaSecondRotation = 0;
                my @axisSecondRotation = (1, 0, 0);
                if($q_c_ab[0] * $q_c_ab[0] + $q_c_ab[1] * $q_c_ab[1] + $q_c_ab[2] * $q_c_ab[2] != 0) 
                {
                    $thetaSecondRotation = acos(dotProduct(\@q_c_ab_u, \@db_c_ab_u));
                    if($thetaSecondRotation != 0)
                    {                    
                        if(abs(pi - $thetaSecondRotation) < 0.00001)
                        {
                            @axisSecondRotation = @q_ab_u;
                        }
                        else
                        {
                            @axisSecondRotation = crossProduct(\@db_c_ab_u, \@q_c_ab_u);
                            @axisSecondRotation = unitVector(\@axisSecondRotation);
                        }
                    }
                }
                else
                {
                    next;
                }
                if($debug)
                {
                    print "            thetaSecondRotation: $thetaSecondRotation\n";
                    print "            axisSecondRotation: @axisSecondRotation\n";
                    print "            q_c_ab_u: @q_c_ab_u\n";
                    print "            db_c_ab: @db_c_ab\n";
                    print "            q_ab_u: @q_ab_u\n";
                    my $dpp = dotProduct(\@q_ab_u, \@axisSecondRotation);
                    print "            dpp: $dpp\n";
                }
                        
                # Rotate localdbMinKar so the third best matches align.
                if($thetaSecondRotation != 0)
                {                    
                    for(my $a = 0; $a < $dbNumAtoms; $a++)
                    {
                        if($a != $db_a && $a != $db_b)
                        {
                            my $x = $localdbMinKar[POSCAR_COORDINATES][$a][0] - $localdbMinKar[POSCAR_COORDINATES][$db_a][0];
                            my $y = $localdbMinKar[POSCAR_COORDINATES][$a][1] - $localdbMinKar[POSCAR_COORDINATES][$db_a][1];
                            my $z = $localdbMinKar[POSCAR_COORDINATES][$a][2] - $localdbMinKar[POSCAR_COORDINATES][$db_a][2];
                            my ($xp, $yp, $zp) = rotatePointAxis($x, $y, $z, $axisSecondRotation[0], $axisSecondRotation[1], $axisSecondRotation[2], -$thetaSecondRotation);
                            $localdbMinKar[POSCAR_COORDINATES][$a][0] = $xp + $localdbMinKar[POSCAR_COORDINATES][$db_a][0];
                            $localdbMinKar[POSCAR_COORDINATES][$a][1] = $yp + $localdbMinKar[POSCAR_COORDINATES][$db_a][1];
                            $localdbMinKar[POSCAR_COORDINATES][$a][2] = $zp + $localdbMinKar[POSCAR_COORDINATES][$db_a][2];
                        }
                    }
                }
                if($debug && 1)
                {
                    appendXYZ(\@localdbMinKar, "$dirIndex-$q_b-$q_c.xyz");

                }
                

    #===============================================================================================
    # Perform torque minimization:
    #-----------------------------------------------------------------------------------------------
                
                # Keep track of each torque vector and theta size.
                my @torques = ();
                my @thetas = ();
                my @steps = ();
                my $torqueMagnitude = 0.777;
                my $forceMagnitude = 0.777;
                my @torqueVelocity = (0.0, 0.0, 0.0);
                my @forceVelocity = (0.0, 0.0, 0.0);
                my $tdt = 0.2;
                my $fdt = 0.2;
                # Minimize the torque and translational force.
                my $iterCount = 0;
#                print "\n\n";
                while($torqueMagnitude + $forceMagnitude > $FineCriteria && $iterCount < $maxIter)
                {
                    $iterCount++;
#                    print "$iterCount: $torqueMagnitude + $forceMagnitude > $FineCriteria\n";
    #                print "$numMatches, $iterCount, $torqueMagnitude, $forceMagnitude\n";
                    my @torque = (0, 0, 0);
                    my @force = (0, 0, 0);
                    @localdbMinCar = @{dclone(\@localdbMinKar)};
                    $localdbMinCar[POSCAR_BASIS] = $Car->[POSCAR_BASIS];
                    $localdbMinCar[POSCAR_LATTICE] = $Car->[POSCAR_LATTICE];
                    kardir($localdbMinCar[0], $localdbMinCar[POSCAR_BASIS], $localdbMinCar[POSCAR_LATTICE], $localdbMinCar[POSCAR_TOTAL_ATOMS]);

                    
                    for(my $a = 0; $a < $dbNumAtoms; $a++)
                    {
                        my $aCoords = [[$Car->[0][$maps{$a}][0], $Car->[0][$maps{$a}][1], $Car->[0][$maps{$a}][2]]];
                        my $bCoords = [[$localdbMinCar[0][$a][0], $localdbMinCar[0][$a][1], $localdbMinCar[0][$a][2]]];
                        my @diff = pbc_difference($bCoords, $aCoords, 1);
                        dirkar(@diff, $Car->[POSCAR_BASIS], $Car->[POSCAR_LATTICE], 1);
		                my($dx, $dy, $dz) = @{$diff[0][0]};
                        my @F = ($dx * -$springK, $dy * -$springK, $dz * -$springK); 
                        my @T = crossProduct([$localdbMinKar[POSCAR_COORDINATES][$a][0] - $localdbMinKar[POSCAR_COORDINATES][$db_a][0], 
                                              $localdbMinKar[POSCAR_COORDINATES][$a][1] - $localdbMinKar[POSCAR_COORDINATES][$db_a][1],
                                              $localdbMinKar[POSCAR_COORDINATES][$a][2] - $localdbMinKar[POSCAR_COORDINATES][$db_a][2]], \@F);
                        $torque[0] += $T[0];
                        $torque[1] += $T[1];
                        $torque[2] += $T[2];
                        $force[0] += $F[0];
                        $force[1] += $F[1];
                        $force[2] += $F[2];
                    }
                    $torqueMagnitude = sqrt($torque[0]*$torque[0] + $torque[1]*$torque[1] + $torque[2]*$torque[2]);
                    $forceMagnitude = sqrt($force[0]*$force[0] + $force[1]*$force[1] + $force[2]*$force[2]);
                    
                    if($torque[0] * $torqueVelocity[0] + $torque[1] * $torqueVelocity[1] + $torque[2] * $torqueVelocity[2] < 0.0)
                    {
                        @torqueVelocity = (0.0, 0.0, 0.0);
                        $tdt *= 0.9;
                    }
                    $torqueVelocity[0] += $torque[0] * $tdt;
                    $torqueVelocity[1] += $torque[1] * $tdt;
                    $torqueVelocity[2] += $torque[2] * $tdt;

                    if($force[0] * $forceVelocity[0] + $force[1] * $forceVelocity[1] + $force[2] * $forceVelocity[2] < 0.0)
                    {
                        @forceVelocity = (0.0, 0.0, 0.0);
                        $fdt *= 0.9;
                    }
                    $forceVelocity[0] += $force[0] * $fdt;
                    $forceVelocity[1] += $force[1] * $fdt;
                    $forceVelocity[2] += $force[2] * $fdt;
                    
                    if($debug && 1)
                    {
                        print "                Torque magnitude: $torqueMagnitude\n";
                    }
                    my $theta = $tdt * -sqrt($torqueVelocity[0]*$torqueVelocity[0] + $torqueVelocity[1]*$torqueVelocity[1] + $torqueVelocity[2]*$torqueVelocity[2]);
                    my $thetaMax = 0.01; #XXX Parameterize (and increase default?) this
                    if($theta > $thetaMax)
                    {
                        $theta = $thetaMax;
                    }
                    if($theta < -$thetaMax)
                    {
                        $theta = -$thetaMax;
                    }
                    my @step = ($forceVelocity[0] * $fdt, $forceVelocity[1] * $fdt, $forceVelocity[2] * $fdt);
                    my $mag = sqrt($step[0]*$step[0] + $step[1]*$step[1] + $step[2]*$step[2]);
                    my $stepMax = 0.1;
                    if($mag > $stepMax)
                    {
                        $step[0] = $stepMax * $step[0] / $mag;
                        $step[1] = $stepMax * $step[1] / $mag;
                        $step[2] = $stepMax * $step[2] / $mag;
                    }
                    
                    my $tx = $localdbMinKar[POSCAR_COORDINATES][$db_a][0];
                    my $ty = $localdbMinKar[POSCAR_COORDINATES][$db_a][1];
                    my $tz = $localdbMinKar[POSCAR_COORDINATES][$db_a][2];
                
                    for(my $a = 0; $a < $dbNumAtoms; $a++)
                    {
                        $localdbMinKar[POSCAR_COORDINATES][$a][0] -= $tx;
                        $localdbMinKar[POSCAR_COORDINATES][$a][1] -= $ty;
                        $localdbMinKar[POSCAR_COORDINATES][$a][2] -= $tz;
                    }

                    $localdbMinKar[POSCAR_COORDINATES] = rotatePointsAxis($localdbMinKar[POSCAR_COORDINATES], $dbNumAtoms, $torque[0], $torque[1], $torque[2], $theta);

                    for(my $a = 0; $a < $dbNumAtoms; $a++)
                    {
                        $localdbMinKar[POSCAR_COORDINATES][$a][0] += $tx;
                        $localdbMinKar[POSCAR_COORDINATES][$a][1] += $ty;
                        $localdbMinKar[POSCAR_COORDINATES][$a][2] += $tz;
                        $localdbMinKar[POSCAR_COORDINATES][$a][0] += $step[0];
                        $localdbMinKar[POSCAR_COORDINATES][$a][1] += $step[1];
                        $localdbMinKar[POSCAR_COORDINATES][$a][2] += $step[2];
                    }

                    push(@torques, \@torque);
                    push(@thetas, $theta);
                    push(@steps, \@step);
                    if($debug && 1)
                    {
                        appendXYZ(\@localdbMinKar, "$dirIndex-$q_b-$q_c.xyz");

                    }
                }
                
    #===============================================================================================
    # Calculate a score for this alignment.
    #-----------------------------------------------------------------------------------------------
                # Create a table of distances from each db atom to atom in @CarKar.
                my @allDistances = ();
                my @costTable = ();
                my @minCar = @{dclone(\@localdbMinKar)};
                $minCar[POSCAR_BASIS] = $Car->[POSCAR_BASIS];
                $minCar[POSCAR_LATTICE] = $Car->[POSCAR_LATTICE];
                kardir($minCar[0], $minCar[POSCAR_BASIS], $minCar[POSCAR_LATTICE], $minCar[POSCAR_TOTAL_ATOMS]);
                foreach my $a(keys %maps)
                {
                    my $dist = interAtomicDistanceDirKar(\@minCar, $Car, $a, $maps{$a}, $Car);
                    $costTable[$a][$maps{$a}] = $dist;
                }

                # Calculate the total score for this angle.
                my $score = 0;
                for(my $i = 0; $i < $dbNumAtoms; $i++)
                {
                    $score = max($score, $costTable[$i][$maps{$i}]);
                }
                
    #===============================================================================================
    # Calculate store this suggestion if the score is good enough.
    #-----------------------------------------------------------------------------------------------
                if($score < $ScoreCutoff)
                {
    #                foreach my $a(keys %maps)            
    #                {
    #                    foreach my $b(keys %maps)            
    #                    {
    #                            my $d = abs($dbMinDistanceTable[$a][$b] - $CarDistanceTable[$maps{$a}][$maps{$b}]);
    #                            printf("\n%10f", $d);
    #                    }
    #                }
    #                print "\n$numMatches: $score, $iterCount\n";#XXX
    #                exit;#XXX
                    # Load the db saddle and mobile atoms.
                    my @dbSaddleKar = loadXYZ2CAR($dir[0] . "/saddle.xyz");
                    
                    if($mirrorFlag)
                    {
                        for(my $i = 1; $i < $dbNumAtoms; $i++)
                        {
                        
                            $dbSaddleKar[POSCAR_COORDINATES][$i][0] += 2.0 * ($dbSaddleKar[POSCAR_COORDINATES][0][0] - $dbSaddleKar[POSCAR_COORDINATES][$i][0]);
                            $dbSaddleKar[POSCAR_COORDINATES][$i][1] += 2.0 * ($dbSaddleKar[POSCAR_COORDINATES][0][1] - $dbSaddleKar[POSCAR_COORDINATES][$i][1]);
                            $dbSaddleKar[POSCAR_COORDINATES][$i][2] += 2.0 * ($dbSaddleKar[POSCAR_COORDINATES][0][2] - $dbSaddleKar[POSCAR_COORDINATES][$i][2]);
                        }
                    }

                    # Translate the mobile atoms.
                    for(my $j = 0; $j < $dbNumAtoms; $j++)
                    {
                        $dbSaddleKar[POSCAR_COORDINATES][$j][0] += $translation[0];
                        $dbSaddleKar[POSCAR_COORDINATES][$j][1] += $translation[1];
                        $dbSaddleKar[POSCAR_COORDINATES][$j][2] += $translation[2];
                    }
                
                    my $rotPointX = $origDbMinKar[POSCAR_COORDINATES][$db_a][0] + $translation[0];
                    my $rotPointY = $origDbMinKar[POSCAR_COORDINATES][$db_a][1] + $translation[1];
                    my $rotPointZ = $origDbMinKar[POSCAR_COORDINATES][$db_a][2] + $translation[2];

                    # Rotate dbSaddleKar according to the first rotation.
                    for(my $a = 0; $a < $dbNumAtoms; $a++)
                    {
                        my $x = $dbSaddleKar[POSCAR_COORDINATES][$a][0] - $rotPointX;
                        my $y = $dbSaddleKar[POSCAR_COORDINATES][$a][1] - $rotPointY;
                        my $z = $dbSaddleKar[POSCAR_COORDINATES][$a][2] - $rotPointZ;
                        my ($xp, $yp, $zp) = rotatePointAxis($x, $y, $z, $axisFirstRotation[0], $axisFirstRotation[1], $axisFirstRotation[2], -$thetaFirstRotation);
                        $dbSaddleKar[POSCAR_COORDINATES][$a][0] = $xp + $rotPointX;
                        $dbSaddleKar[POSCAR_COORDINATES][$a][1] = $yp + $rotPointY;
                        $dbSaddleKar[POSCAR_COORDINATES][$a][2] = $zp + $rotPointZ;
                    }

                    # Rotate dbSaddleKar according to the second rotation.
                    for(my $a = 0; $a < $dbNumAtoms; $a++)
                    {
                        my $x = $dbSaddleKar[POSCAR_COORDINATES][$a][0] - $rotPointX;
                        my $y = $dbSaddleKar[POSCAR_COORDINATES][$a][1] - $rotPointY;
                        my $z = $dbSaddleKar[POSCAR_COORDINATES][$a][2] - $rotPointZ;
                        my ($xp, $yp, $zp) = rotatePointAxis($x, $y, $z, $axisSecondRotation[0], $axisSecondRotation[1], $axisSecondRotation[2], -$thetaSecondRotation);
                        $dbSaddleKar[POSCAR_COORDINATES][$a][0] = $xp + $rotPointX;
                        $dbSaddleKar[POSCAR_COORDINATES][$a][1] = $yp + $rotPointY;
                        $dbSaddleKar[POSCAR_COORDINATES][$a][2] = $zp + $rotPointZ;
                    }

                    # Rotate and translate dbSaddleKar around db_a according to the stored steps, torques, and thetas.
                    for(my $t = 0; $t < @thetas; $t++)
                    {
                        for(my $a = 0; $a < $dbNumAtoms; $a++)
                        {
                            my $x = $dbSaddleKar[POSCAR_COORDINATES][$a][0] - $rotPointX;
                            my $y = $dbSaddleKar[POSCAR_COORDINATES][$a][1] - $rotPointY;
                            my $z = $dbSaddleKar[POSCAR_COORDINATES][$a][2] - $rotPointZ;
                            my ($xp, $yp, $zp) = rotatePointAxis($x, $y, $z, $torques[$t][0], $torques[$t][1], $torques[$t][2], $thetas[$t]);
                            $dbSaddleKar[POSCAR_COORDINATES][$a][0] = $xp + $rotPointX;
                            $dbSaddleKar[POSCAR_COORDINATES][$a][1] = $yp + $rotPointY;
                            $dbSaddleKar[POSCAR_COORDINATES][$a][2] = $zp + $rotPointZ;
                            $dbSaddleKar[POSCAR_COORDINATES][$a][0] += $steps[$t][0];
                            $dbSaddleKar[POSCAR_COORDINATES][$a][1] += $steps[$t][1];
                            $dbSaddleKar[POSCAR_COORDINATES][$a][2] += $steps[$t][2];
                        }
                    }
                
                    # Create a cartesian copy of $Car to be used as a result saddle poscar file.                        print "step: @step\n";

                    my @carSaddleKar = DirKar($Car);
                
                    # Set each of the saddle coordinates in the result saddle file equal to the transformed mobile coords.
                    for(my $j = 0; $j < $dbNumAtoms; $j++)
                    {
                        my $frozen = index($carSaddleKar[6][$maps{$j}], "F");
                        if($frozen < 0)
                        {
                            $carSaddleKar[POSCAR_COORDINATES][$maps{$j}][0] = $dbSaddleKar[POSCAR_COORDINATES][$j][0];
                            $carSaddleKar[POSCAR_COORDINATES][$maps{$j}][1] = $dbSaddleKar[POSCAR_COORDINATES][$j][1];
                            $carSaddleKar[POSCAR_COORDINATES][$maps{$j}][2] = $dbSaddleKar[POSCAR_COORDINATES][$j][2];
                        }
                    }
                    # Convert it back to direct coords.
                    my @result = KarDir(\@carSaddleKar);
                    if($dedupe)
                    {
                        # Make sure it's not duplicated in @saddleMatches.
                        my $isDuplicate = 0;
                        for(my $j = 0; $j < $numMatches; $j++)
                        {
                            for(my $k = 0; $k < $result[4]; $k++)
                            {
                                my $diff = pbc_difference([$result[0][$k]], [$saddleMatches[$j][$k]], 1);
                                dirkar($diff, $result[POSCAR_BASIS], $result[POSCAR_LATTICE], 1); 
                                my $dist = magnitude($diff, 1);
                                if($dist > 0.25) #TODO: Parameterize this value.
                                {
                                    $isDuplicate = 0;
                                    last;
                                }
                                $isDuplicate = 1;
                            }
                            if($isDuplicate)
                            {
                                last;   
                            }                
                        }
                        if($isDuplicate)
                        {
                            next;   
                        }                
                        else
                        {
                            push @saddleMatches, $result[0];
                        }
                    }
                    
                    # Save it to ./kdbmatches/SADDLE_$numMatches
                    write_poscar($result[POSCAR_COORDINATES], $result[POSCAR_BASIS], $result[POSCAR_LATTICE], 
                                 $result[POSCAR_NUM_ATOMS], $result[POSCAR_TOTAL_ATOMS], 
                                 $result[POSCAR_SELECTIVE_FLAG], $result[POSCAR_SELECTIVE], 
                                 $result[POSCAR_DESCRIPTION], "./kdbmatches/SADDLE_" . $numMatches);
                                
                    # Load the mode data.
                    my @mode = ();
                    for(my $j = 0; $j < $NumCarAtoms; $j++)
                    {
                        $mode[$j] = [0.0, 0.0, 0.0];
                    }
                    open(MODEFILE, "<$dir[0]/mode");
                    for(my $j = 0; $j < $dbNumAtoms; $j++)
                    {
                        my $line = <MODEFILE>;
                        chomp $line;
                        my @data = split(" ", $line);
                        for(my $k = 0; $k < $CarKar[POSCAR_TOTAL_ATOMS]; $k++)
                        {
                            if($maps{$j} == $k)
                            {
                                $mode[$maps{$j}] = [$data[0], $data[1], $data[2]];
                            }
                        }
                    }
                    close MODEFILE;

                    # Rotate the modes.
                    for(my $j = 0; $j < $dbNumAtoms; $j++)
                    {
                        my $x = $mode[$maps{$j}][0];
                        my $y = $mode[$maps{$j}][1];
                        my $z = $mode[$maps{$j}][2];
                        my ($xp, $yp, $zp) = rotatePointAxis($x, $y, $z, $axisFirstRotation[0], $axisFirstRotation[1], $axisFirstRotation[2], -$thetaFirstRotation);
                        $mode[$maps{$j}][0] = $xp;
                        $mode[$maps{$j}][1] = $yp;
                        $mode[$maps{$j}][2] = $zp;
                    }
                    for(my $j = 0; $j < $dbNumAtoms; $j++)
                    {
                        my $x = $mode[$maps{$j}][0];
                        my $y = $mode[$maps{$j}][1];
                        my $z = $mode[$maps{$j}][2];
                        my ($xp, $yp, $zp) = rotatePointAxis($x, $y, $z, $axisSecondRotation[0], $axisSecondRotation[1], $axisSecondRotation[2], -$thetaSecondRotation);
                        $mode[$maps{$j}][0] = $xp;
                        $mode[$maps{$j}][1] = $yp;
                        $mode[$maps{$j}][2] = $zp;
                    }
                    for(my $t = 0; $t < @thetas; $t++)
                    {
                        for(my $j = 0; $j < $dbNumAtoms; $j++)
                        {
                            my $x = $mode[$maps{$j}][0];
                            my $y = $mode[$maps{$j}][1];
                            my $z = $mode[$maps{$j}][2];
                            my ($xp, $yp, $zp) = rotatePointAxis($x, $y, $z, $torques[$t][0], $torques[$t][1], $torques[$t][2], $thetas[$t]);
                            $mode[$maps{$j}][0] = $xp;
                            $mode[$maps{$j}][1] = $yp;
                            $mode[$maps{$j}][2] = $zp;
                        }
                    }
                    # Normalize the rotated modecar.
                    my $total = 0;
                    for(my $j = 0; $j < @mode; $j++)
                    {
                        $total += $mode[$j][0] * $mode[$j][0];
                        $total += $mode[$j][1] * $mode[$j][1];
                        $total += $mode[$j][2] * $mode[$j][2];
                    }
                    $total = sqrt($total);
                    if($total != 0.0)
                    {
                        for(my $j = 0; $j < @mode; $j++)
                        {
                            $mode[$j][0] /= $total;
                            $mode[$j][1] /= $total;
                            $mode[$j][2] /= $total;
                        }
                    }
                    
                    # Save the modecar in ./kdbmatches/MODE_XXX
                    open(MODE, ">./kdbmatches/MODE_" . $numMatches);
                    for(my $j = 0; $j < $NumCarAtoms; $j++)
                    {
                        printf MODE "    % .16f    % .16f    % .16f\n", $mode[$j][0], $mode[$j][1], $mode[$j][2];
                    }
                    close MODE;

                    system("touch ./kdbmatches/.done_$numMatches");
                    
                    $numMatches += 1;

                }            
                
            } # End mirror pair.

        } # End mappings.

    } # End loop over process directories.

    print "\nDone.\n";
    
    
#    print "$debugCount\n";


#===================================================================================================
# max and min between two numbers. sqrt is built into perl, but not this??
#---------------------------------------------------------------------------------------------------
    sub max
    {
        my $a = shift;
        my $b = shift;
        if($a > $b)
        {
            return $a;
        }
        else
        {
            return $b;
        }
    }
#---------------------------------------------------------------------------------------------------

    sub min
    {
        my $a = shift;
        my $b = shift;
        if($a < $b)
        {
            return $a;
        }
        else
        {
            return $b;
        }
    }

#===================================================================================================
# Subroutine interAtomicVectorDirKar($data1, $data2, $a1, $a2, $reference).
#---------------------------------------------------------------------------------------------------
#   DESCRIPTION:    Returns a vector in full cartesian coordinates from atom $a1 in @data1 to atom $a2
#					in @data2 using @reference for the BASIS and LATTICE information.  This function
#					takes direct coordinates and returns cartesian coordinates. 
#
#   INPUT:          see DESCRIPTION
#
#   OUTPUT:         a 1x3 vector in full cartesian coordinates
#---------------------------------------------------------------------------------------------------
    sub interAtomicVectorDirKar
    {
        my ($data1, $data2, $a1, $a2, $reference) = @_;
        #my @data1 = @{dclone($data1)};
        #my @data2 = @{dclone($data2)};
        my @reference = @$reference;
        my @aCoords = $data1->[0][$a1];
        my @bCoords = $data2->[0][$a2];
        my @diff = pbc_difference(\@bCoords, \@aCoords, 1);
        dirkar(@diff, $reference[POSCAR_BASIS], $reference[POSCAR_LATTICE], 1);
		@diff = @{$diff[0][0]};
        return @diff;
    }

#===================================================================================================
# Subroutine interAtomicVectorKarKar($data1, $data2, $a1, $a2, $reference).
#---------------------------------------------------------------------------------------------------
#   DESCRIPTION:    Returns a vector in full cartesian coordinates from atom $a1 in @data1 to atom $a2
#					in @data2 using @reference for the BASIS and LATTICE information.  This function
#					takes cartesian coordinates and returns cartesian coordinates. 
#
#   INPUT:          see DESCRIPTION
#
#   OUTPUT:         a 1x3 vector in full cartesian coordinates
#---------------------------------------------------------------------------------------------------
    sub interAtomicVectorKarKar
    {
        my ($data1, $data2, $a1, $a2, $reference) = @_;
#        my @data1 = @{dclone($data1)};
#        my @data2 = @{dclone($data2)};
        my @reference = @$reference;
        my $aCoords = [[$data1->[0][$a1][0], $data1->[0][$a1][1], $data1->[0][$a1][2]]];
        my $bCoords = [[$data2->[0][$a2][0], $data2->[0][$a2][1], $data2->[0][$a2][2]]];
        kardir($aCoords, $reference[POSCAR_BASIS], $reference[POSCAR_LATTICE], 1);
        kardir($bCoords, $reference[POSCAR_BASIS], $reference[POSCAR_LATTICE], 1);
        my @diff = pbc_difference($bCoords, $aCoords, 1);
        dirkar(@diff, $reference[POSCAR_BASIS], $reference[POSCAR_LATTICE], 1);
		@diff = @{$diff[0][0]};
        return @diff;
    }


#===================================================================================================
# Subroutine interAtomicDistanceDirKar($data1, $data2, $a1, $a2, $reference).
#---------------------------------------------------------------------------------------------------
#   DESCRIPTION:    Returns a distance in full cartesian units from atom $a1 in @data1 to atom $a2
#					in @data2 using @reference for the BASIS and LATTICE information.  This function
#					takes direct coordinates and returns cartesian units. 
#
#   INPUT:          see DESCRIPTION
#
#   OUTPUT:         a scalar representing the distance between the two atoms in full cartesian units
#---------------------------------------------------------------------------------------------------
    sub interAtomicDistanceDirKar
    {
        my ($data1, $data2, $a1, $a2, $reference) = @_;
        #my @data1 = @{dclone($data1)};
        #my @data2 = @{dclone($data2)};
        my @aCoords = $data1->[0][$a1];
        my @bCoords = $data2->[0][$a2];
        my @diff = pbc_difference(\@bCoords, \@aCoords, 1);
        dirkar(@diff, $reference->[POSCAR_BASIS], $reference->[POSCAR_LATTICE], 1);
        my $x = $diff[0][0][0];
        my $y = $diff[0][0][1];
        my $z = $diff[0][0][2];
        return sqrt($x * $x + $y * $y + $z * $z);
    }

#===================================================================================================
# Subroutine intraAtomicVectorDirKar($data, $a1, $a2).
#---------------------------------------------------------------------------------------------------
#   DESCRIPTION:    Returns a vector in full cartesian coordinates from atom $a1 in to atom $a2
#					in @data. This function	takes direct coordinates and returns cartesian coordinates. 
#
#   INPUT:          see DESCRIPTION
#
#   OUTPUT:         a 1x3 vector in full cartesian coordinates
#---------------------------------------------------------------------------------------------------
    sub intraAtomicVectorDirKar
    {
        my ($data, $a1, $a2) = @_;
        #my @data = @{dclone($data)};
        my @aCoords = $data->[0][$a1];
        my @bCoords = $data->[0][$a2];
        my @diff = pbc_difference(\@bCoords, \@aCoords, 1);
        dirkar(@diff, $data->[POSCAR_BASIS], $data->[POSCAR_LATTICE], 1);
		@diff = @{$diff[0][0]};
        return @diff;
    }



