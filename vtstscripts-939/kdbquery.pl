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
use Carp 'verbose';
$SIG{'SIGINT'} = sub { Carp::confess( @_ ) };

#===================================================================================================
# Global Constants
#---------------------------------------------------------------------------------------------------

    my $debug = '';
    my $BarrierCutoff = 3.0;
    my $ScoreCutoff = 0.25;
    my $RHScoreCutoff = 0.5;
    my $RHCutoff = 5.0;
    my $NeighborFudgePercent = 0.3;
    my $RHBinSize = 5.0;
    my $RHScale = 1.0;
    my $MAPPING = '';
    GetOptions("debug"=>\$debug, 
               "mapping"=>\$MAPPING,
               "score=f", \$ScoreCutoff, 
               "rhsco=f", \$RHScoreCutoff, 
               "rhc=f", \$RHCutoff, 
               "rhbs=f", \$RHBinSize, 
               "rhs=f", \$RHScale, 
               "nfp=f", \$NeighborFudgePercent);
        
    my $debugCount = 0;
    
    
    # Length in Angstroms within which to say two distances are equal.
    use constant ANGSTROM_FUDGE => 0.25;

    # The number of bins in the radial distribution function, given the cutoff and bin size.
    my $RHNumBins = $RHCutoff / $RHBinSize;
    
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
    
    # The maximum number of iterations to minimize a torque.
    my $maxIter = 16;

    # The spring constant used for torque minimization.
    my $springK = 1.0;    
    
    
    
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

    # Store the list of elements inside $Desc.
    my $Desc = $Car->[POSCAR_DESCRIPTION];

    my @CarElements = CarElementTable($Car);

    # Create a copy of $Car in cartesian format.
    my @CarKar = DirKar($Car);

#===================================================================================================
# Main algorithm.
#---------------------------------------------------------------------------------------------------

    # Create the radial histograms for the POSCAR we are trying to match. 
    my @carRHs = ();
    for(my $atom = 0; $atom < $NumCarAtoms; $atom++)
    {
        $carRHs[$atom] = getRadialHistogram($atom, $Car, \@CarDistanceTable);
    }
    
    
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

    my $numTriads = 0;
    my $numDistedTriads = 0;
    my $numGoodTriads = 0;

    # Loop over each directory.
    my $dirIndex = 0;
    foreach my $dirPointer (@AllDirs)
    {

        # Dereference $dirPointer into @dir.
        my @dir = @{$dirPointer};
        
        # Notify the user we are checking the minimum.
        print "Checking $dir[1]";

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
                $minMatches[$i][$j] = matchDescription($minElements[$i], $CarElements[$j]);
            }
        }

        # Create a distance table for the database minimum.
        my @dbMinDistanceTable = atomicDistanceTableDir(\@dbMinCar);
        
        my $smallestMatchCount = 1e300;
        my @q_as = ();
        my $db_a = $dbMinKar[8]{$mobile[0]};
        for(my $mob = 0; $mob < @mobile; $mob++)
        {
            # Create the RDFs for the mobile atoms of the database minimum.
            my $db_aRH = getRadialHistogram($dbMinKar[8]{$mobile[$mob]}, \@dbMinCar, \@dbMinDistanceTable);


            # Create a list of scores for each rh.
            my @minScores = ();
            for(my $b = 0; $b < $NumCarAtoms; $b++)
            {
                $minScores[$b] = rh_score($db_aRH, $carRHs[$b]);

            }
            
            my @carAtoms = ();
            
            # Get the atoms in Car which best match $db_a.
            my $elementA = nthAtomElement($dbMinKar[8]{$mobile[$mob]}, \@dbMinKar);
            for(my $b = 0; $b < $NumCarAtoms; $b++)
            {
                if (index($Car->[POSCAR_SELECTIVE][$b], "F") != -1)
                {
                    next;
                }
                my $elementB = nthAtomElement($b, $Car);
                if(matchDescription($elementA, $elementB))
                {
                    if($minScores[$b] < $RHScoreCutoff)
                    {
                        push @carAtoms, $b;
                    }
                }
            }    
            
            if(@carAtoms < $smallestMatchCount)
            {
                @q_as = @carAtoms;
                $smallestMatchCount = @carAtoms;
                $db_a = $dbMinKar[8]{$mobile[$mob]};
            }
        }
        
        # Calculate the radius of the db minimum as the distance from db_a to the 
        # atom furthest from it.
        my $db_radius = 0.0;
        for(my $i = 0; $i < $dbNumAtoms; $i++)
        {
            if($i != $db_a)
            {
                my $d = $dbMinDistanceTable[$db_a][$i];
                if($d > $db_radius)
                {
                    $db_radius = $d;
                }
            }
        }

        if(@q_as == 0)
        {
            if($debug)
            {
                print "    No good Car atoms, nexting.\n";
            }
            if(!$debug)
            {
                print "\n";
            }
            next;
        }
        

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
        my $dbPairRadius = 99999999999;
        my @db_ab = intraAtomicVectorDirKar(\@dbMinCar, $db_a, $db_b); 
        my @db_ab_u = unitVector(\@db_ab);
        for(my $i = 0; $i < $dbNumAtoms; $i++)
        {
            if($i != $db_a && $i != $db_b)
            {
                my $d = $dbMinDistanceTable[$db_a][$i];
                if($d < $dbPairRadius)
                {
                    my $c = $i;
                    my @db_ac = intraAtomicVectorDirKar(\@dbMinCar, $db_a, $c);        
                    my @db_ac_u = unitVector(\@db_ac);
                    my $dp = abs(dotProduct(\@db_ab_u, \@db_ac_u));
                    # Make sure the three atoms are not in a line.
                    if($dp < 0.9)
                    {
                        $db_c = $i;
                        $dbPairRadius = $d;
                    }
                }
            }
        }
        
        
        if($debug)
        {
            print "\n    db_abc: $db_a, $db_b, $db_c\n";
        }        
        
        if($debug)
        {
            my $foo = @q_as;
            print "    q_as: @q_as\n";
        }
        
        
        if($MAPPING)
        {}
            
        
        
        QA_LOOP: for my $q_a(@q_as)
        {


######### Make a list of possible mappings of neighbors around db_a and q_a. ############################################
my @mappings;
if($MAPPING)
{
            # Make a list of nearest neighbors around db_a and q_a
            my @dbNeighbors = ();
            my @qNeighbors = ();
            
            my $dbType = nthAtomElement($db_a, \@dbMinCar);
            my $dbd = $covalentRadii{$dbType};
            for(my $i = 0; $i < $dbNumAtoms; $i++)
            {
                if ($i == $db_a)
                {
                    next;
                }
                
                my $d = $dbMinDistanceTable[$db_a][$i];
                my $iType = nthAtomElement($i, \@dbMinCar);
                my $id = $covalentRadii{$iType};
                
                if ($d < ($id + $dbd) * (1.0 + $NeighborFudgePercent))
                {
                    push @dbNeighbors, $i;
                }
            }

            my $qType = nthAtomElement($q_a, $Car);
            my $qd = $covalentRadii{$qType};
            for(my $i = 0; $i < $NumCarAtoms; $i++)
            {
                if ($i == $q_a)
                {
                    next;
                }
                
                my $d = $CarDistanceTable[$q_a][$i];
                my $iType = nthAtomElement($i, $Car);
                my $id = $covalentRadii{$iType};
                
                if ($d < ($id + $dbd) * (1.1 + $NeighborFudgePercent)) #THIS IS A LITTLE MORE LOOSE THAN ABOVE.
                {
                    push @qNeighbors, $i;
                }
            }

            @mappings = ({$db_a => $q_a});
            my @newMappings = ();
    
            foreach my $dbz(@dbNeighbors)
            {
                foreach my $qz(@qNeighbors)
                {
                    MAPPING: for (my $i = 0; $i < @mappings; $i++)
                    {
                        my %hash = %{$mappings[$i]};
                        
                        foreach my $key(keys %hash)
                        {
                            if(abs($dbMinDistanceTable[$key][$dbz] - $CarDistanceTable[$hash{$key}][$qz]) > ANGSTROM_FUDGE)
                            {
                                next MAPPING;
                            }
                        }
                        
                        my %map = %hash;
                        $map{$dbz} = $qz;
                        push @newMappings, \%map;
                    }                
                }
                
                @mappings = @newMappings;
                @newMappings = ();
            }

            my $l = @mappings;
            print " $l";

            if ($l < 1)
            {
                next QA_LOOP;
            }    
}            
#########################################################################################################################


            if(!$debug)
            {
                printf(".");
            }
            
            if($debug) 
            {
                print "    q_a: $q_a\n";
            }

            # Create a list of Car atoms within $dbPairRadius + ANGSTROM_FUDGE of $q_a.
            my @carPairAtoms = ();
            for(my $i = 0; $i < $NumCarAtoms; $i++)
            {
                if($i != $q_a)
                {
                    if($CarDistanceTable[$q_a][$i] < ANGSTROM_FUDGE + $dbPairRadius)
                    {
                        push(@carPairAtoms, $i)
                    }
                }
            }
            
            
            # Create a list of Car atoms within $db_radius + ANGSTROM_FUDGE of $q_a.
            my @carRadiusAtoms = ();
            for(my $i = 0; $i < $NumCarAtoms; $i++)
            {
                if($CarDistanceTable[$q_a][$i] < ANGSTROM_FUDGE + $db_radius)
                {
                    push @carRadiusAtoms, $i;
                }
            }
            
            
            if($debug) 
            {
                print "        carPairAtoms: @carPairAtoms\n";
            }

            
            # Get a list of viable pairs of atoms.
            my @carAtomPairs = ();
            my $IA = $dbMinDistanceTable[$db_a][$db_b];
            my $IB = $dbMinDistanceTable[$db_a][$db_c];
            my $AB = $dbMinDistanceTable[$db_b][$db_c];
            for(my $i = 0; $i < @carPairAtoms; $i++)
            {
                for(my $j = 0; $j < @carPairAtoms; $j++)
                {
                    if($i != $j)
                    {
                        if(abs($CarDistanceTable[$carPairAtoms[$i]][$carPairAtoms[$j]] - $AB) < ANGSTROM_FUDGE)
                        {
                            if(abs($CarDistanceTable[$q_a][$carPairAtoms[$i]] - $IA) < ANGSTROM_FUDGE)
                            {
                                if(abs($CarDistanceTable[$q_a][$carPairAtoms[$j]] - $IB) < ANGSTROM_FUDGE)
                                {
                                    push(@carAtomPairs, [$carPairAtoms[$i], $carPairAtoms[$j]]);
                                }
                            }
                        }
                    }
                }   
            }


#########################################################################################################################
if($MAPPING)
{            
            @carAtomPairs = ();
            for (my $i = 0; $i < @mappings; $i++)
            {
                my %hash = %{$mappings[$i]};
                push(@carAtomPairs, [$hash{$db_b}, $hash{$db_c}]);
            }
}
#########################################################################################################################
            
            if($debug)
            {
                my $foo = @carAtomPairs;
                print "        carAtomPairs: $foo\n";
            }
            
            
            
            LOOP_ATOM_PAIR: for(my $atomPair = 0; $atomPair < @carAtomPairs; $atomPair++)
            {                 
            
                $numTriads += 1;
            
                if($debug)
                {
                    print "        atomPair: $carAtomPairs[$atomPair]->[0], $carAtomPairs[$atomPair]->[1]\n"
                }
                
                my @localdbMinKar = @{dclone(\@dbMinKar)};
                my @localdbMinCar = @{dclone(\@dbMinCar)};
                my $q_b = $carAtomPairs[$atomPair]->[0];
                my $q_c = $carAtomPairs[$atomPair]->[1];
                my @translation = (0, 0, 0);
                
                # Find a translation for the minimum.
		        @translation = interAtomicVectorDirKar(\@localdbMinCar, $Car, $db_a, $q_a, $Car);
		        
		        if($debug)
		        {
		            print "            translation: @translation\n";
		        }                        


                # Shift the new coords by @translation.
                for(my $a = 0; $a < $localdbMinKar[POSCAR_TOTAL_ATOMS]; $a++)
                {
                    $localdbMinKar[POSCAR_COORDINATES][$a][0] += $translation[0];
                    $localdbMinKar[POSCAR_COORDINATES][$a][1] += $translation[1];
                    $localdbMinKar[POSCAR_COORDINATES][$a][2] += $translation[2];
                }
                
                #@localdbMinKar = @{setPBCKarKar(\@localdbMinKar, $Car)};
                #@CarKar = @{setPBCKarKar(\@CarKar, $Car)};


                # If $debugging,start the .xyz file.
                if($debug)
                {
                    writeXYZ(\@localdbMinKar, "$dirIndex-$q_b-$q_c.xyz");
                }
                

                
                # Determine the unit vector from the best db atom to the second db best atom.
                my @bestDBUnit = intraAtomicVectorDirKar(\@localdbMinCar, $db_a, $db_b);
                @bestDBUnit = unitVector(\@bestDBUnit);
                # Determine the unit vector from the best car atom to the second best car atom.
                my @bestCarUnit = intraAtomicVectorDirKar($Car, $q_a, $q_b);
                @bestCarUnit = unitVector(\@bestCarUnit);

                # Calculate the rotation axis and angle needed to rotate bestDBUnit onto bestCarUnit
                my @crossProduct = crossProduct(\@bestDBUnit, \@bestCarUnit);
                
                # If the cross product is the zero vector, then the vectors are the same vector and we don't need to do the first rotation.
                my @axisFirstRotation = ();
                my $thetaFirstRotation = 0;
                if ($crossProduct[0] != 0 or $crossProduct[1] != 0 or $crossProduct[2] != 0)
                {
                    @axisFirstRotation = unitVector(\@crossProduct);
                    my $dp = dotProduct(\@bestDBUnit, \@bestCarUnit);
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

                # Rotate localdbMinKar so bestDBUnit aligns with bestCarUnit.
                for(my $a = 0; $a < $localdbMinKar[POSCAR_TOTAL_ATOMS]; $a++)
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
                # Calculate the angle for matching the third-best-matched atom.
                my @q_ab = intraAtomicVectorDirKar($Car, $q_a, $q_b);        
                my @q_ab_u = unitVector(\@q_ab);
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
                    for(my $a = 0; $a < $localdbMinKar[POSCAR_TOTAL_ATOMS]; $a++)
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

                if($debug)
                {
                    appendXYZ(\@localdbMinKar, "$dirIndex-$q_b-$q_c.xyz");
                }
                
                
                if($debug && 1)
                {
                    print "            Creating distance table...\n";
                }
                
                #my $distType = "brute";
                my $distType = "quick";

                my %takenDB;

                $numDistedTriads += 1;

                if ($distType eq "brute")
                {
                    # Create a table of distances from each db atom to atom in @CarKar.
                    my @allDistances = ();
                    my @costTable;
                    my @minCar = @{dclone(\@localdbMinKar)};
                    $minCar[POSCAR_BASIS] = $Car->[POSCAR_BASIS];
                    $minCar[POSCAR_LATTICE] = $Car->[POSCAR_LATTICE];
                    kardir($minCar[0], $minCar[POSCAR_BASIS], $minCar[POSCAR_LATTICE], $minCar[POSCAR_TOTAL_ATOMS]);
                    for(my $i = 0; $i < $minCar[POSCAR_TOTAL_ATOMS]; $i++)
                    {
                        my $bestDistance = 777.0;
                        foreach(@carRadiusAtoms)
                        {
                            if($minMatches[$i][$_])
                            {
                                my $dist = interAtomicDistanceDirKar(\@minCar, $Car, $i, $_, $Car);
                                if ($dist < $bestDistance)
                                {
                                    $bestDistance = $dist;
                                }
                                $costTable[$i][$_] = $dist;
                                push @allDistances, [$i, $_, $dist];
                            }
                        }
                        # If we don't have a good distance for this atom, skip to the next pair.
                        if ($bestDistance > 0.7)
                        {
                            next LOOP_ATOM_PAIR;
                        }
                    }
                    
                    if($debug && 1)
                    {
                        print "            Determining best distances...\n";
                    }
                    
                    # Get the best distance for each db atom.
                    my %takenCar;
                    my @sorted = sort {$a->[2] <=> $b->[2]} @allDistances;
                    my $matchCount = 0;
                    my $index = 0;
                    while($matchCount < $localdbMinCar[POSCAR_TOTAL_ATOMS])
                    {
                        # If we don't have a pretty good match for all atoms, skip to the next pair. 
                        if($sorted[$index][2] > 0.7)
                        {
                            next LOOP_ATOM_PAIR;
                        }
                        my $a = $sorted[$index][0];
                        my $b = $sorted[$index][1];
                        if(!exists $takenDB{$a} && !exists $takenCar{$b})
                        {
                            $takenDB{$a} = $b;
                            $takenCar{$b} = $a;
                            $matchCount++;
                        }
                        $index++;
                    }
                }
                
                
                
                # Create a naive mapping from each atom in minCar to each atom in CarKar
                if ($distType eq "quick")
                {

                    my @localdbMinCar = @{dclone(\@localdbMinKar)};
                    $localdbMinCar[POSCAR_BASIS] = $Car->[POSCAR_BASIS];
                    $localdbMinCar[POSCAR_LATTICE] = $Car->[POSCAR_LATTICE];
                    kardir($localdbMinCar[0], $localdbMinCar[POSCAR_BASIS], $localdbMinCar[POSCAR_LATTICE], $localdbMinCar[POSCAR_TOTAL_ATOMS]);
                    my %localCarRadiusAtoms = ();
                    foreach(@carRadiusAtoms)
                    {
                        $localCarRadiusAtoms{$_} = 1;
                    }
                    for(my $i = 0; $i < $localdbMinCar[POSCAR_TOTAL_ATOMS]; $i++)
                    {
                        my $bestDistance = 777.0;
                        my $bestID = -1;
                        for my $j (keys %localCarRadiusAtoms)
                        {
                            if($minMatches[$i][$j])
                            {
                                my $dist = interAtomicDistanceDirKar(\@localdbMinCar, $Car, $i, $j, $Car);
                                if ($dist < $bestDistance)
                                {
                                    $bestDistance = $dist;
                                    $bestID = $j;
                                    if($bestDistance < 0.1) #TODO: Parameterize.
                                    {
                                        last;
                                    }
                                }
                            }
                        }
                        # If we don't have a good distance for this atom, skip to the next pair.
                        if ($bestDistance > 0.7) #TODO: Parameterize.
                        {
                            next LOOP_ATOM_PAIR;
                        }
                        $takenDB{$i} = $bestID;
                        delete $localCarRadiusAtoms{$bestID};
                    }
                }


                if($debug && 1)
                {
                    print "            Minimizing torque...\n";
                }
                
                            
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
                while($torqueMagnitude + $forceMagnitude > 0.00001 && $iterCount < $maxIter)
                {
                    $iterCount++;
                    my @torque = (0, 0, 0);
                    my @force = (0, 0, 0);
                    for(my $a = 0; $a < $localdbMinKar[POSCAR_TOTAL_ATOMS]; $a++)
                    {
                        my($dx, $dy, $dz) = interAtomicVectorKarKar(\@CarKar, \@localdbMinKar, $takenDB{$a}, $a, $Car);
                        $debugCount += 1;
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
                    
                    if($debug && 0)
                    {
                        print "                Torque magnitude: $torqueMagnitude\n";
                    }
                    my $theta = $tdt * -sqrt($torqueVelocity[0]*$torqueVelocity[0] + $torqueVelocity[1]*$torqueVelocity[1] + $torqueVelocity[2]*$torqueVelocity[2]);
                    my $thetaMax = 0.01;
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
                        $step[0] = $stepMax * $step[0]/$mag;
                        $step[1] = $stepMax * $step[1]/$mag;
                        $step[2] = $stepMax * $step[2]/$mag;
                    }
                    for(my $a = 0; $a < $dbMinKar[POSCAR_TOTAL_ATOMS]; $a++)
                    {
                        my $x = $localdbMinKar[POSCAR_COORDINATES][$a][0] - $localdbMinKar[POSCAR_COORDINATES][$db_a][0];
                        my $y = $localdbMinKar[POSCAR_COORDINATES][$a][1] - $localdbMinKar[POSCAR_COORDINATES][$db_a][1];
                        my $z = $localdbMinKar[POSCAR_COORDINATES][$a][2] - $localdbMinKar[POSCAR_COORDINATES][$db_a][2];
                        my ($xp, $yp, $zp) = rotatePointAxis($x, $y, $z, $torque[0], $torque[1], $torque[2], $theta);
                        $localdbMinKar[POSCAR_COORDINATES][$a][0] = $xp + $localdbMinKar[POSCAR_COORDINATES][$db_a][0];
                        $localdbMinKar[POSCAR_COORDINATES][$a][1] = $yp + $localdbMinKar[POSCAR_COORDINATES][$db_a][1];
                        $localdbMinKar[POSCAR_COORDINATES][$a][2] = $zp + $localdbMinKar[POSCAR_COORDINATES][$db_a][2];
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

                if($debug && 1)
                {
                    print "            Creating distance table...\n";
                }

                # Create a table of distances from each db atom to atom in @CarKar.
                my @allDistances = ();
                my @costTable = ();
                my @minCar = @{dclone(\@localdbMinKar)};
                $minCar[POSCAR_BASIS] = $Car->[POSCAR_BASIS];
                $minCar[POSCAR_LATTICE] = $Car->[POSCAR_LATTICE];
                kardir($minCar[0], $minCar[POSCAR_BASIS], $minCar[POSCAR_LATTICE], $minCar[POSCAR_TOTAL_ATOMS]);
                foreach my $a(keys %takenDB)
                {
                    my $dist = interAtomicDistanceDirKar(\@minCar, $Car, $a, $takenDB{$a}, $Car);
                    $costTable[$a][$takenDB{$a}] = $dist;
                }

                if($debug && 0)
                {
                    print "            Calculating score...\n";
                }
                
                # Calculate the total score for this angle.
                my $score = 0;
                for(my $i = 0; $i < $localdbMinCar[POSCAR_TOTAL_ATOMS]; $i++)
                {
                    $score = max($score, $costTable[$i][$takenDB{$i}]);
                }

                if($debug)
                {
                    print "            score: $score\n";
                }

                if($score < $ScoreCutoff)
                {
                
                    $numGoodTriads += 1;
                
                    # Load the db saddle and mobile atoms.
                    my @dbSaddleKar = loadXYZ2CAR($dir[0] . "/saddle.xyz");

                    # Translate the mobile atoms.
                    for(my $j = 0; $j < $dbSaddleKar[POSCAR_TOTAL_ATOMS]; $j++)
                    {
                        $dbSaddleKar[POSCAR_COORDINATES][$j][0] += $translation[0];
                        $dbSaddleKar[POSCAR_COORDINATES][$j][1] += $translation[1];
                        $dbSaddleKar[POSCAR_COORDINATES][$j][2] += $translation[2];
                    }
                
                    my $rotPointX = $dbMinKar[POSCAR_COORDINATES][$db_a][0] + $translation[0];
                    my $rotPointY = $dbMinKar[POSCAR_COORDINATES][$db_a][1] + $translation[1];
                    my $rotPointZ = $dbMinKar[POSCAR_COORDINATES][$db_a][2] + $translation[2];

                    # Rotate dbSaddleKar according to the first rotation.
                    for(my $a = 0; $a < $dbSaddleKar[POSCAR_TOTAL_ATOMS]; $a++)
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
                    for(my $a = 0; $a < $dbSaddleKar[POSCAR_TOTAL_ATOMS]; $a++)
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
                        for(my $a = 0; $a < $dbSaddleKar[POSCAR_TOTAL_ATOMS]; $a++)
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
                    for(my $j = 0; $j < $dbSaddleKar[POSCAR_TOTAL_ATOMS]; $j++)
                    {
                        my $frozen = index($carSaddleKar[6][$takenDB{$j}], "F");
                        if($frozen < 0)
                        {
                            $carSaddleKar[POSCAR_COORDINATES][$takenDB{$j}][0] = $dbSaddleKar[POSCAR_COORDINATES][$j][0];
                            $carSaddleKar[POSCAR_COORDINATES][$takenDB{$j}][1] = $dbSaddleKar[POSCAR_COORDINATES][$j][1];
                            $carSaddleKar[POSCAR_COORDINATES][$takenDB{$j}][2] = $dbSaddleKar[POSCAR_COORDINATES][$j][2];
                        }
                    }
                    # Convert it back to direct coords.
                    my @result = KarDir(\@carSaddleKar);
                
#                    # Make sure it's not duplicated in @saddleMatches.
#                    my $isDuplicate = 0;
#                    for(my $j = 0; $j < $numMatches; $j++)
#                    {
#                        my $diff = pbc_difference($result[0], $saddleMatches[$j], $result[4]);
#                        dirkar($diff, $result[POSCAR_BASIS], $result[POSCAR_LATTICE], $result[POSCAR_TOTAL_ATOMS]); 
#                        my $dist = magnitude($diff, $result[4]);
#                        if($dist < 0.25) #TODO: Parameterize this value.
#                        {
#                            $isDuplicate = 1;
#                            last;
#                        }
#                    }
#                    if($isDuplicate)
#                    {
#                        next;   
#                    }                
#                    else
#                    {
#                        push @saddleMatches, $result[0];
#                    }
                    
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
                    for(my $j = 0; $j < $dbSaddleKar[POSCAR_TOTAL_ATOMS]; $j++)
                    {
                        my $line = <MODEFILE>;
                        chomp $line;
                        my @data = split(" ", $line);
                        for(my $k = 0; $k < $CarKar[POSCAR_TOTAL_ATOMS]; $k++)
                        {
                            if($takenDB{$j} == $k)
                            {
                                $mode[$takenDB{$j}] = [$data[0], $data[1], $data[2]];
                            }
                        }
                    }
                    close MODEFILE;

                    # Rotate the modes.
                    for(my $j = 0; $j < $dbSaddleKar[POSCAR_TOTAL_ATOMS]; $j++)
                    {
                        my $x = $mode[$takenDB{$j}][0];
                        my $y = $mode[$takenDB{$j}][1];
                        my $z = $mode[$takenDB{$j}][2];
                        my ($xp, $yp, $zp) = rotatePointAxis($x, $y, $z, $axisFirstRotation[0], $axisFirstRotation[1], $axisFirstRotation[2], -$thetaFirstRotation);
                        $mode[$takenDB{$j}][0] = $xp;
                        $mode[$takenDB{$j}][1] = $yp;
                        $mode[$takenDB{$j}][2] = $zp;
                    }
                    for(my $j = 0; $j < $dbSaddleKar[POSCAR_TOTAL_ATOMS]; $j++)
                    {
                        my $x = $mode[$takenDB{$j}][0];
                        my $y = $mode[$takenDB{$j}][1];
                        my $z = $mode[$takenDB{$j}][2];
                        my ($xp, $yp, $zp) = rotatePointAxis($x, $y, $z, $axisSecondRotation[0], $axisSecondRotation[1], $axisSecondRotation[2], -$thetaSecondRotation);
                        $mode[$takenDB{$j}][0] = $xp;
                        $mode[$takenDB{$j}][1] = $yp;
                        $mode[$takenDB{$j}][2] = $zp;
                    }
                    for(my $t = 0; $t < @thetas; $t++)
                    {
                        for(my $j = 0; $j < $dbSaddleKar[POSCAR_TOTAL_ATOMS]; $j++)
                        {
                            my $x = $mode[$takenDB{$j}][0];
                            my $y = $mode[$takenDB{$j}][1];
                            my $z = $mode[$takenDB{$j}][2];
                            my ($xp, $yp, $zp) = rotatePointAxis($x, $y, $z, $torques[$t][0], $torques[$t][1], $torques[$t][2], $thetas[$t]);
                            $mode[$takenDB{$j}][0] = $xp;
                            $mode[$takenDB{$j}][1] = $yp;
                            $mode[$takenDB{$j}][2] = $zp;
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

            } # End loop over carAtomPairs.

        } # End loop over @q_as.        

        $dirIndex++;    
        print "\n";

    } # End loop over process directories.

    print "$numTriads, $numDistedTriads, $numGoodTriads\n";

    print "Done.\n\n";
    
    
#    print "$debugCount\n";
    
#===================================================================================================
# Subroutine getRadialHistogram
#---------------------------------------------------------------------------------------------------

    sub getRadialHistogram
    {
        my ($atomid, $configuration, $distanceTable) = @_;
        my @rh = ();

        # Zero out the bins.
        for(my $element = 0; $element < 23; $element++)
        {
            for(my $bin = 0; $bin < $RHNumBins; $bin++)
            {
                $rh[$element][$bin] = 0;
            }
        }

        # Get the center atom element.
        my $centerElement = nthAtomElement($atomid, $configuration);
        
        for(my $i = 0; $i < $configuration->[POSCAR_TOTAL_ATOMS]; $i++)
        {
            if($i != $atomid)
            {
                my $atomElement = nthAtomElement($i, $configuration);   
                my $d = $distanceTable->[$atomid][$i];
                if($d < ($covalentRadii{$centerElement} + $covalentRadii{$atomElement}) * (1.0 + $NeighborFudgePercent))
                {

                    my $element = 0;
                    if(exists($NON_METALS{$atomElement}))
                    {
                        $element = $NON_METALS{$atomElement};
                    }
                    for(my $bin = 0; $bin < $RHNumBins; $bin++)
                    {
                        # Determine the left and right limits of this bin.
                        my $a = $bin * $RHBinSize;
                        my $b = ($bin + 1.0) * $RHBinSize;
                        # Increment the bin by the gaussian-weighted amount.
                        $rh[$element][$bin] += rh_integral($a, $b, $d);
                    }                     
                }
            }
        }
        return \@rh;            
    }
    
#===================================================================================================
# Subroutine rh_integral($a, $b, $center)
#---------------------------------------------------------------------------------------------------
#   DESCRIPTION:    Returns the value of the definite integral between $a and $b of the 
#                   gaussian-like function centered at $center.
#
#   INPUT:          $a: The limts for which we are evaluating the gaussian-like function.
#                   $center: The center point of this particular gaussian-like function.
#
#   OUTPUT:         The value of the definite integral between $a and $b.
#
#---------------------------------------------------------------------------------------------------
    sub rh_integral
    {
        my $a = shift;
        my $b = shift;
        my $center = shift;
        return (0.5 * tanh($RHScale * ($b - $center))) - (0.5 * tanh($RHScale * ($a - $center)));
    }

#===================================================================================================
# Subroutine rh_score($rh1, $rh2)
#---------------------------------------------------------------------------------------------------
#   DESCRIPTION:    Returns the distance, or score, between two rh's. Lower score is equivalent to
#                   a better match.
#
#   INPUT:          $rh1: The first rh. (from kdb)
#                   $rh2: The second rh. (the query)
#
#   OUTPUT:         The distance, or score, between rh1 and rh2.
#
#---------------------------------------------------------------------------------------------------
    sub rh_score
    {
        my $rh1Ref = shift;
        my $rh2Ref = shift;
        my @rh1 = @$rh1Ref;
        my @rh2 = @$rh2Ref;
        my $sum = 0;
        my $diff;
        # 23 is the number of nonmetals + 1 for metals.
        for(my $element = 0; $element < 23; $element++)
        {
            for(my $bin = 0; $bin < $RHNumBins; $bin++)
            {
                $diff = $rh1[$element][$bin] - $rh2[$element][$bin];
                $sum += abs($diff);# * $diff;
            }
        }
        return $sum;
    }


#===================================================================================================
# Subroutine rh_score_2($rh1, $rh2)
#---------------------------------------------------------------------------------------------------
#   DESCRIPTION:    Returns the distance, or score, between two rdfs. Lower score is equivalent to
#                   a better match. Extra atoms in system 2 do not count against the score.
#
#   INPUT:          $rh1: The first rh. (from kdb)
#                   $rh2: The second rh. (the query)
#
#   OUTPUT:         The distance, or score, between rh1 and rh2.
#
#---------------------------------------------------------------------------------------------------
    sub rh_score_2
    {
        my $rh1Ref = shift;
        my $rh2Ref = shift;
        my @rh1 = @$rh1Ref;
        my @rh2 = @$rh2Ref;
        my $sum = 0;
        my $diff;
        # 23 is the number of nonmetals + 1 for metals.
        for(my $element = 0; $element < 23; $element++)
        {
            for(my $bin = 0; $bin < $RHNumBins; $bin++)
            {
                $diff = $rh1[$element][$bin] - $rh2[$element][$bin];
                $sum += max(0, $diff);
            }
        }
        return $sum;
    }


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
# Subroutine setPBCKarKar($data, $reference).
#---------------------------------------------------------------------------------------------------
#   DESCRIPTION:    Returns a poscar structure with PBCs applied to $data according to BASIS and 
#                   LATTICE in $reference.  Takes and returns cartesian coordinates.
#
#   INPUT:          see DESCRIPTION
#
#   OUTPUT:         a reference to the PBC'd structure
#---------------------------------------------------------------------------------------------------
    sub setPBCKarKar
    {
        my($data, $reference) = @_;
        $data = dclone($data);
        kardir($data->[POSCAR_COORDINATES], $reference->[POSCAR_BASIS], $reference->[POSCAR_LATTICE], $data->[POSCAR_TOTAL_ATOMS]);
        set_bc($data->[POSCAR_COORDINATES], $data->[POSCAR_TOTAL_ATOMS]);
        dirkar($data->[POSCAR_COORDINATES], $reference->[POSCAR_BASIS], $reference->[POSCAR_LATTICE], $data->[POSCAR_TOTAL_ATOMS]);
        return $data;
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
# Subroutine coordString($data, $a)
#---------------------------------------------------------------------------------------------------
#   DESCRIPTION:    returns the coordinates of atom $a in configuration $data in a string.
#
#---------------------------------------------------------------------------------------------------
    sub coordString
    {
        my ($data, $a) = @_;
        my @aCoords = $data->[0][$a];
        return "($aCoords[0][0], $aCoords[0][1], $aCoords[0][2])";
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
# Subroutine interAtomicDistanceKarKar($data1, $data2, $a1, $a2, $reference).
#---------------------------------------------------------------------------------------------------
#   DESCRIPTION:    Returns a distance in full cartesian units from atom $a1 in @data1 to atom $a2
#					in @data2 using @reference for the BASIS and LATTICE information.  This function
#					takes cartesian coordinates and returns cartesian units. 
#
#   INPUT:          see DESCRIPTION
#
#   OUTPUT:         a scalar representing the distance between the two atoms in full cartesian units
#---------------------------------------------------------------------------------------------------
    sub interAtomicDistanceKarKar
    {
        my @v = interAtomicVectorKarKar(@_);
        return sqrt($v[0]*$v[0] + $v[1]*$v[1] + $v[2]*$v[2]);
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

#===================================================================================================
# Subroutine intraAtomicVectorKarKar($data, $a1, $a2).
#---------------------------------------------------------------------------------------------------
#   DESCRIPTION:    Returns a vector in full cartesian coordinates from atom $a1 to atom $a2
#					in @data. This function	takes cartesian coordinates and returns cartesian coordinates. 
#
#   INPUT:          see DESCRIPTION
#
#   OUTPUT:         a 1x3 vector in full cartesian coordinates
#---------------------------------------------------------------------------------------------------
    sub intraAtomicVectorKarKar
    {
        my ($data, $a1, $a2) = @_;
        #my @data = @{dclone($data)};
        my @aCoords = $data->[0][$a1];
        my @bCoords = $data->[0][$a2];
        return ($bCoords[0][0] - $aCoords[0][0], $bCoords[0][1] - $aCoords[0][1], $bCoords[0][2] - $aCoords[0][2]);
    }

#===================================================================================================
# Subroutine intraAtomicDistanceDirKar($data, $a1, $a2).
#---------------------------------------------------------------------------------------------------
#   DESCRIPTION:	Returns the distance between two atoms in full cartesian units.     
#
#   INPUT:          $data: 	the poscar structure in direct coordinates
#					$a1: 	the index of atom1
#					$a2: 	the index of atom2
#
#   OUTPUT:         a scalar value representing the distance in full cartesian units 
#---------------------------------------------------------------------------------------------------
    sub intraAtomicDistanceDirKar
    {
        my ($data, $a1, $a2) = @_;
        #my @data = @{dclone($data)};
        my @aCoords = $data->[0][$a1];
        my @bCoords = $data->[0][$a2];
        my @diff = pbc_difference(\@bCoords, \@aCoords, 1);
        dirkar(@diff, $data->[POSCAR_BASIS], $data->[POSCAR_LATTICE], 1);
        my $x = $diff[0][0][0];
        my $y = $diff[0][0][1];
        my $z = $diff[0][0][2];
        return sqrt($x * $x + $y * $y + $z * $z);
    }

#===================================================================================================
## Subroutine atomicDistanceTableDir($data).
##---------------------------------------------------------------------------------------------------
##   DESCRIPTION:    This routine returns a table with the distances between every atom in the Direct 
##                   Coordinate poscar data.
##
##   INPUT:          $data: the data list from the read_poscar
##
##   OUTPUT:         The distance table: $table[1][13] is the distance between atom 1 and 13.
##                                       $table[13][1] is the same thing.
##---------------------------------------------------------------------------------------------------
#    sub atomicDistanceTableDir
#    {
#        my $data = shift;
#        my @data = @$data;
#        my @table = ();
#        for($a = 0; $a < $data[POSCAR_TOTAL_ATOMS] - 1; $a++)
#        {
#            $table[$a][$a] = 0.00000;
#            for($b = $a + 1; $b < $data[POSCAR_TOTAL_ATOMS]; $b++)
#            {
#                my @aCoords = $data[0][$a];
#                my @bCoords = $data[0][$b];
#                my @diff = pbc_difference(\@aCoords, \@bCoords, 1);
#                dirkar(@diff, $data[POSCAR_BASIS], $data[POSCAR_LATTICE], 1);
#                my $x = $diff[0][0][0];
#                my $y = $diff[0][0][1];
#                my $z = $diff[0][0][2];
#                my $dist = sqrt($x * $x + $y * $y + $z * $z);
#                $table[$a][$b] = $dist;
#                $table[$b][$a] = $dist;
#            }
#        }
#        $table[$data[POSCAR_TOTAL_ATOMS] - 1][$data[POSCAR_TOTAL_ATOMS] - 1] = 0.0000;
#        return @table;
#    }

##===================================================================================================
## Subroutine atomicDistanceTableKar($data).
##---------------------------------------------------------------------------------------------------
##   DESCRIPTION:    This routine returns a table with the distances between every atom in the 
##                   Cartesian Coordinate poscar data.
##
##   INPUT:          $data: the data list from the read_poscar
##
##   OUTPUT:         The distance table: $table[1][13] is the distance between atom 1 and 13.
##                                       $table[13][1] is the same thing.
##---------------------------------------------------------------------------------------------------
#    sub atomicDistanceTableKar
#    {
#        my $data = shift;
#        my @data = @$data;
#        my @table = ();
#        for($a = 0; $a < $data[POSCAR_TOTAL_ATOMS]; $a++)
#        {
#            for($b = $a; $b < $data[POSCAR_TOTAL_ATOMS]; $b++)
#            {
#                my @aCoords = $data[0][$a];
#                my @bCoords = $data[0][$b];
#                my $x = $data[POSCAR_COORDINATES][$a][0] - $data[POSCAR_COORDINATES][$b][0];
#                my $y = $data[POSCAR_COORDINATES][$a][1] - $data[POSCAR_COORDINATES][$b][1];
#                my $z = $data[POSCAR_COORDINATES][$a][2] - $data[POSCAR_COORDINATES][$b][2];
#                my $dist = sqrt($x * $x + $y * $y + $z * $z);
#                $table[$a][$b] = $dist;
#                $table[$b][$a] = $dist;
#            }
#        }
#        return @table;
#    }

#===================================================================================================
# Subroutine hasElement($e, $data).
#---------------------------------------------------------------------------------------------------
#   DESCRIPTION:    This routine tells you if the poscar data has a certain element in it
#
#   INPUT:          $e: the element, by atomic symbol.
#                   $data: the data list from the read_poscar
#
#   OUTPUT:         1 or 0
#---------------------------------------------------------------------------------------------------
    sub hasElement
    {
        my $e = shift;
        my $data = shift;
        my @data = @$data;
        my @elements = split(/ /, $data[7]);
        foreach(@elements)
        {
            if($_ eq $e)
            {
                return 1;
            }
        }
        return 0;
    }
#===================================================================================================
# Subroutine arrayContains($array, $value)
#---------------------------------------------------------------------------------------------------
#   DESCRIPTION:    Utility that tells you whether or not a given array contains a given value.
#
#---------------------------------------------------------------------------------------------------
    sub arrayContains
    {
        my ($array, $value) = @_;
        my $len = scalar @$array;
        for(my $i = 0; $i < $len; $i++)
        {
            if($array->[$i] == $value)
            {
                return 1;
            }
        }
        return 0;
    }
    
    
#===================================================================================================
# Subroutine printCar($data)
#---------------------------------------------------------------------------------------------------
#   DESCRIPTION:    takes a pointer to a poscar structure and returns a string representing the data
#
#---------------------------------------------------------------------------------------------------
    sub carString
    {
        my ($data) = @_;
        my $str = "";
        for(my $i = 0; $i < $data->[POSCAR_TOTAL_ATOMS]; $i++)
        {
			$str .= "$i: $data->[POSCAR_COORDINATES][$i][0], $data->[POSCAR_COORDINATES][$i][1], $data->[POSCAR_COORDINATES][$i][2]\n";
        }
		$str .= "basis-x: $data->[POSCAR_BASIS][0][0], $data->[POSCAR_BASIS][0][1], $data->[POSCAR_BASIS][0][2]\n";
		$str .= "basis-y: $data->[POSCAR_BASIS][1][0], $data->[POSCAR_BASIS][1][1], $data->[POSCAR_BASIS][1][2]\n";
		$str .= "basis-z: $data->[POSCAR_BASIS][2][0], $data->[POSCAR_BASIS][2][1], $data->[POSCAR_BASIS][2][2]\n";
		$str .= "lattice constant: $data->[POSCAR_LATTICE]\n";
		$str .= "total atoms: $data->[POSCAR_TOTAL_ATOMS]\n";
        return $str;
    }
    
    
    
    
    
