#===================================================================================================

package kdbutil;

#===================================================================================================
# Exports
#---------------------------------------------------------------------------------------------------

use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(stdatStatus stdatQuality stdatBarrier stdatMin1 stdatMin2 crossProduct dotProduct);
push @EXPORT, qw(matchDescription carDropAtom atomicDistanceTableKar);
push @EXPORT, qw(%covalentRadii %NON_METALS %INVERSE_NON_METALS appendXYZ appendXYZs);
push @EXPORT, qw(rotatePointXYZ rotatePointAxis rotatePointsAxis CarElementTable elementType atomicDistanceTableDir);
push @EXPORT, qw(getProcessDirs kdbHome DirKar KarDir loadXYZ2CAR writeXYZ carElementCountHash);
push @EXPORT, qw(nthAtomElement centerOfMass carNumComponents carComponentCount nthAtomComponent);
push @EXPORT, qw(@rotationVectors unitVector neighborTypeCount neighborTypeCountDifference);
push @EXPORT, qw(POSCAR_COORDINATES POSCAR_BASIS POSCAR_LATTICE POSCAR_NUM_ATOMS POSCAR_TOTAL_ATOMS);
push @EXPORT, qw(POSCAR_SELECTIVE_FLAG POSCAR_SELECTIVE POSCAR_DESCRIPTION);


#===================================================================================================
# Imports
#---------------------------------------------------------------------------------------------------

use Storable qw(dclone);
#use Carp 'verbose';
#$SIG{ __DIE__ } = sub { Carp::confess( @_ ) };
use Vasp;

#===================================================================================================
# Constants
#---------------------------------------------------------------------------------------------------

# The indices of the read_poscar data format.
use constant POSCAR_COORDINATES      => 0;
use constant POSCAR_BASIS            => 1;
use constant POSCAR_LATTICE          => 2;
use constant POSCAR_NUM_ATOMS        => 3;
use constant POSCAR_TOTAL_ATOMS      => 4;
use constant POSCAR_SELECTIVE_FLAG   => 5;
use constant POSCAR_SELECTIVE        => 6;
use constant POSCAR_DESCRIPTION      => 7;

#===================================================================================================
# Data
#---------------------------------------------------------------------------------------------------

# A list of covalent radii.
our %covalentRadii = (Ac => 1.88, Ag => 1.59, Al => 1.35, Am => 1.51, Ar => 1.51, As => 1.21, At => 1.21, Au => 1.5, B => 0.83, Ba => 1.34, Be => 0.35, Bh => 1.5, Bi => 1.54, Bk => 1.54, Br => 1.21, C => 0.68, Ca => 0.99, Cd => 1.69, Ce => 1.83, Cf => 1.83, Cl => 0.99, Cm => 0.99, Co => 1.33, Cr => 1.35, Cs => 1.67, Cu => 1.52, Db => 1.5, Ds => 1.5, Dy => 1.75, Er => 1.73, Es => 1.5, Eu => 1.99, F => 0.64, Fe => 1.34, Fm => 1.5, Fr => 1.5, Ga => 1.22, Gd => 1.79, Ge => 1.17, H => 0.23, He => 1.5, Hf => 1.57, Hg => 1.7, Ho => 1.74, Hs => 1.5, I => 1.4, In => 1.63, Ir => 1.32, K => 1.33, Kr => 1.5, La => 1.87, Li => 0.68, Lr => 1.5, Lu => 1.72, Md => 1.5, Mg => 1.1, Mn => 1.35, Mo => 1.47, Mt => 1.5, N => 0.68, Na => 0.97, Nb => 1.48, Nd => 1.81, Ne => 1.5, Ni => 1.5, No => 1.5, Np => 1.55, O => 0.68, Os => 1.37, P => 1.05, Pa => 1.61, Pb => 1.54, Pd => 1.5, Pm => 1.8, Po => 1.68, Pr => 1.82, Pt => 1.5, Pu => 1.53, Ra => 1.9, Rb => 1.47, Re => 1.35, Rf => 1.5, Rh => 1.45, Rn => 1.5, Ru => 1.4, S => 1.02, Sb => 1.46, Sc => 1.44, Se => 1.22, Sg => 1.5, Si => 1.2, Sm => 1.8, Sn => 1.46, Sr => 1.12, Ta => 1.43, Tb => 1.76, Tc => 1.35, Te => 1.47, Th => 1.79, Ti => 1.47, Tl => 1.55, Tm => 1.72, U => 1.58, V => 1.33, W => 1.37, Xe => 1.5, Y => 1.78, Yb => 1.94, Zn => 1.45, Zr => 1.56);

 # A list of nonmetals and its inverse. 
our %NON_METALS = ('H' => 1, 'He' => 2, 'B' => 3, 'C' => 4, 'N' => 5, 'O' => 6, 'F' => 7, 'Ne' => 8, 'Si' => 9, 'P' => 10, 'S' => 11, 'Cl' => 12, 'Ar' => 13, 'As' => 14, 'Se' => 15, 'Br' => 16, 'Kr' => 17, 'Te' => 18, 'I' => 19, 'Xe' => 20, 'At' => 21, 'Rn' => 22);
our %INVERSE_NON_METALS = (1 => 'H', 2 => 'He', 3 => 'B', 4 => 'C', 5 => 'N', 6 => 'O', 7 => 'F', 8 => 'Ne', 9 => 'Si', 10 => 'P', 11 => 'S', 12 => 'Cl', 13 => 'Ar', 14 => 'As', 15 => 'Se', 16 => 'Br', 17 => 'Kr', 18 => 'Te', 19 => 'I', 20 => 'Xe', 21 => 'At', 22 => 'Rn');

#===================================================================================================
# Subroutine atomicDistanceTableDir($data).
#---------------------------------------------------------------------------------------------------
#   DESCRIPTION:    This routine returns a table with the distances between every atom in the Direct 
#                   Coordinate poscar data.
#
#   INPUT:          $data: the data list from the read_poscar
#
#   OUTPUT:         The distance table: $table[1][13] is the distance between atom 1 and 13.
#                                       $table[13][1] is the same thing.
#---------------------------------------------------------------------------------------------------
    sub atomicDistanceTableDir
    {
        my $data = shift;
        my @data = @$data;
        my @table = ();
        for($a = 0; $a < $data[POSCAR_TOTAL_ATOMS] - 1; $a++)
        {
            $table[$a][$a] = 0.00000;
            for($b = $a + 1; $b < $data[POSCAR_TOTAL_ATOMS]; $b++)
            {
                my @aCoords = $data[0][$a];
                my @bCoords = $data[0][$b];
                my @diff = pbc_difference(\@aCoords, \@bCoords, 1);
                dirkar(@diff, $data[POSCAR_BASIS], $data[POSCAR_LATTICE], 1);
                my $x = $diff[0][0][0];
                my $y = $diff[0][0][1];
                my $z = $diff[0][0][2];
                my $dist = sqrt($x * $x + $y * $y + $z * $z);
                $table[$a][$b] = $dist;
                $table[$b][$a] = $dist;
            }
        }
        $table[$data[POSCAR_TOTAL_ATOMS] - 1][$data[POSCAR_TOTAL_ATOMS] - 1] = 0.0000;
        return @table;
    }


#===================================================================================================
# Subroutine atomicDistanceTableKar($data).
#---------------------------------------------------------------------------------------------------
#   DESCRIPTION:    This routine returns a table with the distances between every atom in the Direct 
#                   Coordinate poscar data.
#
#   INPUT:          $data: the data list from the read_poscar
#
#   OUTPUT:         The distance table: $table[1][13] is the distance between atom 1 and 13.
#                                       $table[13][1] is the same thing.
#---------------------------------------------------------------------------------------------------
    sub atomicDistanceTableKar
    {
        my $data = shift;
        my @data = @$data;
        my @table = ();
        for($a = 0; $a < $data[POSCAR_TOTAL_ATOMS] - 1; $a++)
        {
            $table[$a][$a] = 0.00000;
            for($b = $a + 1; $b < $data[POSCAR_TOTAL_ATOMS]; $b++)
            {
                my @aCoords = $data[0][$a];
                my @bCoords = $data[0][$b];
                my $x = $data[POSCAR_COORDINATES][$a][0] - $data[POSCAR_COORDINATES][$b][0];
                my $y = $data[POSCAR_COORDINATES][$a][1] - $data[POSCAR_COORDINATES][$b][1];
                my $z = $data[POSCAR_COORDINATES][$a][2] - $data[POSCAR_COORDINATES][$b][2];
                my $dist = sqrt($x * $x + $y * $y + $z * $z);
                $table[$a][$b] = $dist;
                $table[$b][$a] = $dist;
            }
        }
        $table[$data[POSCAR_TOTAL_ATOMS] - 1][$data[POSCAR_TOTAL_ATOMS] - 1] = 0.0000;
        return @table;
    }

#===================================================================================================
# Subroutine neighborTypeCountDifference($data).
#---------------------------------------------------------------------------------------------------
    sub neighborTypeCountDifference
    {
        my ($tc1, $tc2) = @_;
        my $diff = 0;
        if(length(keys(%$tc1)) != length(keys(%$tc2)))
        {
            return 1000;
        }
        for my $t1(keys %$tc1)
        {
            if ($t1 ne 'total')
            {
                if (exists $tc2->{$t1})
                {
                    $diff += abs($tc2->{$t1} - $tc1->{$t1})
                }
                else
                {
                    return 1000;
                }
            }
        }
        return $diff;
    }


#===================================================================================================
# Subroutine neighborTypeCount($data).
#---------------------------------------------------------------------------------------------------
    sub neighborTypeCount
    {
        my ($car, $distances, $fudge) = @_;
        my $typeCount = {};
        my @types = CarElementTable($car);
        for (my $a = 0; $a < $car->[POSCAR_TOTAL_ATOMS]; $a++)
        {
            if (index($car->[POSCAR_SELECTIVE][$a], "F") != -1) { next; }
            $typeCount->{$a} = {};
            $typeCount->{$a}->{'total'} = 0;
            for(my $b = 0; $b < $car->[POSCAR_TOTAL_ATOMS]; $b++)
            {
                if($b == $a) { next; }
                if ($distances->[$a][$b] < ($covalentRadii{$types[$a]} + $covalentRadii{$types[$b]}) * (1.0 + $fudge))
                {
                    if (exists $typeCount->{$a}->{$types[$b]})
                    {
                        $typeCount->{$a}->{$types[$b]}++;
                    }
                    else
                    {
                        $typeCount->{$a}->{$types[$b]} = 1;
                    }
                    $typeCount->{$a}->{'total'}++;
                }
            }
        }
        return $typeCount;
    }

#===================================================================================================
# Subroutine unitVector($data).
#---------------------------------------------------------------------------------------------------
#   DESCRIPTION:    Returns a 3-component unit vector
#
#   INPUT:          $data: the vector to normalize
#
#   OUTPUT:         a unit vector
#---------------------------------------------------------------------------------------------------
    sub unitVector
    {
        my $data = shift;
        my @data = @$data;
        my $x = $data[0];
        my $y = $data[1];
        my $z = $data[2];
        my $mag = sqrt($x * $x + $y * $y + $z * $z);
        if($mag == 0)
        {
            #confess "divide by zero imminent.";
            return ($x, $y, $z);
        }
        $x /= $mag;
        $y /= $mag;
        $z /= $mag;
		return ($x, $y, $z);
    }


#===================================================================================================
# Subroutine rotatePointAxis(x, y, z, xa, ya, za, theta)
#---------------------------------------------------------------------------------------------------
#   Rotate a point (x, y, z) around the (xa, ya, za) axis by theta.
#---------------------------------------------------------------------------------------------------
    sub rotatePointAxis
    {
        my $x = shift;
        my $y = shift;
        my $z = shift;
        my $xa = shift;
        my $ya = shift;
        my $za = shift;
        my $theta = shift;
        my $magnitude = sqrt($xa*$xa + $ya*$ya + $za*$za);
        if($magnitude > 0)
        {
            $xa /= $magnitude;
            $ya /= $magnitude;
            $za /= $magnitude;
            my $t = 1.0 - cos($theta);
            my $c = cos($theta);
            my $s = sin($theta);
            my @R = ();
            my $txaxa = $t * $xa * $xa;
            my $txaya = $t * $xa * $ya;
            my $txaza = $t * $xa * $za;
            my $tyaya = $t * $ya * $ya;
            my $tyaza = $t * $ya * $za;
            my $tzaza = $t * $za * $za;
            $R[0][0] = $txaxa + $c;
            $R[0][1] = $txaya + $s * $za;
            $R[0][2] = $txaza - $s * $ya;
            $R[1][0] = $txaya - $s * $za;
            $R[1][1] = $tyaya + $c;
            $R[1][2] = $tyaza + $s * $xa;
            $R[2][0] = $txaza + $s * $ya;
            $R[2][1] = $tyaza - $s * $xa;
            $R[2][2] = $tzaza + $c;
            my $xp = $R[0][0] * $x + $R[0][1] * $y + $R[0][2] * $z;
            my $yp = $R[1][0] * $x + $R[1][1] * $y + $R[1][2] * $z;
            my $zp = $R[2][0] * $x + $R[2][1] * $y + $R[2][2] * $z;
            return ($xp, $yp, $zp);
        }
        else
        {
            return ($x, $y, $z);
        }
    }

#===================================================================================================
# Subroutine rotatePointsAxis(xyzs, n, xa, ya, za, theta)
#---------------------------------------------------------------------------------------------------
#   Rotate n points (xyz) around the (xa, ya, za) axis by theta.
#---------------------------------------------------------------------------------------------------
    sub rotatePointsAxis
    {
        my ($xyzs, $n, $xa, $ya, $za, $theta) = @_;
#        print "\n$n $xa $ya $za $theta\n";
        my $xyzs = dclone($xyzs);
        my $magnitude = sqrt($xa*$xa + $ya*$ya + $za*$za);
        if($magnitude > 0)
        {
            $xa /= $magnitude;
            $ya /= $magnitude;
            $za /= $magnitude;
            my $t = 1.0 - cos($theta);
            my $c = cos($theta);
            my $s = sin($theta);
            my @R = ();
            my $txaxa = $t * $xa * $xa;
            my $txaya = $t * $xa * $ya;
            my $txaza = $t * $xa * $za;
            my $tyaya = $t * $ya * $ya;
            my $tyaza = $t * $ya * $za;
            my $tzaza = $t * $za * $za;
            $R[0][0] = $txaxa + $c;
            $R[0][1] = $txaya + $s * $za;
            $R[0][2] = $txaza - $s * $ya;
            $R[1][0] = $txaya - $s * $za;
            $R[1][1] = $tyaya + $c;
            $R[1][2] = $tyaza + $s * $xa;
            $R[2][0] = $txaza + $s * $ya;
            $R[2][1] = $tyaza - $s * $xa;
            $R[2][2] = $tzaza + $c;
            for(my $i = 0; $i < $n; $i++)
            {
                my $xp = $R[0][0] * $xyzs->[$i][0] + $R[0][1] * $xyzs->[$i][1] + $R[0][2] * $xyzs->[$i][2];
                my $yp = $R[1][0] * $xyzs->[$i][0] + $R[1][1] * $xyzs->[$i][1] + $R[1][2] * $xyzs->[$i][2];
                my $zp = $R[2][0] * $xyzs->[$i][0] + $R[2][1] * $xyzs->[$i][1] + $R[2][2] * $xyzs->[$i][2];
                $xyzs->[$i][0] = $xp;
                $xyzs->[$i][1] = $yp;
                $xyzs->[$i][2] = $zp;
            }
            return $xyzs;
        }
        else
        {
#            print "\nERROR\n";
            return $xyzs;
        }
    }

#===================================================================================================
# Subroutine rotatePointXYZ(x, y, z, phi, theta, psi)
#---------------------------------------------------------------------------------------------------
#   Rotate a point (x, y, z) around the X, then Y, then Z axes, by phi, theta, and psi, respectively.
#---------------------------------------------------------------------------------------------------
    sub rotatePointXYZ
    {
        my $x = shift;
        my $y = shift;
        my $z = shift;
        my $phi = shift;
        my $theta = shift;
        my $psi = shift;
        # Rotate around the x-axis.
        my $zp = $z * cos($phi) - $y * sin($phi);
        my $yp = $z * sin($phi) + $y * cos($phi);
        $z = $zp;
        $y = $yp;
        # Rotate around the y-axis.
        my $xp = $x * cos($theta) - $z * sin($theta);
        my $zp = $x * sin($theta) + $z * cos($theta);
        $z = $zp;
        $x = $xp;
        # Rotate around the z-axis.
        my $xp = $x * cos($psi) - $y * sin($psi);
        my $yp = $x * sin($psi) + $y * cos($psi);
        $x = $xp;
        $y = $yp;
        # Update the coordinates.
        return ($x, $y, $z);
    }
    
#===================================================================================================
# Subroutine crossProduct(a, b).
#---------------------------------------------------------------------------------------------------
#   Returns a 3D vector representing the cross product between vectors a and b.
#
#---------------------------------------------------------------------------------------------------
    sub crossProduct
    {
        my $a = shift;
        my @a = @$a;
        my $b = shift;
        my @b = @$b;
        return ($a[1] * $b[2] - $a[2] * $b[1], $a[2] * $b[0] - $a[0] * $b[2], $a[0] * $b[1] - $a[1] * $b[0]);
    }

#===================================================================================================
# Subroutine dotProduct(a, b).
#---------------------------------------------------------------------------------------------------
#   Returns the dot product between vectors a and b.
#
#---------------------------------------------------------------------------------------------------
    sub dotProduct
    {
        my $a = shift;
        my @a = @$a;
        my $b = shift;
        my @b = @$b;
        return $a[0] * $b[0] + $a[1] * $b[1] + $a[2] * $b[2];
    }

#===================================================================================================
# Subroutine writeXYZ(data, fileName).
#---------------------------------------------------------------------------------------------------
#   Write out an xyz file of a config already in cartesian coordinates.
#
#---------------------------------------------------------------------------------------------------
    sub writeXYZ
    {
        my $data = shift;
        my @data = @$data;
        my $filename = shift;
        open (xyz, ">$filename");
        print xyz $data[POSCAR_TOTAL_ATOMS] . "\n";
        print xyz $data[POSCAR_DESCRIPTION] . "\n";
        for(my $a = 0; $a < $data[POSCAR_TOTAL_ATOMS]; $a++)
        {
            print xyz nthAtomElement($a, \@data) . " ";
            print xyz $data[POSCAR_COORDINATES][$a][0] . " ";
            print xyz $data[POSCAR_COORDINATES][$a][1] . " ";
            print xyz $data[POSCAR_COORDINATES][$a][2] . "\n";
        }
        close xyz;
    }

#===================================================================================================
# Subroutine appendXYZ(data, fileName).
#---------------------------------------------------------------------------------------------------
#   Append to an xyz file of a config already in cartesian coordinates.
#
#---------------------------------------------------------------------------------------------------
    sub appendXYZ
    {
        my $data = shift;
        my @data = @$data;
        my $filename = shift;
        open (xyz, ">>$filename");
        print xyz $data[POSCAR_TOTAL_ATOMS] . "\n";
        print xyz $data[POSCAR_DESCRIPTION] . "\n";
        for(my $a = 0; $a < $data[POSCAR_TOTAL_ATOMS]; $a++)
        {
            print xyz nthAtomElement($a, \@data) . " ";
            print xyz $data[POSCAR_COORDINATES][$a][0] . " ";
            print xyz $data[POSCAR_COORDINATES][$a][1] . " ";
            print xyz $data[POSCAR_COORDINATES][$a][2] . "\n";
        }
        close xyz;
    }

#===================================================================================================
# Subroutine appendXYZs(data, fileName, comment).
#---------------------------------------------------------------------------------------------------
#   Append to an xyz file of a config already in cartesian coordinates.
#
#---------------------------------------------------------------------------------------------------
    sub appendXYZs
    {
        ($data, $filename, $comment) = @_;
        @data = @$data;
        open (xyz, ">>$filename");
        $totalatoms = 0;
        for(my $i = 0; $i < @data; $i++)
        {
            $totalatoms += $data[i][POSCAR_TOTAL_ATOMS];
        }
        print xyz "$totalatoms\n";
        print xyz "$comment\n";
        for $dr (@data)
        {
            @d = @$dr;
            for(my $a = 0; $a < $d[POSCAR_TOTAL_ATOMS]; $a++)
            {
                print xyz nthAtomElement($a, \@d) . " ";
                print xyz $d[POSCAR_COORDINATES][$a][0] . " ";
                print xyz $d[POSCAR_COORDINATES][$a][1] . " ";
                print xyz $d[POSCAR_COORDINATES][$a][2] . "\n";
            }
        }
        close xyz;
    }

#===================================================================================================
# Subroutine carElementCountHash(@data)                                                             
#---------------------------------------------------------------------------------------------------
#   Returns a hash that contains each element of the poscar pointing to the number of atoms of that
#   element.
#
#---------------------------------------------------------------------------------------------------
    sub carElementCountHash
    {
        my $data = shift;
        my @data = @$data;
        my %hash = ();
        for(my $i = 0; $i < $data[POSCAR_TOTAL_ATOMS]; $i++)
        {
            my $e = nthAtomElement($i, \@data);
            if($hash{$e})
            {
                $hash{$e} += 1;
            }
            else
            {
                $hash{$e} = 1;
            }
        }
        return %hash;
    }

#===================================================================================================
# Subroutine carDropAtom(selected, data).
#---------------------------------------------------------------------------------------------------
#   Remove selected atom from a POSCAR.
#
#---------------------------------------------------------------------------------------------------
    sub carDropAtom
    {
        my $selected = shift;
        my $data = dclone(shift);
        my @data = @$data;
        my $component = nthAtomComponent($selected, \@data);
        my $element = nthAtomElement($selected, \@data);
        my @desc = carDescArray(\@data);
        # Fix the coords.
        splice(@{@data[POSCAR_COORDINATES]}, $selected, 1);
        # Fix the frozen atoms.
        splice(@{@data[POSCAR_SELECTIVE]}, $selected, 1);
        # Fix num_atoms.
        if(carComponentCount($component, \@data) <= 1)
        {
            # Remove the component.
            splice(@{$data[POSCAR_NUM_ATOMS]}, $component, 1);
            
            # Fix the desc.
            splice(@desc, $component, 1);
            $data[POSCAR_DESCRIPTION] = join(" ", @desc);
        }
        else
        {
            # Reduce the component count.
            $data[POSCAR_NUM_ATOMS][$component] -= 1;
        }
        # Reduce the total atom count.
        $data[POSCAR_TOTAL_ATOMS] -= 1;
        return @data;        
    }

#===================================================================================================
# Subroutine centerOfMass(data).
#---------------------------------------------------------------------------------------------------
#   Get the center of mass of all the atoms.
#
#---------------------------------------------------------------------------------------------------
    sub centerOfMass
    {
        my $data = shift;
        my @data = @$data;
        my $xTotal = 0;
        my $yTotal = 0;
        my $zTotal = 0;
        for(my $a = 0; $a < $data[POSCAR_TOTAL_ATOMS]; $a++)
        {
            $xTotal += $data[POSCAR_COORDINATES][$a][0];
            $yTotal += $data[POSCAR_COORDINATES][$a][1];
            $zTotal += $data[POSCAR_COORDINATES][$a][2];
        }
        my $x = $xTotal / $data[POSCAR_TOTAL_ATOMS];
        my $y = $yTotal / $data[POSCAR_TOTAL_ATOMS];
        my $z = $zTotal / $data[POSCAR_TOTAL_ATOMS];
        return ($x, $y, $z);
    }


#===================================================================================================
# Subroutine loadXYZ2CAR(fileName)
#---------------------------------------------------------------------------------------------------
#       Returns an XYZ loaded into POSCAR format.
#---------------------------------------------------------------------------------------------------
    sub loadXYZ2CAR
    {
        # The xyz file filename.
        my $fileName = shift;
        
        # Open the file and get the number of atoms and comment.
        open (XYZ, $fileName);
        my $numAtoms = <XYZ>; chomp $numAtoms;
        my $comment = <XYZ>; chomp $comment;
        
        # Create a square (cartesian) basis.
        my @basis = ([1, 0, 0], [0, 1, 0], [0, 0, 1]);

        # Get the atom type and coordinates.
        my @types = ();
        my @coordinates = ();
        for(my $i = 0; $i < $numAtoms; $i++)
        {
            my $line = <XYZ>;
            chomp $line;
            my @data = split(" ", $line);
            $types[$i] = $data[0];
            $coordinates[$i] = [$data[1], $data[2], $data[3]];
        }
    
        # Create a hash of atom types and the number of atoms of each.
        my %atomTypes = ();
        my $numAtomTypes = 0;
        for(my $i = 0; $i < $numAtoms; $i++)
        {
            if(!$atomTypes{$types[$i]})
            {
                $atomTypes{$types[$i]} = 1;
                $numAtomTypes++;
            }
            else
            {
                $atomTypes{$types[$i]}++;
            }
        }

        # Create an array containing unique descriptions, i.e., ['Cu', 'C', 'O', 'H'].
        my @descArray = ();
        foreach my $key(keys %atomTypes)
        {
            push @descArray, $key;
        }
        
        # Sort @descArray and store it in $desc.
        @descArray = sort @descArray;
        my $desc = join(" ", @descArray);
        
        # Create an array of the number of atoms of each element sorted by @descArray.
        my @numAtomsEachType = ();
        foreach my $type(@descArray)
        {
            push @numAtomsEachType, $atomTypes{$type};
        }
        
        # Sort the coordinates by $descArray and create a hash that points from the old atom number
        # to the new atom number.
        my @tempCoords = ();
        my $tempIndex = 0;
        my %switch = ();
        foreach my $type(@descArray)
        {
            for(my $i = 0; $i < $numAtoms; $i++)
            {
                if($types[$i] eq $type)
                {
                    push @tempCoords, $coordinates[$i];
                    $switch{$i} = $tempIndex;
                    $tempIndex++;
                }
            }
        }
        @coordinates = @tempCoords;
        
        # Fill selective array with zeros.
        my @selective = ();
        for(my $i = 0; $i < $numAtoms; $i++)
        {
            $selective[$i] = 0;
        }
            
        # Return the POSCAR array.
        return (\@coordinates, \@basis, 1.0, \@numAtomsEachType, $numAtoms, 0, \@selective, $desc, \%switch);
    }        

#===================================================================================================
# Subroutine kdbHome()
#---------------------------------------------------------------------------------------------------
#       Returns the directory. The environment variable VTST_KDB holds this if it is defined, 
#       otherwise the default of "userhome/kdb" is retured.
#---------------------------------------------------------------------------------------------------
    sub kdbHome
    {
        # If the VTST_KDB environment variable is not set, return the default (userhome/kdb).
        if (!$ENV{'VTST_KDB'}) 
        {
            my $userdir = $ENV{'HOME'};
            if(! -d "$userdir/kdb")
            {
                mkdir("$userdir/kdb", 0777);
            }
            return "$userdir/kdb";
        }
        # Otherwise, return the value of the VTST_KDB environment variable.
        else
        {
            return $ENV{'VTST_KDB'};
        }
    }
    
#===================================================================================================
# Subroutine getProcessDirs()
#---------------------------------------------------------------------------------------------------
#   DESCRIPTION:    Get the list of kdb process subdirectories that match the global $Desc 
#                   description of the process we're trying to match.
#   INPUT:          The kdb directory (from kdbHome).
#                   The description string we are trying to match, i.e., "C Cu H O".
#   OUTPUT:         All matching directories.
#---------------------------------------------------------------------------------------------------
    sub getProcessDirs
    {
        my $kdbDir = shift;
        my $desc = shift;
        my @procDirs;
        my @allSystemDirs = glob("$kdbDir/*");
        foreach my $dir (@allSystemDirs) 
        {
            if(-d $dir)
            {
                if (matchDescription($desc, dir2desc($dir)))
                {
                    my $i = 0;
                    while (-e "$dir/$i")
                    {
                        push(@procDirs, "$dir/$i");
                        $i++;
                    }
                }
            }
        }
        return @procDirs;
    }

#===================================================================================================
# Subroutine dir2desc()
#---------------------------------------------------------------------------------------------------
#   DESCRIPTION:    Converts kdb subdirectory into a description, e.g. "/kdb/C_Cu_H_O" returns 
#                   "C Cu H O"
#---------------------------------------------------------------------------------------------------
    sub dir2desc
    {
        my $desc = shift;
        $desc = substr($desc, rindex($desc, "/"));
        $desc =~ s/_/ /g;
        $desc =~ s#/##;
        return $desc;
    }

#===================================================================================================
# Subroutine elementType (desc)
#---------------------------------------------------------------------------------------------------
#   DESCRIPTION:    Returns the atom type for nonmetals, and "Metal" for metals.
#
#---------------------------------------------------------------------------------------------------
    sub elementType
    {
        my $desc = shift;
        if($NON_METALS{$desc})
        {
            return $desc;
        }
        else
        {
            return "Metal";
        }
    }


#===================================================================================================
# Subroutine matchDescription (desc1, desc2)
#---------------------------------------------------------------------------------------------------
#   DESCRIPTION:    Determines whether or not two kdb descriptions match according to 
#                   THE RULES.
#
#   INPUT:          Two strings representing kdb descriptions, i.e., ("H O C Ag", "H O C Cu")
#
#   OUTPUT:         A 1 or 0 representing true or false for a positive or negative match.
#
#   THE RULES:      The number of elements must be the same.
#                   The nonmetal elements must match in each description.
#---------------------------------------------------------------------------------------------------
    sub matchDescription
    {
        # Get the descriptions.
        my $desc1 = shift;
        my $desc2 = shift;
        # Get the individual elements.
        my @atoms1 = split(/ /, $desc1);
        my @atoms2 = split(/ /, $desc2);
        # Get the number of elements in each description.
        my $len1 = @atoms1;
        my $len2 = @atoms2;
        # If the number of elements does not match, return false (no match).
        if($len1 != $len2)
        {
            return 0;
        }
        # Store the elements inside a hash so that we can perform an array.contains() operation. 
        my %isAtoms1;
        my %isAtoms2;
        foreach(@atoms1)
        {
            $isAtoms1{$_} = 1;
        }
        foreach(@atoms2)
        {
            $isAtoms2{$_} = 1;
        }
        # Check to see if each nonmetal in desc1 matches a nonmetal in desc2.
        foreach(@atoms1)
        {
            if($NON_METALS{$_})
            {
                if(!$isAtoms2{$_})
                {
                    return 0;
                }
            }
        }
        # Check to see if each nonmetal in desc2 matches a nonmetal in desc1.
        # This is repeated in case there are more NonMetals in desc2 than in desc1.
        foreach(@atoms2)
        {
            if($NON_METALS{$_})
            {
                if(!$isAtoms1{$_})
                {
                    return 0;
                }
            }
        }
        # Return true (successful match).
        return 1;
    }

#===================================================================================================
# Subroutine CarElementTable(@data).
#---------------------------------------------------------------------------------------------------
#   DESCRIPTION:    This routine returns an array where the nth element is the element of atom n.
#
#   INPUT:          @data: the data list from the read_poscar
#
#   OUTPUT:         I already told you. 
#---------------------------------------------------------------------------------------------------
    sub CarElementTable
    {
        my $data = shift;
        my @data = @$data;
        my @table = ();
        for(my $i = 0; $i < $data[4]; $i++)
        {
            push(@table, nthAtomElement($i, \@data));
        }
        return @table;
    }

#===================================================================================================
# Subroutine nthAtomElement($a, $data).
#---------------------------------------------------------------------------------------------------
#   DESCRIPTION:    This routine returns the nth atom element name from a read_poscar data 
#                   structure.
#
#   INPUT:          $a: the integer representing the atom desired
#                   $data: the data list from the read_poscar
#
#   OUTPUT:         $type: A string representing the symbol of this atom's element, e.g., "H", 
#---------------------------------------------------------------------------------------------------
    sub nthAtomElement
    {
        my $a = shift;
        my $data = shift;
        my @desc = carDescArray($data);
        return $desc[nthAtomComponent($a, $data)];
    }

#===================================================================================================
# Subroutine carDescArray(\@data).
#---------------------------------------------------------------------------------------------------
#   DESCRIPTION:    This routine returns the POSCAR desc in an array of strings
#
#   INPUT:          $data: the data list from the read_poscar
#
#   OUTPUT:         @desc: the desc in an array, e.g., ["Mg", "O", "Cu"];
#---------------------------------------------------------------------------------------------------
    sub carDescArray
    {
        my $data = shift;
        my @data = @$data;
        my $description = $data[POSCAR_DESCRIPTION];
        my @types = split(" ", $description);
        return @types;
    }

#===================================================================================================
# Subroutine carComponentCount($component, \@data).
#---------------------------------------------------------------------------------------------------
#   DESCRIPTION:    This routine returns the number of atoms of component $component.
#
#   INPUT:          $component: the component we want the count of.
#                   $data: the data list from the read_poscar
#
#   OUTPUT:         The number of atoms in component $component.kdbaddpr.pl REACTANT PRODUCT SADDLE MODE --debug=movie
#---------------------------------------------------------------------------------------------------
    sub carComponentCount
    {
        my $component = shift;
        my $data = shift;
        my @data = @$data;
        return $data[POSCAR_NUM_ATOMS][$component];
    }

#===================================================================================================
# Subroutine carNumComponents(\@data).
#---------------------------------------------------------------------------------------------------
#   DESCRIPTION:    This routine returns the number of components.
#
#   INPUT:          $data: the data list from the read_poscar
#
#   OUTPUT:         the number of components
#---------------------------------------------------------------------------------------------------
    sub carNumComponents
    {
        my $data = shift;
        my @data = @$data;
        return int(scalar(@{$data[POSCAR_NUM_ATOMS]}));
    }

#===================================================================================================
# Subroutine nthAtomComponent($a, $data).
#---------------------------------------------------------------------------------------------------
#   DESCRIPTION:    This routine returns the nth atom component number from a read_poscar data 
#                   structure.
#
#   INPUT:          $a: the integer representing the atom desired
#                   $data: the data list from the read_poscar
#
#   OUTPUT:         $component: the component number of this atom.
#---------------------------------------------------------------------------------------------------
    sub nthAtomComponent
    {
        my $a = shift;
        my $data = shift;
        my @data = @$data;
        my @num_atoms = @{$data[3]};
        my $len = carNumComponents(\@data);
        my $atomTotal = 0;
        my $component = 0;
        for(my $component = 0; $component < $len; $component++)
        {
            $atomTotal += $num_atoms[$component];
            if($a < $atomTotal)
            {
                return $component;
            }
        }
        die "\nFailed in nthAtomComponent.";
    }

#=========================================================================================================
# Subroutine DirKar($data).
#---------------------------------------------------------------------------------------------------------
#   DESCRIPTION:    This converts a poscar data structure to cartesian coordinates from direct 
#                   coordinates.  Does not alter incoming data.
#
#   INPUT:          $data: the data list from the read_poscar
#                   
#   OUTPUT:         The poscar data, now in cartesian coords.
#---------------------------------------------------------------------------------------------------------
    sub DirKar
    {
        my $data = dclone(shift);
        my @data = @$data;
        dirkar($data[0], $data[POSCAR_BASIS], $data[POSCAR_LATTICE], $data[POSCAR_TOTAL_ATOMS]);
        return @data;
    }

#=========================================================================================================
# Subroutine KarDir($data).
#---------------------------------------------------------------------------------------------------------
#   DESCRIPTION:    This converts a poscar data structure to direct coordinates from cartesian 
#                   coordinates.  Does not alter incoming data.
#
#   INPUT:          $data: the data list from the read_poscar
#                   
#   OUTPUT:         The poscar data, now in direct coords.
#---------------------------------------------------------------------------------------------------------
    sub KarDir
    {
        my $data = dclone(shift);
        my @data = @$data;
        kardir($data[0], $data[POSCAR_BASIS], $data[POSCAR_LATTICE], $data[POSCAR_TOTAL_ATOMS]);
        return @data;
    }

#=========================================================================================================
# Subroutine KarDirRef($data, $reference).
#---------------------------------------------------------------------------------------------------------
#   DESCRIPTION:    This converts a poscar data structure to direct coordinates from cartesian 
#                   coordinates.  Does not alter incoming data.
#
#   INPUT:          $data: the data list from the read_poscar
#                   $reference: the poscar structure that contains the BASIS and LATTICE information
#                   
#   OUTPUT:         The poscar data, now in direct coords.
#---------------------------------------------------------------------------------------------------------
    sub KarDirRef
    {
        my ($data, $reference) = @_;
        $data = dclone($data);
        $reference = dclone($reference);
        my @data = @$data;
        my @reference = @$reference;
        kardir($data[0], $reference[POSCAR_BASIS], $reference[POSCAR_LATTICE], $data[POSCAR_TOTAL_ATOMS]);
        return @data;
    }

#===================================================================================================
# Deal with akmc st.dat

sub stdatStatus
{
    my $fileName = shift;
    my $process = shift;
    eval
    {
        open(file, $fileName);
    };
    if($@)
    {
        return -1;
    }
    while(<file>)
    {
        my @values = split(" ", $_);
        if($values[0] eq $process)
        {
            return $values[1];
        }
    }
    return -1;
}

#---------------------------------------------------------------------------------------------------
sub stdatQuality
{
    my $fileName = shift;
    my $process = shift;
    eval
    {
        open(file, $fileName);
    };
    if($@)
    {
        return -1;
    }
    while(<file>)
    {
        my @values = split(' ', $_);
        if($values[0] eq $process)
        {
            return $values[2];
        }
    }
    return -1;
}

#---------------------------------------------------------------------------------------------------
sub stdatBarrier
{
    my $fileName = shift;
    my $process = shift;
    eval
    {
        open(file, $fileName);
    };
    if($@)
    {
        return;
    }
    while(<file>)
    {
        my @values = split(' ', $_);
        if($values[0] eq $process)
        {
            return $values[3];
        }
    }
    return -1;
}

#---------------------------------------------------------------------------------------------------
sub stdatMin1
{
    my $fileName = shift;
    my $process = shift;
    eval
    {
        open(file, $fileName);
    };
    if($@)
    {
        return;
    }
    while(<file>)
    {
        my @values = split(' ', $_);
        if($values[0] eq $process)
        {
            return $values[4];
        }
    }
    return -1;
}

#---------------------------------------------------------------------------------------------------
sub stdatMin2
{
    my $fileName = shift;
    my $process = shift;
    eval
    {
        open(file, $fileName);
    };
    if($@)
    {
        return -1;
    }
    while(<file>)
    {
        my @values = split(' ', $_);
        if($values[0] eq $process)
        {
            return $values[5];
        }
    }
    return -1;
}


#===================================================================================================
# Perl modules must return 1.
    1;
    
    
