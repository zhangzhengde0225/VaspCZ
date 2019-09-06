#!/usr/bin/env perl
#;-*- Perl -*-

# Created by dano 4-16-09
# LAST MODIFIED by dano: 4-29-10
# Reads chgcar files  returns value at r_ions
# Reads locpot files (electrostatic potential) returns value at r_ions
######################################################
# Still need to add periodic boundy cond
######################################################

use FindBin qw($Bin);
use lib "$Bin";
use Vasp;
use Math::Trig;
$fact = 180/pi;

print "\n";
@args=@ARGV;
(@args>=1 || @args>=2) || die "usage: chgvalue.pl <LOCPOT or CHGCAR file> <method>\n";
$inputfilename = $args[0];
if(@args == 2){
    $method = @args[1];
}else{
    $method = 1;
}

$unittype = "angstroms and eV";
$l_units = 1.0;
$e_units = 1;
$unit_flag = 1; # Ang

if ($method==1){$method_flag= "linear"};
if ($method==2){$method_flag= "quadratic"};
if ($method==3){$method_flag= "cubic"};
if ($method==4){$method_flag= "minvalue"};

print "  USING '$unittype' FOR OUTPUT UNITS interpolation medthod ='$method_flag'\n";

$inputfile = "";
open (IN,$inputfilename);
while (<IN>) {$_=~s/^\s+//g;$inputfile.=$_;}
close (IN);

@inputfile = split(/\n/,$inputfile);

# check for a comment line
$line = $inputfile[0] ;
$line =~ s/^\s+//;
@line = split(/\s+/,$line);
if ($line[0] =~ /^\d+(\.\d+)?/) {
    $index = 0;
} else {
    chop($description = $inputfile[0]);
    $index = 1;
}
$lattice = $inputfile[$index];
#print "lattice constant "; print "$lattice\n";

#Check for vasp5 format
$line = $inputfile[$index+4] ;
$line =~ s/^\s+//;
@line = split(/\s+/,$line);
#  if ($line=~/^s/i) {
if ($line[0] =~ /^\d+$/) {
  $vasp5 = 0;
} else {
  chop($description = $inputfile[$index+4]);
  $vasp5 = 1;
}
#print "vasp5_flag ";print "$vasp5\n";

# read total atoms
$num_atoms_ = $inputfile[$index+4+$vasp5];
$num_atoms_ =~ s/^\s+//;
@num_atoms = split(/\s+/,$num_atoms_);
for ($i=0;$i<@num_atoms;$i++) {
    $num_atoms->[$i] = $num_atoms[$i];
    $total_atoms += $num_atoms[$i];
}

for ($i=0; $i<3; $i++) {
    $line = $inputfile[$i+$index+1];
    $line =~ s/^\s+//;
    @line = split(/\s+/,$line);
    # This is how Vasp reads in the basis
    for ($j=0; $j<3; $j++) {
        $basis->[$j][$i] = $line[$j]*$lattice*$l_units;
    }
}

$magv1 = sqrt($basis->[0][0]**2+$basis->[1][0]**2+$basis->[2][0]**2);
$magv2 = sqrt($basis->[0][1]**2+$basis->[1][1]**2+$basis->[2][1]**2);
$magv3 = sqrt($basis->[0][2]**2+$basis->[1][2]**2+$basis->[2][2]**2);

# this is wrong (except for orthogonal box)
#$volume_ang=$magv1*$magv2*$magv3;
# or expansion by minors! 

# correct volume for non-orthognal box
$vol = $basis->[0][0]*($basis->[1][1]*$basis->[2][2] - $basis->[2][1]*$basis->[1][2])
     - $basis->[1][0]*($basis->[0][1]*$basis->[2][2] - $basis->[2][1]*$basis->[0][2])
     + $basis->[2][0]*($basis->[0][1]*$basis->[1][2] - $basis->[1][1]*$basis->[0][3]);

# find inverse basis
$inverse = &inverse($basis); # call function inverse from Vasp
#print "$inverse->[0][0] $inverse->[0][2] $inverse->[0][1]\n";

# make sure volume is +
# this volume is in units selected 
$vol = abs($vol);
#print "volume "; print "$vol\n";

# read in coordinates
$index = $index + $vasp5 + 6;
#print "$index\n";
for ($i=$index; $i<$index+$total_atoms; $i++) {
    $imputfile[$i] =~ s/^\s+//;
    @line = split(/\s+/,$inputfile[$i]);
    for ($j=0; $j<3; $j++) {
        $coordinates->[$i-$index][$j] = $line[$j]; }
    }
#print "$coordinates->[0][0] $coordinates->[0][1] $coordinates->[0][2]\n";
#print "$basis->[0][0] $basis->[0][1] $basis->[0][2]\n";

# change coordinates to cartesian
for ($i=0; $i<$total_atoms; $i++) {
  $v1 = $coordinates->[$i][0]*$basis->[0][0]+$coordinates->[$i][1]*$basis->[0][1]+$coordinates->[$i][2]*$basis->[0][2];
  $v2 = $coordinates->[$i][0]*$basis->[1][0]+$coordinates->[$i][1]*$basis->[1][1]+$coordinates->[$i][2]*$basis->[1][2];
  $v3 = $coordinates->[$i][0]*$basis->[2][0]+$coordinates->[$i][1]*$basis->[2][1]+$coordinates->[$i][2]*$basis->[2][2];
  $coordinates->[$i][0] = $v1;
  $coordinates->[$i][1] = $v2;
  $coordinates->[$i][2] = $v3;
  }

#print "$coordinates->[0][0] $coordinates->[0][1] $coordinates->[0][2]\n";

# read number of grid-points
$index = $index + $total_atoms;
$fft_grid_ = $inputfile[$index];
$fft_grid_ =~ s/^\s+//;
@fft_grid = split(/\s+/,$fft_grid_);
$grid_points = 1;
for ($i=0; $i<@fft_grid; $i++) {
    $fft_grid->[$i] = $fft_grid[$i];
    $grid_points *= $fft_grid[$i];
}

#print "$fft_grid[1]\n";
#print "$grid_points\n";
$data_lines = int($grid_points/5+1);

####READ IN THE RHO DATA####
$index += 1;
$ii = 0;
#print $index;
for ($i=$index; $i<$index+$data_lines; $i++) {
#for ($i=$index;$i<$index+3; $i++) {
    $imputfile[$i] =~ s/^\s+//;
    @line = split(/\s+/,$inputfile[$i]);
    for ($j=0; $j<@line; $j++) {
        # density is read in x-inner,y-middle,z-outer loop
        #$density->[($i-$index)*@line+$j]=$line[$j];
        $density->[$ii] = $line[$j];
        #print "$density->[$ii]\n";
        $ii++;
    }
}
 
#print "$density->[0]\n";
#print "$density->[1]\n";
#print "$density->[2]\n";
#print "$density->[3]\n";

# del inputfile
undef @inputfile;
undef $inputfile;


# it would be nice if grid data was in t x,y,z format
$ii = 0;
for ($i=0; $i<$fft_grid->[2]; $i++) {
    for ($j=0; $j<$fft_grid->[1]; $j++) {
        for ($k=0; $k<$fft_grid->[0]; $k++) {
            $rho->[$k][$j][$i] = $density->[$ii];
            #  print "$density->[$ii]\n";
            $ii++;
        }
    }
}


# find value of chg-file at posions of ions
for ($i=0; $i<$total_atoms; $i++) {
    # change xyz to grid value (function)
    $ijk=&xyz2grid($coordinates->[$i][0],$coordinates->[$i][1],$coordinates->[$i][2]);
    #print "$ijk->[0] $ijk->[1] $ijk->[2]\n";
    #interpolate for density value
    if ($method==1){
        $rho_R->[$i] = &lineinterp($ijk->[0],$ijk->[1],$ijk->[2]);}
    if ($method==2){
        $rho_R->[$i] = &quadinterp($ijk->[0],$ijk->[1],$ijk->[2]);}
    if ($method==3){
        $rho_R->[$i] = &cubicinterp($ijk->[0],$ijk->[1],$ijk->[2]);}
    if ($method==4){
        $rho_R->[$i] = &minvalue($ijk->[0],$ijk->[1],$ijk->[2]);}
    #print "$rho_R->[$i]\n";
}

########## output values at atom posions ########################
if(@args>=3){$outputfilename=$args[2];}
#else{$outputfilename=$inputfilename.".chempot";}
else{$outputfilename="out.chempot";}
open (OUT,">$outputfilename");

#print OUT "Cartesian Coordinates (Bohr)\n";
#print OUT "   Zval       x            y            z        u(x,y,z) hartree\n";
# print data
#print "\n";
#print "Volume DATA R_ions (hartree)\n";
#for ($i=0;$i<$total_atoms;$i++) {
#  printf OUT "%5i %12.6f %12.6f %12.6f %12.6f\n",$atom_type->[$i],$coordinates->[$i][0],$coordinates->[$i][1],$coordinates->[$i][2],$rho_R->[$i];
#  print "Volume DATA R_ions (hartree)\n"
#  print "  $rho_R->[$i]\n";
#}

print OUT "\n";
print OUT "\n";
print OUT "Cartesian Coordinates (Angstroms)\n";
print OUT "   Zval       x            y            z        u(x,y,z) eV\n";
# print data
print "\n";
print "Volume DATA R_ions (eV)\n";

for ($i=0; $i<$total_atoms; $i++) {
    $rho_R_eV->[$i] = $rho_R->[$i];
    printf OUT "%5i %12.6f %12.6f %12.6f %12.6f\n",$atom_type->[$i],$coordinates->[$i][0],$coordinates->[$i][1],$coordinates->[$i][2],$rho_R_eV->[$i];
    print "  $rho_R_eV->[$i]\n";
}
print OUT "\n";
close(OUT);
print "DONE!\n";


#############################################################################
# subroutine to find xyz coordinates of point on the grid
#############################################################################
# assign cartensin coordiantes to grid points (call it center of the volume)
sub grid2xyz {
    my($i,$j,$k) = @_;
    my($xyz);
    my($gx) = $fft_grid[0];
    my($gy) = $fft_grid[1];
    my($gz) = $fft_grid[2];
    $xyz->[0] = ($i/$gx)*$basis->[0][0] + ($j/$gy)*$basis->[0][1] + ($k/$gz)*$basis->[0][2];
    $xyz->[1] = ($i/$gx)*$basis->[1][0] + ($j/$gy)*$basis->[1][1] + ($k/$gz)*$basis->[1][2];
    $xyz->[2] = ($i/$gx)*$basis->[2][0] + ($j/$gy)*$basis->[2][1] + ($k/$gz)*$basis->[2][2];
    return($xyz);
}
##############################################################################

#############################################################################
# subroutine to exact grid point for xyz input (does not return int)
#############################################################################
# (call it center of the volume of grid)
sub xyz2grid {
    my($x,$y,$z) = @_;
    my($ijk);
    my($gx) = $fft_grid[0];
    my($gy) = $fft_grid[1];
    my($gz) = $fft_grid[2];
    $ijk->[0] = $gx*$x*$inverse->[0][0] + $gy*$y*$inverse->[1][0] + $gz*$z*$inverse->[2][0];
    $ijk->[1] = $gx*$x*$inverse->[0][1] + $gy*$y*$inverse->[1][1] + $gz*$z*$inverse->[2][1];
    $ijk->[2] = $gx*$x*$inverse->[0][2] + $gy*$y*$inverse->[1][2] + $gz*$z*$inverse->[2][2];
    return($ijk);
}
##############################################################################

#############################################################################
# linear interpolation for finding density value on grid
#############################################################################
sub lineinterp {
    my($i,$j,$k) = @_;
    # gridpoint prev and next
    my($ip);
    my($in);
    my($jp);
    my($jn);
    my($kp);
    my($kn);
    # distance from prev gridpoint
    my($di);
    my($dj);
    my($dk);
    # interpolation matrix
    my($interp);
    my($val);
    # find gridpoints around point
    # check to see if outside grid bounries
    if ($i>0 and $i+1<$fft_grid->[0]) {
        $ip = int($i);
        $in = int($i+1);
        $di = $i - $ip;
    } else {
        $ip = $fft_grid->[0];
        $in = 0;
        if ($i>0) {
            $di=1+$i;
        } else {
            $di=$i-$ip;
        }
    }
    if ($j>0 and $j+1<$fft_grid->[1]) {
        $jp = int($j);
        $jn = int($j+1);
        $dj = $j-$jp;
    } else {
        $jp = $fft_grid->[1];
        $jn = 0;
        if ($j>0) {
            $dj = 1 + $j;
        } else {
            $dj = $j - $jp;
        }
    }
    if ($k>0 and $k+1<$fft_grid->[2]) {
        $kp = int($k);
        $kn = int($k+1);
        $dk = $k - $kp;}
    else {
        $kp = $fft_grid->[2];
        $kn = 0;
        if ($k>0) {
            $dk = 1 + $k;
        } else {
            $dk = $k - $kp;
        }
    }

    # get 8 values around the real
    $interp->[0][0][0] = $rho->[$ip][$jp][$kp];
    $interp->[2][0][0] = $rho->[$in][$jp][$kp];
    $interp->[0][2][0] = $rho->[$ip][$jn][$kp];
    $interp->[2][2][0] = $rho->[$in][$jn][$kp];
    $interp->[0][0][2] = $rho->[$ip][$jp][$kn];
    $interp->[2][0][2] = $rho->[$in][$jp][$kn];
    $interp->[0][2][2] = $rho->[$ip][$jn][$kn];
    $interp->[2][2][2] = $rho->[$in][$jn][$kn];

    ##### interpolate 12 edges #######
    #$interp->[1][0][0] = (1-$di)*$interp->[0][0][0] + $di*$interp->[2][0][0]; 
    $interp->[0][1][0] = (1-$dj)*$interp->[0][0][0] + $dj*$interp->[0][2][0]; 
    #$interp->[0][0][1] = (1-$dk)*$interp->[0][0][0] + $dk*$interp->[0][0][2];

    $interp->[2][1][0] = (1-$dj)*$interp->[2][0][0] + $dj*$interp->[2][2][0]; 
    #$interp->[2][0][1] = (1-$dk)*$interp->[2][0][0] + $dk*$interp->[2][0][2]; 

    #$interp->[1][2][0] = (1-$di)*$interp->[0][2][0] + $di*$interp->[2][2][0]; 
    #$interp->[0][2][1] = (1-$dk)*$interp->[0][2][0] + $dk*$interp->[0][2][2]; 

  #$interp->[1][0][2] = (1-$di)*$interp->[0][0][2] + $di*$interp->[2][0][2]; 
  $interp->[0][1][2] = (1-$dj)*$interp->[0][0][2] + $dj*$interp->[0][2][2]; 
  
  #$interp->[1][2][2] = (1-$di)*$interp->[0][2][2] + $di*$interp->[2][2][2]; 
  $interp->[2][1][2] = (1-$dj)*$interp->[2][0][2] + $dj*$interp->[2][2][2]; 
  #$interp->[2][2][1] = (1-$dk)*$interp->[2][2][0] + $dk*$interp->[2][2][2];

  #### interpolate 6 faces ########
  #$interp->[1][1][0] = (1-$dj)*$interp->[1][0][0] + $dj*$interp->[1][2][0];
  #$interp->[1][0][1] = (1-$dk)*$interp->[1][0][0] + $dk*$interp->[1][0][2];
  $interp->[0][1][1] = (1-$dk)*$interp->[0][1][0] + $dk*$interp->[0][1][2];

  #$interp->[1][1][2] = (1-$dj)*$interp->[1][0][2] + $dj*$interp->[1][2][2];
  #$interp->[1][2][1] = (1-$dk)*$interp->[1][2][0] + $dk*$interp->[1][2][2];
  $interp->[2][1][1] = (1-$dk)*$interp->[2][1][0] + $dk*$interp->[2][1][2];

  #### interpolate to the center #####
  $interp->[1][1][1] = (1-$di)*$interp->[0][1][1] + $di*$interp->[2][1][1];
  $val = $interp->[1][1][1];
  return($val);
}
##############################################################################
#############################################################################
# use the min value surrounding gridpoints 
#############################################################################
sub minvalue {
    my($i,$j,$k) = @_;
    # gridpoint prev and next
    my($ip);
    my($in);
    my($jp);
    my($jn);
    my($kp);
    my($kn);
    # distance from prev gridpoint
    my($di);
    my($dj);
    my($dk);
    # interpolation matrix
    my($interp);
    my($val);
    # find gridpoints around point
    # check to see if outside grid bounries
    if ($i>0 and $i+1<$fft_grid->[0]) {
        $ip = int($i);
        $in = int($i+1);
        $di = $i-$ip;
    } else {
        $ip = $fft_grid->[0];
        $in = 0;
        if ($i>0) {
            $di = 1+$i;
        } else {
            $di = $i-$ip;
        }
    }
    if ($j>0 and $j+1<$fft_grid->[1]) {
        $jp = int($j);
        $jn = int($j+1);
        $dj = $j-$jp;}
    else {
        $jp = $fft_grid->[1];
        $jn = 0;
        if ($j>0) {
            $dj = 1+$j;
        } else {
            $dj = $j-$jp;
        }
    }
    if ($k>0 and $k+1<$fft_grid->[2]) {
        $kp = int($k);
        $kn = int($k+1);
        $dk = $k-$kp;
    } else {
        $kp = $fft_grid->[2];
        $kn = 0;
        if ($k>0) { 
            $dk = 1+$k;
        } else {
            $dk = $k-$kp;
        }
    }

    # get 8 values around the real
    $interp->[0][0][0] = $rho->[$ip][$jp][$kp];
    $interp->[2][0][0] = $rho->[$in][$jp][$kp];
    $interp->[0][2][0] = $rho->[$ip][$jn][$kp];
    $interp->[2][2][0] = $rho->[$in][$jn][$kp];
    $interp->[0][0][2] = $rho->[$ip][$jp][$kn];
    $interp->[2][0][2] = $rho->[$in][$jp][$kn];
    $interp->[0][2][2] = $rho->[$ip][$jn][$kn];
    $interp->[2][2][2] = $rho->[$in][$jn][$kn];

    $val = $interp->[0][0][0];
    if ($val>$interp->[2][0][0]){$val=$interp->[2][0][0];}
    if ($val>$interp->[0][2][0]){$val=$interp->[0][2][0];}
    if ($val>$interp->[2][2][0]){$val=$interp->[2][2][0];}
    if ($val>$interp->[0][0][2]){$val=$interp->[0][0][2];}
    if ($val>$interp->[2][0][2]){$val=$interp->[2][0][2];}
    if ($val>$interp->[0][2][2]){$val=$interp->[0][2][2];}
    if ($val>$interp->[2][2][2]){$val=$interp->[2][2][2];}

    return($val);
}
##############################################################################
#############################################################################
# quadradic interpolation for finding density value on grid
#############################################################################
sub quadinterp {
    my($i,$j,$k) = @_;
    # gridpoint prev and next
    my($ip);
    my($in);
    my($i2n);
    my($jp);
    my($jn);
    my($j2n);
    my($kp);
    my($kn);
    my($k2n);
    # distance from prev gridpoint
    my($di);
    my($dj);
    my($dk);
    # interpolation matrix
    my($interp);
    my($val);

    # find gridpoints around point
    if ($i-int($i) >= 0.5) {
        $ip = int($i);
    } else {
        $ip = int($i-1);
    }
    $in = int($ip+1);
    $i2n = int($ip+2);
    $di = $i-$ip;
    if ($j-int($j) >= 0.5) {
        $jp = int($j);
    } else {
        $jp = int($j-1);
    }
    $jn = int($jp+1);
    $j2n = int($jp+2);
    $dj = $j-$jp;
    if ($k-int($k) >= 0.5) {
        $kp = int($k);
    } else {
        $kp = int($k-1);
    }
    $kn = int($kp+1);
    $k2n = int($kp+2);
    $dk = $k-$kp;

  # get 27 values around the point of interest
    $interp->[0][0][0] = $rho->[$ip][$jp][$kp];
    $interp->[2][0][0] = $rho->[$in][$jp][$kp];
    $interp->[3][0][0] = $rho->[$i2n][$jp][$kp];
    $interp->[0][2][0] = $rho->[$ip][$jn][$kp];
    $interp->[2][2][0] = $rho->[$in][$jn][$kp];
    $interp->[3][2][0] = $rho->[$i2n][$jn][$kp];
    $interp->[0][3][0] = $rho->[$ip][$j2n][$kp];
    $interp->[2][3][0] = $rho->[$in][$j2n][$kp];
    $interp->[3][3][0] = $rho->[$i2n][$j2n][$kp];
    $interp->[0][0][2] = $rho->[$ip][$jp][$kn];
    $interp->[2][0][2] = $rho->[$in][$jp][$kn];
    $interp->[3][0][2] = $rho->[$i2n][$jp][$kn];
    $interp->[0][2][2] = $rho->[$ip][$jn][$kn];
    $interp->[2][2][2] = $rho->[$in][$jn][$kn];
    $interp->[3][2][2] = $rho->[$i2n][$jn][$kn];
    $interp->[0][3][2] = $rho->[$ip][$j2n][$kn];
    $interp->[2][3][2] = $rho->[$in][$j2n][$kn];
    $interp->[3][3][2] = $rho->[$i2n][$j2n][$kn];
    $interp->[0][0][3] = $rho->[$ip][$jp][$k2n];
    $interp->[2][0][3] = $rho->[$in][$jp][$k2n];
    $interp->[3][0][3] = $rho->[$i2n][$jp][$k2n];
    $interp->[0][2][3] = $rho->[$ip][$jn][$k2n];
    $interp->[2][2][3] = $rho->[$in][$jn][$k2n];
    $interp->[3][2][3] = $rho->[$i2n][$jn][$k2n];
    $interp->[0][3][3] = $rho->[$ip][$j2n][$k2n];
    $interp->[2][3][3] = $rho->[$in][$j2n][$k2n];
    $interp->[3][3][3] = $rho->[$i2n][$j2n][$k2n];

    ###### subroutine to fit a parabla to and interpolate y at x #####
    sub parabfit {
        #### intput is three x,f(x) points and value x for to solve for
        my($x1,$y1,$x2,$y2,$x3,$y3,$x) = @_;
        my($dx12);
        my($sq_dx12);
        my($dx13);
        my($sq_dx13);
        my($dy12);
        my($dy13);
        my($a);
        my($b);
        my($c);
        my($y);

        $dx12 = ($x1-$x2);
        $sq_dx12 = ($x1**2-$x2**2);
        $dx13 = ($x1-$x3);
        $sq_dx13 = ($x1**2-$x3**2);
        $dy12 = ($y1-$y2);
        $dy13 = ($y1-$y3);

        $b = ($dy13/$sq_dx13-$dy12/$sq_dx12)/($dx13/$sq_dx13-$dx12/$sq_dx12);
        $a = ($dy12-$b*$dx12)/$sq_dx12;
        $c = $y1-$a*$x1**2-$b*$x1;
        $y = $a*$x**2+$b*$x+$c;
        return($y);
    };

    ##### quadraticly interpolate along the x direction #######
    $interp->[1][0][0] = &parabfit(0.,$interp->[0][0][0],1.,$interp->[2][0][0],2.,$interp->[3][0][0],$di);
    $interp->[1][2][0] = &parabfit(0.,$interp->[0][2][0],1.,$interp->[2][2][0],2.,$interp->[3][2][0],$di);
    $interp->[1][3][0] = &parabfit(0.,$interp->[0][3][0],1.,$interp->[2][3][0],2.,$interp->[3][3][0],$di);
    $interp->[1][0][2] = &parabfit(0.,$interp->[0][0][2],1.,$interp->[2][0][2],2.,$interp->[3][0][2],$di);
    $interp->[1][2][2] = &parabfit(0.,$interp->[0][2][2],1.,$interp->[2][2][2],2.,$interp->[3][2][2],$di);
    $interp->[1][3][2] = &parabfit(0.,$interp->[0][3][2],1.,$interp->[2][3][2],2.,$interp->[3][3][2],$di);
    $interp->[1][0][3] = &parabfit(0.,$interp->[0][0][3],1.,$interp->[2][0][3],2.,$interp->[3][0][3],$di);
    $interp->[1][2][3] = &parabfit(0.,$interp->[0][2][3],1.,$interp->[2][2][3],2.,$interp->[3][2][3],$di);
    $interp->[1][3][3] = &parabfit(0.,$interp->[0][3][3],1.,$interp->[2][3][3],2.,$interp->[3][3][3],$di);

    ##### quadraticly interpolate along the y direction #######
    $interp->[1][1][0] = &parabfit(0.,$interp->[1][0][0],1.,$interp->[1][2][0],2.,$interp->[1][3][0],$dj);
    $interp->[1][1][2] = &parabfit(0.,$interp->[1][0][2],1.,$interp->[1][2][2],2.,$interp->[1][3][2],$dj);
    $interp->[1][1][3] = &parabfit(0.,$interp->[1][0][3],1.,$interp->[1][2][3],2.,$interp->[1][3][3],$dj);

    ##### quadraticly interpolate along the z direction #######
    $interp->[1][1][1] = &parabfit(0.,$interp->[1][1][0],1.,$interp->[1][1][2],2.,$interp->[1][1][3],$dk);

    $val = $interp->[1][1][1];
    return($val);
}
##############################################################################
#############################################################################
# cubic interpolation for finding density value on grid
#############################################################################
sub cubicinterp {
    my($i,$j,$k) = @_;
    # gridpoint prev and next
    my($i2p);
    my($ip);
    my($in);
    my($i2n);
    my($j2p);
    my($jp);
    my($jn);
    my($j2n);
    my($k2p);
    my($kp);
    my($kn);
    my($k2n);
    # distance from prev gridpoint
    my($di);
    my($dj);
    my($dk);
    # interpolation matrix
    my($interp);
    my($val);

    # find gridpoints around point
    $ip = int($i);
    $i2p = int($ip-1);
    $in = int($ip+1);
    $i2n = int($ip+2);
    $di = $i-$i2p;
    $jp = int($j);
    $j2p = int($jp-1);
    $jn = int($jp+1);
    $j2n = int($jp+2);
    $dj = $j-$j2p;
    $kp = int($k);
    $k2p = int($kp-1);
    $kn = int($kp+1);
    $k2n = int($kp+2);
    $dk = $k-$k2p;

    # get 48 values around the point of interest
    $interp->[0][0][0] = $rho->[$i2p][$j2p][$k2p];
    $interp->[1][0][0] = $rho->[$ip][$j2p][$k2p];
    $interp->[3][0][0] = $rho->[$in][$j2p][$k2p];
    $interp->[4][0][0] = $rho->[$i2n][$j2p][$k2p];
    $interp->[0][1][0] = $rho->[$i2p][$jp][$k2p];
    $interp->[1][1][0] = $rho->[$ip][$jp][$k2p];
    $interp->[3][1][0] = $rho->[$in][$jp][$k2p];
    $interp->[4][1][0] = $rho->[$i2n][$jp][$k2p];
    $interp->[0][3][0] = $rho->[$i2p][$jn][$k2p];
    $interp->[1][3][0] = $rho->[$ip][$jn][$k2p];
    $interp->[3][3][0] = $rho->[$in][$jn][$k2p];
    $interp->[4][3][0] = $rho->[$i2n][$jn][$k2p];
    $interp->[0][4][0] = $rho->[$i2p][$j2n][$k2p];
    $interp->[1][4][0] = $rho->[$ip][$j2n][$k2p];
    $interp->[3][4][0] = $rho->[$in][$j2n][$k2p];
    $interp->[4][4][0] = $rho->[$i2n][$j2n][$k2p];
    $interp->[0][0][1] = $rho->[$i2p][$j2p][$kp];
    $interp->[1][0][1] = $rho->[$ip][$j2p][$kp];
    $interp->[3][0][1] = $rho->[$in][$j2p][$kp];
    $interp->[4][0][1] = $rho->[$i2n][$j2p][$kp];
    $interp->[0][1][1] = $rho->[$i2p][$jp][$kp];
    $interp->[1][1][1] = $rho->[$ip][$jp][$kp];
    $interp->[3][1][1] = $rho->[$in][$jp][$kp];
    $interp->[4][1][1] = $rho->[$i2n][$jp][$kp];
    $interp->[0][3][1] = $rho->[$i2p][$jn][$kp];
    $interp->[1][3][1] = $rho->[$ip][$jn][$kp];
    $interp->[3][3][1] = $rho->[$in][$jn][$kp];
    $interp->[4][3][1] = $rho->[$i2n][$jn][$kp];
    $interp->[0][4][1] = $rho->[$i2p][$j2n][$kp];
    $interp->[1][4][1] = $rho->[$ip][$j2n][$kp];
    $interp->[3][4][1] = $rho->[$in][$j2n][$kp];
    $interp->[4][4][1] = $rho->[$i2n][$j2n][$kp];
    $interp->[0][0][3] = $rho->[$i2p][$j2p][$kn];
    $interp->[1][0][3] = $rho->[$ip][$j2p][$kn];
    $interp->[3][0][3] = $rho->[$in][$j2p][$kn];
    $interp->[4][0][3] = $rho->[$i2n][$j2p][$kn];
    $interp->[0][1][3] = $rho->[$i2p][$jp][$kn];
    $interp->[1][1][3] = $rho->[$ip][$jp][$kn];
    $interp->[3][1][3] = $rho->[$in][$jp][$kn];
    $interp->[4][1][3] = $rho->[$i2n][$jp][$kn];
    $interp->[0][3][3] = $rho->[$i2p][$jn][$kn];
    $interp->[1][3][3] = $rho->[$ip][$jn][$kn];
    $interp->[3][3][3] = $rho->[$in][$jn][$kn];
    $interp->[4][3][3] = $rho->[$i2n][$jn][$kn];
    $interp->[0][4][3] = $rho->[$i2p][$j2n][$kn];
    $interp->[1][4][3] = $rho->[$ip][$j2n][$kn];
    $interp->[3][4][3] = $rho->[$in][$j2n][$kn];
    $interp->[4][4][3] = $rho->[$i2n][$j2n][$kn];
    $interp->[0][0][4] = $rho->[$i2p][$j2p][$k2n];
    $interp->[1][0][4] = $rho->[$ip][$j2p][$k2n];
    $interp->[3][0][4] = $rho->[$in][$j2p][$k2n];
    $interp->[4][0][4] = $rho->[$i2n][$j2p][$k2n];
    $interp->[0][1][4] = $rho->[$i2p][$jp][$k2n];
    $interp->[1][1][4] = $rho->[$ip][$jp][$k2n];
    $interp->[3][1][4] = $rho->[$in][$jp][$k2n];
    $interp->[4][1][4] = $rho->[$i2n][$jp][$k2n];
    $interp->[0][3][4] = $rho->[$i2p][$jn][$k2n];
    $interp->[1][3][4] = $rho->[$ip][$jn][$k2n];
    $interp->[3][3][4] = $rho->[$in][$jn][$k2n];
    $interp->[4][3][4] = $rho->[$i2n][$jn][$k2n];
    $interp->[0][4][4] = $rho->[$i2p][$j2n][$k2n];
    $interp->[1][4][4] = $rho->[$ip][$j2n][$k2n];
    $interp->[3][4][4] = $rho->[$in][$j2n][$k2n];
    $interp->[4][4][4] = $rho->[$i2n][$j2n][$k2n];

    ###### subroutine to fit a cubic to 4 points and interpolate y at x #####
    sub cubefit {
        #### intput is three x,f(x) points and value x for to solve for
        my($x1,$y1,$x2,$y2,$x3,$y3,$x4,$y4,$x) = @_;
        my($dx12);
        my($sq_dx12);
        my($cu_dx12);
        my($dx13);
        my($sq_dx13);
        my($cu_dx13);
        my($dx14);
        my($sq_dx14);
        my($cu_dx14);
        my($dy12);
        my($dy13);
        my($dy14);
        my($a);
        my($b);
        my($c);
        my($d);
        my($a_dnom);
        my($y);

        $dx12 = ($x1-$x2);
        $sq_dx12 = ($x1**2-$x2**2);
        $cu_dx12 = ($x1**3-$x2**3);
        $dx13 = ($x1-$x3);
        $sq_dx13 = ($x1**2-$x3**2);
        $cu_dx13 = ($x1**3-$x3**3);
        $dx14 = ($x1-$x4);
        $sq_dx14 = ($x1**2-$x4**2);
        $cu_dx14 = ($x1**3-$x4**3);
        $dy12 = ($y1-$y2);
        $dy13 = ($y1-$y3);
        $dy14 = ($y1-$y4);
    
        $a = ($dy12*$dx13-$dy13*$dx12)*($sq_dx12*$dx14-$sq_dx14*$dx12);
        $a = $a-($dy12*$dx14-$dy14*$dx12)*($sq_dx12*$dx13-$sq_dx13*$dx12);
        $a_dnom = ($cu_dx12*$dx14-$cu_dx14*$dx12)*($sq_dx12*$dx13-$sq_dx13*$dx12);
        $a_dnom = $a_dnom-($cu_dx12*$dx13-$cu_dx13*$dx12)*($sq_dx12*$dx14-$sq_dx14*$dx12);
        $a = $a/$a_dnom;
        $b = $dy12*$dx13-$dy13*$dx12-$a*($cu_dx12*$dx13-$cu_dx13*$dx12);
        $b = $b/($sq_dx12*$dx13-$sq_dx13*$dx12);
        $c = ($dy12-$a*$cu_dx12-$b*$sq_dx12)/$dx12;
        $d = $y1-$a*$x1**3-$b*$x1**2-$c*$x1;
        $y = $a*$x**3+$b*$x**2+$c*$x+$d;
        return($y);
    };

    ##### cubicly interpolate along the x direction #######
    $interp->[2][0][0] = &cubefit(0.,$interp->[0][0][0],1.,$interp->[1][0][0],2.,$interp->[3][0][0],3.,$interp->[3][0][0],$di);
    $interp->[2][1][0] = &cubefit(0.,$interp->[0][1][0],1.,$interp->[1][1][0],2.,$interp->[3][1][0],3.,$interp->[3][1][0],$di);
    $interp->[2][3][0] = &cubefit(0.,$interp->[0][3][0],1.,$interp->[1][3][0],2.,$interp->[3][3][0],3.,$interp->[3][3][0],$di);
    $interp->[2][4][0] = &cubefit(0.,$interp->[0][4][0],1.,$interp->[1][4][0],2.,$interp->[3][4][0],3.,$interp->[3][4][0],$di);
    $interp->[2][0][1] = &cubefit(0.,$interp->[0][0][1],1.,$interp->[1][0][1],2.,$interp->[3][0][1],3.,$interp->[3][0][1],$di);
    $interp->[2][1][1] = &cubefit(0.,$interp->[0][1][1],1.,$interp->[1][1][1],2.,$interp->[3][1][1],3.,$interp->[3][1][1],$di);
    $interp->[2][3][1] = &cubefit(0.,$interp->[0][3][1],1.,$interp->[1][3][1],2.,$interp->[3][3][1],3.,$interp->[3][3][1],$di);
    $interp->[2][4][1] = &cubefit(0.,$interp->[0][4][1],1.,$interp->[1][4][1],2.,$interp->[3][4][1],3.,$interp->[3][4][1],$di);
    $interp->[2][0][3] = &cubefit(0.,$interp->[0][0][3],1.,$interp->[1][0][3],2.,$interp->[3][0][3],3.,$interp->[3][0][3],$di);
    $interp->[2][1][3] = &cubefit(0.,$interp->[0][1][3],1.,$interp->[1][1][3],2.,$interp->[3][1][3],3.,$interp->[3][1][3],$di);
    $interp->[2][3][3] = &cubefit(0.,$interp->[0][3][3],1.,$interp->[1][3][3],2.,$interp->[3][3][3],3.,$interp->[3][3][3],$di);
    $interp->[2][4][3] = &cubefit(0.,$interp->[0][4][3],1.,$interp->[1][4][3],2.,$interp->[3][4][3],3.,$interp->[3][4][3],$di);
    $interp->[2][0][4] = &cubefit(0.,$interp->[0][0][4],1.,$interp->[1][0][4],2.,$interp->[3][0][4],3.,$interp->[3][0][4],$di);
    $interp->[2][1][4] = &cubefit(0.,$interp->[0][1][4],1.,$interp->[1][1][4],2.,$interp->[3][1][4],3.,$interp->[3][1][4],$di);
    $interp->[2][3][4] = &cubefit(0.,$interp->[0][3][4],1.,$interp->[1][3][4],2.,$interp->[3][3][4],3.,$interp->[3][3][4],$di);
    $interp->[2][4][4] = &cubefit(0.,$interp->[0][4][4],1.,$interp->[1][4][4],2.,$interp->[3][4][4],3.,$interp->[3][4][4],$di);
    ##### cubically interpolate along the y direction #######
    $interp->[2][2][0] = &cubefit(0.,$interp->[2][0][0],1.,$interp->[2][1][0],2.,$interp->[2][3][0],3.,$interp->[2][4][0],$dj);
    $interp->[2][2][1] = &cubefit(0.,$interp->[2][0][1],1.,$interp->[2][1][1],2.,$interp->[2][3][1],3.,$interp->[2][4][1],$dj);
    $interp->[2][2][3] = &cubefit(0.,$interp->[2][0][3],1.,$interp->[2][1][3],2.,$interp->[2][3][3],3.,$interp->[2][4][3],$dj);
    $interp->[2][2][4] = &cubefit(0.,$interp->[2][0][4],1.,$interp->[2][1][4],2.,$interp->[2][3][4],3.,$interp->[2][4][4],$dj);
    ##### quadratically interpolate along the z direction #######
    $interp->[2][2][2] = &cubefit(0.,$interp->[2][2][0],1.,$interp->[2][2][1],2.,$interp->[2][2][3],3.,$interp->[2][2][4],$dk);

    $val = $interp->[2][2][2];
    return($val);
}
##############################################################################
