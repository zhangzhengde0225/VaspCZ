#!/usr/bin/env perl

use FindBin qw($Bin);
use lib "$Bin";

$Options::usageString = "[POSCAR] [Bin Size] (Density) (Outfile)";
$Options::docString = 'Outputs a radial distribution of the given POSCAR.';
use Options;
my %options = &Options::parseOptions("h");

@ARGV>=2 || die &Options::usage("Too few options given.");
$poscarfile = $ARGV[0];
$bin_size = $ARGV[1];
$rho0 = 0.0;
if(@ARGV>2) { $rho0 = $ARGV[2]; }
$outfile = "pdf.dat";
if(@ARGV>3) { $outfile = $ARGV[3]; }

use Vasp;

# find the pair-distance distribution
#sub pdf{
#  my ($poscarfile,$bin_size)=@_;
#  my ($k,$i,$j, $difference,$mag_difference,$bins,$index,$number);
#  my ($coordinates,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective);
#  my ($PI,$radius,$)
if($bin_size == 0) { die "Error in coord_number, binsize is $bin_size.\n"; }
$pi = 4*atan2 1, 1;

# read in the poscar file
($coordinates,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective)
  = read_poscar($poscarfile);

# find half the box length
$mindist = 0;
for($i=0; $i<3; $i++) {
    for($j=0; $j<3; $j++) {
        $v->[0][$j] = $basis->[$i][$j];
    }
    $vmag = magnitude($v,1);
    if($mindist==0) { $mindist = $vmag; }
    elsif($mindist>$vmag) { $mindist = $vmag; }
}
$mindist /= 2.0;

# count the pair distance for each atom
for ($k=0; $k<$total_atoms-1; $k++) {
    for ($i=$k+1; $i<$total_atoms; $i++) {
        for ($j=0; $j<3; $j++) {
            $difference->[0][$j] = pbc($coordinates->[$k][$j] - $coordinates->[$i][$j]);
        }
        $difference = dirkar($difference,$basis,$lattice,1);
        $mag_difference = magnitude($difference,1);
        if($mag_difference < $mindist) {
            $index = int($mag_difference/$bin_size);
            if(!(exists($bins->{$index}))){
                $bins->{$index} = 0;
            }
            $bins->{$index} = $bins->{$index}+2;
            if($index == 0){
                print "$k $i $mag_difference\n";
                for ($j=0; $j<3; $j++) {
                    print $coordinates->[$k][$j]," ",$coordinates->[$i][$j],"\n";
                }
            }
        }
    }
}

  # calculate the center of the particle
# GH: this will cause problems for a particle crossing the cell boundary
#  for ($i=0;$i<$total_atoms;$i++) {
#    for ($j=0;$j<3;$j++) {
#      $center[$j]=$center[$j]+$coordinates->[$i][$j];
#    }
#  }
#  for ($j=0;$j<3;$j++) {
#    $center[$j]=$center[$j]/$total_atoms;
#  }

print "--------------------------------------------------\n";
print " Bin size of $bin_size\n";

# max pair distance
@a = (sort {$a<=>$b} keys %$bins);
$bin_max = $a[-1];

# average density
$rmax = ($bin_max+0.5)*$bin_size;
#  $vol_tot=4./3.*$pi*$rmax**3;
#  $vol_tot=4./3.*$pi*($rmax/2)**3;
$vol_tot = volume($basis);
if($rho0 == 0) { $rho0 = $total_atoms/$vol_tot; }
#  $bulkrho0=0.061238;
#  $rho0=$bulkrho0;
#  print "rho0: ",$rho0,"\n";

print " Distance cutoff (1/2 cell): $mindist\n";
print " Maximum pair distance: $rmax\n";
print " Volume: $vol_tot\n";
print " Average density: $rho0\n";

# Calculate PDF

open(OUT,">$outfile");
#  $bin_tot=$bin_max+10;
$bin_tot = $bin_max;
$dr = $bin_size;
$dr3 = $dr**3;
$natoms_tot = 0;  # this is really N^2-N

for ($bin=0; $bin<=$bin_tot; $bin++) {
    if(exists($bins->{$bin})) {
        $natoms_tot += $bins->{$bin};
    }
}

for ($bin=0; $bin<=$bin_tot; $bin++) {
    $r = ($bin+0.5)*$dr;
    $natoms = 0;
    if(exists($bins->{$bin})) {
        $natoms = $bins->{$bin};
    }

# GH: old calculation
#    $atom_inside=0;
#    for ($bin=1;$bin<$bin_total;$bin++) {
#       # count the atom inside spherical r and calculate rho_0
#       for ($k=0;$k<$total_atoms-1;$k++) {
#           for ($j=0;$j<3;$j++) {
#              $difference->[0][$j]=pbc($coordinates->[$k][$j]-$center[$j]);
#            }
#            $difference=dirkar($difference,$basis,$lattice,1);
#            $mag_difference=magnitude($difference,1);
#            if($mag_difference < $radius){
#               $atom_inside=$atom_inside+1;
#            }
#        }
#       $rho_0=$atom_inside/(4/3*3.1415926*($radius)**3);

    # calculate pdf_val

#    $bin_vol=4./3.*$pi*(($bin+1)**3-$bin**3)*$dr3;
    $bin_vol = 4./3.*$pi*(($r+$dr/2.)**3-($r-$dr/2.)**3);
    $navg = $bin_vol*$rho0;
#    $rho=$natoms/$bin_vol;

#    $gr=$natoms*$vol_tot/$natoms_tot;   # natoms_tot is the number of 2*npairs
#    $gr=$natoms*$vol_tot/$total_atoms**2;
    $gr = $natoms/($total_atoms*$navg);

#    $pdf_val=$natoms/$radius-4*3.1415926*$radius*$rho_0;
#    $pdf_val=4.*$pi*$r*($rho-$rho0);
    $pdf_val = 4.*$pi*$r*$rho0*($gr-1.);

    printf OUT "%8.4f     %8.4f\n", $r,$pdf_val;
#    print "$r $natoms  $pdf_val  $gr  $navg $bin_vol\n";
}

print " PDF written to pdf.dat\n";
print "--------------------------------------------------\n";

#  return $bins;
#}
