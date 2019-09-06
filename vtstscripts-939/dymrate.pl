#!/usr/bin/env perl
#;-*- Perl -*-

# Setup and plot comparison in the "effective" activation energy using 
# classical results, zero-point energy corrected, quasi-quantum, and 
# quasi-quantum with tunneling approxamation      i.e.
# quantum partition functions for the vibrational modes.

# ----------------------------------------------------------------------------#

# Physical constants
$kB = 8.61738573e-5;            # Boltzmann's constant
$hbar = 6.46538e-2;             # Planck's constant over 2*pi
$pi = 3.141592653589793238464;  # pi
$c = 29979245800;               # speed of light cm
$w2v = 15633304363477.891;      # sqrt(eV/(Ang^2 amu) to s^-1
$w_invcm = 521.47;              # sqrt(eV/(Ang^2 amu) to cm^-1

# Read in the necessary info
die "  INPUT THE INITIAL AND FINAL TEMPERATURE, CLASSICAL ENERGY GAP(eV) 
    AND FILES WITH FREQUENCIES FOR THE MINIMUM AND SADDLE ! \n", if @ARGV < 5;

$ti = shift @ARGV;
$tf = shift @ARGV;
$dV = shift @ARGV;
$f1 = shift @ARGV;
$f2 = shift @ARGV;

# Read the file with the frequencies (angular**2) for the minimum
open FMIN , $f1 ;
while(<FMIN>) { $w1 .= $_; }
close FMIN;
@w1 = split /\n/ , $w1;
$n = @w1;
$sw1 = 0.0;
for($i=0; $i<$n; $i++) {
    if($w1[$i] < 0.0) {
        die "ONE MODE FOR THE MINIMUM IS NEGATIVE ! \n";
    } else {
        $w1[$i] = sqrt($w1[$i]);
        $sw1 += $w1[$i];
    }
}

# Read the file with the frequencies (angular**2) for the saddle
open FSAD , $f2 ;
while(<FSAD>) { $w2 .= $_; }
close FSAD;
@w2 = split /\n/ , $w2;
$n2 = @w2;
die "NOT THE SAME NUMBER OF MODES FOR MINIMUM AND SADDLE ! \n", if $n2 != $n;

$n2 = 0;
$sw2 = 0.0;
for($i=0; $i<$n; $i++) {
    if($w2[$i] < 0.0) {
        $n2++;
    } else {
        $w2[$i] = sqrt($w2[$i]);
        $sw2 += $w2[$i];
    }
}
  die "MORE THAN ONE UNSTABLE MODE FOR THE MINIMUM ! \n", if $n2 > 1;

# The zero-point offset in energy is
$zpe = 0.5*$hbar*($sw2 - $sw1);
$dVz = $dV + $zpe;

# Calculate the classical prefactor
  # first the minimum
#  $vcla = 1.0;
#  for($j=0; $j<$n; $j++) {
#      $vcla *= ($w1[$j]);
#  }
  # now the saddle
  # Assume here that the negative eigenvalue is first ... jump over it
$vcl = 1.0;
for($j=0; $j<$n; $j++) {
    if($w2[$j] < 0.0) {
        $vcl *= $w1[$j]*$w2v;
        next; }
    $vcl *= ($w1[$j]/$w2[$j]);
    }
$vcl = $vcl;

# Set the number of temperature values to 100 for now
$nt = 100;
$ib = 1.0/$ti;
$ie = 1.0/$tf;
$dt = ($ie-$ib)/($nt-1);
open OUT , ">eff_ea.dat";
open ROUT ,">rate.dat";
$t = $ib;

# Loop over temperature values
for($i=0; $i<$nt; $i++) {
# Quantum contribution for the minimum
    $a = 1;
    for($j=0; $j<$n; $j++) {
        $x = $hbar*$w1[$j]/(2*$kB)*$t ;
        $a *= (exp($x)-exp(-$x))/(2*$x);
    }
# Quantum contribution for the saddle
# Assume here that the negative eigenvalue is first ... jump over it
    $b = 1;
    for($j=1; $j<$n; $j++) {
        if($w2[$j] < 0.0) { next; }
        $x = $hbar*$w2[$j]/(2*$kB)*$t;
        $b *= (exp($x) - exp(-$x))/(2*$x);
    }
    $xtun = abs($hbar*sqrt(abs($w2[0]))/(2*$kB)*$t);

    $dVwig = $dV + (-$kB/$t*log($a/$b));
    if(sin($xtun) > 0.0) { 
        $dVtun = $dV + (-$kB/$t*(log($a/$b)+log($xtun/sin($xtun))));
        $TcFlag = 1.0;
    }else {
        $dVtun = 0.0;
        $TcFlag = 2**5120;
    }

    $kcl = $vcl*exp(-$dV*$t/$kB);
    $kclz = $vcl*exp(-$dVz*$t/$kB);
    $kwig = $vcl*exp(-$dVwig*$t/$kB);
    $ktun = $vcl*exp(-$dVtun*$t/$kB)*$TcFlag;

    printf OUT "%14.8f %14.8f %14.8f %14.8f %14.8f %14.8f \n",1/$t,$t*1000,$dV,$dVz,$dVwig,$dVtun;
    printf ROUT "%14.8f %14.8f %14.8e %14.8e %14.8e %14.8e \n",1/$t,$t*1000,$kcl,$kclz,$kwig,$ktun;
    $t += $dt;
  }
  close OUT;
  close ROUT;


