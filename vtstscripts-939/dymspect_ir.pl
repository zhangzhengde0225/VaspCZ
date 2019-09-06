#!/usr/bin/env perl
#;-*- Perl -*-

# Program uses the freq.dat to extract the normal modes and outputs the modes with the amplitude, to display the spectra.

@args = @ARGV;
@args <= 4 || die "usage: dymspect.pl <results.txt> <sigma: variance> <minimum frequency> <maximum frequency>\n";

#defaults
$freqfile = "results.txt";
$sigma = 5;
$x_min = 0;
$x_max = 3500;

if(@args>0){ $freqfile = $args[0]; }
elsif(@args>1){ $sigma = $args[1]; }
elsif(@args>2){ $x_min = $args[2]; }
else{ $x_max = $args[3]; }

$x_min >=0 || die "Note: minimum frequency must be >= 0.\n";
$x_max > $x_min || die "Note: maximum frequency must be > minimum frequency.\n";
$sigma >0 || die "Note: sigma (variance) must be > 0.\n";

# Removing all spaces to one blank
open (IN,$freqfile);
while (<IN>) {
    $_ =~ s/^\s+//g;
    $file1 .= $_;
}
close (IN);

# Placing all characters separated by a new line into array
@freqfile = split(/\n/,$file1);

$numimg = 0;
@freqValue;

# Placing each freqfile array value separated by spaces into new array
for ($i=0; $i<@freqfile; $i++) {
    @line1 = split(/\s+/,$freqfile[$i]);
#    if($line1[3] == 1) {
#        $numimg++;
#    	next;
#    }
#    $freqValue[$i-$numimg] = $line1[0];
    $freqValue[$i] = $line1[1];
    $intenseValue[$i] = $line1[2];
}

#print " Note: freq.dat has $numimg imaginary frequencies.\n";

 # Calculating the gaussian for each frequency data point.
@xValue;
for ($i=0; $i<@freqValue; $i++) {
    for ($j=$x_min; $j<=$x_max; $j++) {
        $x = ($j - $freqValue[$i]);
        $y = ( $x / $sigma)**2;
        $z = exp(-0.5*$y);
        $xValue[$j] += $z * $intenseValue[$i];
    }
}

print " Writing spectra to spect.dat.\n";
# print the x_max values and amplitude to spect.dat
open (OUT, '>spect.dat');
for ($i=0; $i<@xValue; $i++) {
    printf OUT "%-4d %9.5f\n", $i, $xValue[$i];
#    print OUT "\t $i \t $xValue[$i]\n";
} 
close (OUT); 

