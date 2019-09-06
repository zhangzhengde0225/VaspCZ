#!/usr/bin/env perl
#;-*- Perl -*-

# This program extracts the dynamical matrix eigenvalues from the OUTCAR

@args = @ARGV;
@args <=1 || die "usage: dymgetfreq.pl <OUTCAR>\n";

$outfile = "OUTCAR";
if(@args==1){
    $outfile = $args[0];
}

$freqrawdata = `grep THz $outfile`;
@freqdata = split(/\n/,$freqrawdata);

open (OUT,">freq.dat");

$numfreq = scalar(@freqdata);
print "Writing $numfreq frequency values to freq.dat\n";

for ($i=0; $i<@freqdata; $i++) {
    $line = $freqdata[$i];
    $line =~ s/^\s+//g;
    @line = split(/\s+/,$line);
    if($line[1] eq "f"){
        printf OUT "%12.6f cm^{-1} ... 0 \n", $line[7];
    }else{
        printf OUT "%12.6f cm^{-1} ... 1 \n", $line[6];
    }
}

close (OUT);
