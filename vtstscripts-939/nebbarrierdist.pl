#!/usr/bin/env perl
#;-*- Perl -*-

use FindBin qw($Bin) ;
use lib "$Bin";
use Vasp;

@args = @ARGV;
if (@args == 0) {
    opendir(DIR,".") or die "couldn't open . ($!)\n";
    @list = readdir(DIR);
    closedir(DIR);

    @directories = grep {-d && /^[0-9][0-9]$/i} @list;
    @directories = sort {$a<=>$b} @directories;
} else {
    @directories = @args;
}

$numdir = @directories;
#print "#Directories found: ".join(" ",@directories)."\n";

$dist_cum = 0;
for ($i=0; $i<@directories; $i++) {
    (-e "$directories[$i]/OUTCAR") || die "No OUTCAR in $directories[$i]!\n";

    $energy = `grep 'energy  w' $directories[$i]/OUTCAR|tail -1`;
    $energy =~ s/\s+$//g;
    @energy = split(/\s+/,$energy);
    $energy = $energy[@energy-1];
    if($i==0) { $energy0 = $energy; }
    $energy -= $energy0;

    if($i<($numdir-1)) {
        $dist = `grep 'NEB: distance' $directories[$i]/OUTCAR|tail -1`;
        $dist =~ s/\s+$//g;
        @dist = split(/\s+/,$dist);
        $dist = $dist[@dist-3];
    } else {
        $dist = $dist[@dist-2];
    }
    $dist_cum += $dist;

    $force = `grep 'NEB: projections' $directories[$i]/OUTCAR|tail -1`;
    @force = split(/\s+/,$force);
    $force = $force[@force-1];

    printf "%3i %12.6f %12.6f %12.6f %3i\n",$i,$dist_cum,$energy,$force,$directories[$i];
}
