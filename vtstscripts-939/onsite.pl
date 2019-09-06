#!/usr/bin/env perl
#;-*- Perl -*-

if(@ARGV == undef) {
    open(OUTCAR, "OUTCAR") || die "Cannot open OUTCAR!";
} else {
    open(OUTCAR, $ARGV[0]) || die "Usage: onsite.pl OUTCAR or the path to OUTCAR";
}

# Identify the begin and end lines for onsite density matrix data and band occupancies
@beginn = split(/:/, `grep -n Iteration OUTCAR | tail -1`);
@endd = split(/:/, `grep -n 'Free energy of the ion-electron system (eV)' OUTCAR | tail -1`);

$lineno = 0;
while($line = <OUTCAR>) {
    $lineno += 1;
    if($beginn[0] <= $lineno && $lineno < $endd[0]) {
        if($line =~ /atom/ && $line !~ /atomic/) {
            @atomdata = split(/\s+/,$line);
            print join(" ",@atomdata[0,2])."\n";
        }
        if($line =~ /spin/) {
            $lineno += 7;
            @linedata = split(/\s+/,$line);
            print join(" ",@linedata).":  ";
            $line = <OUTCAR>;
            $totalspin = 0;
            for ($i = 1; $i<=5; $i++){
                $line = <OUTCAR>;
                @spindens = split(/\s+/,$line);
                print $spindens[$i]." ";
                $totalspin += $spindens[$i]**2;
            }
            # print " total=".int($totalspin+0.5)."\n";
            print " total=".$totalspin."\n";
            $line = <OUTCAR>;
        }
    }
}
close(OUTCAR);

