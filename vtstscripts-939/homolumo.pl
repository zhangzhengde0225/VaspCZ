#!/usr/bin/env perl

if (@ARGV == undef){
open (OUTCAR, "OUTCAR") || die "Cannot open OUTCAR!";
}else{
open (OUTCAR, $ARGV[0]) || die " Usage: homo-lumo.pl OUTCAR or path to OUTCAR";
}

@beginn = split(/:/, `grep -n 'spin component 1' OUTCAR | tail -1`);
@endd = split(/:/, `grep -n 'soft charge-density along one line, spin component           1' OUTCAR | tail -1`);
@bands = split(/\s+/, `grep NBANDS OUTCAR`);
$NBANDS = $bands[$bands-1];
print "Number of Bands = " . $NBANDS . "\n" ;  
$lineno = 0;
$test = False;
while ($line = <OUTCAR>) {
    $lineno += 1;
    if (($beginn[0] <= $lineno) && ($lineno < $endd[0])) {
        if ($line =~ /k-point/) {
            chomp($line);
            print $line;
            $line = <OUTCAR>;
            $lineno += 1;
            $occ = 1.00000;
            $En = 0;
            $band = 0;
            for ($i = 1; $i <= $NBANDS;  $i++) {
                $line = <OUTCAR>;
                $lineno += 1;
                @line = split(/\s+/, $line);
                if (($line[3]== $occ)||($line[3] > 0.5)){
                    $occ = $line[3];
                    $En = $line[2];
                    $band = $line[1];
                } else {
                    $HOMO = $band;
                    $LUMO = $line[1];
                    $BandGap = $En - $line[2];
                    print " " ."\t" . "HOMO = ". $HOMO . "\t"."E = ".$En."\t" ."LUMO = ". $LUMO . "\t". "E = ".$line[2]."\t" ."BandGap = " .$BandGap . "\n";
                    $occ = $line[3];
                    $En = $line[2];
                    $band = $line[1];
                }
            }
        }
    }
}

