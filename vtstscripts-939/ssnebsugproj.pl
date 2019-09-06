#!/usr/bin/env perl

use FindBin qw($Bin);
use lib "$Bin";

sub vdot {
    my $vec1 = shift;
    my $vec2 = shift;

    my $il;
    my $returnval = 0;
    my $len1;
    if (@{$vec1} < @{$vec2}) {
        $len1 = scalar(@{$vec1});
    } else {
        $len1 = scalar(@{$vec2});
    }
    for ($i1=0; $i1<$len1; $i1++) {
        $returnval += $vec1->[$i1]*$vec2->[$i1];
    }
    return $returnval;
}

sub vadd {
    my $vec1 = shift;
    my $vec2 = shift;

    my $i1;
    my $returnvec = [];
    for ($i1=0; $i1<@{$vec1}; $i1++) {
        $returnvec->[$i1] = $vec1->[$i1] + $vec2->[$i1];
    }
    return $returnvec;
}

sub vsub {
    my $vec1 = shift;
    my $vec2 = shift;

    my $i1;
    my $returnvec = [];
    for ($i1=0; $i1<@{$vec1}; $i1++) {
        $returnvec->[$i1] = $vec1->[$i1] - $vec2->[$i1];
    }
    return $returnvec;
}

sub vmag {
    my $vec1 = shift;
    return sqrt(vdot($vec1,$vec1));
}

sub matrixmult {
    my $matrix = shift;
    my $vec = shift;

    my $returnvec = [];
    my $i1,$j1;
    for ($i1=0;$i1 < @{$matrix}; $i1++) {
        $value = 0;
        for ($j1=0;$j1< @{$matrix->[0]}; $j1++) {
            $value += $matrix->[$i1]->[$j1]*$vec->[$j1];
        }
        $returnvec->[$i1] = $value;
    }
    return $returnvec;
}

open output, ">ssneb.sproj";

@args = @ARGV;
$angthreshold = @args[0];

if (!defined($angthreshold)) {
    $angthreshold = 45;
}
$costhreshold = cos($angthreshold* 3.14159/180.0);

opendir(DIR,".") or die "couldn't open . ($!)\n";
@list=readdir(DIR);
closedir(DIR);
@directories = grep {-d && /^[0-9][0-9]$/i} @list;
@directories = sort {$a<=>$b} @directories;

$basislist = [];
for ($i=0;$i<@directories;$i++) {
    if (-e "$directories[$i]/CONTCAR") {
        open POSCAR, "<$directories[$i]/CONTCAR";
    } elsif (-e "$directories[$i]/POSCAR") {
        open POSCAR, "<$directories[$i]/POSCAR";
    } else {
        die "No valid POSCAR or CONTCAR found in $directories[$i]\n";
    }

    scalar(<POSCAR>); 
    $_ = <POSCAR>;
    $_ = <POSCAR>;
    $_ = <POSCAR>;
    $_ = <POSCAR>;

    $_ = <POSCAR>; s/\n//;
    @linesplit = split(" ",$_);
    $_ = $linesplit[0];
    s/\D//g;
    if ($_ == 0) {
        $_ = <POSCAR>; s/\n//;
        @linesplit = split(" ",$_);
    }
    $numatoms = 0;
    for ($j = 0; $j<scalar(@linesplit); $j++) {
        $numatoms = $numatoms + $linesplit[$j];
    }
    scalar(<POSCAR>);
    for ($j = 0; $j<$numatoms;$j++) {
        $_ = <POSCAR>; s/\n//;
        @linesplit = split(" ",$_);
        $basislist->[$i]->[$j*3] = $linesplit[0];
        $basislist->[$i]->[$j*3+1] = $linesplit[1];
        $basislist->[$i]->[$j*3+2] = $linesplit[2];
    }
    close(POSCAR);
}

# Make sure there are no wrap around artifacts in the basis
$difflist = [];
for ($i=0; $i < @{$basislist}-1; $i++) {
    $diff = vsub($basislist->[$i+1],$basislist->[$i]);
    for ($j=0; $j < @{$diff}; $j++) {
        if (abs($diff->[$j]) > .5) {
            if ($basislist->[$i+1]->[$j] > $basislist->[$i]->[$j]) {
                $basislist->[$i+1]->[$j] = $basislist->[$i+1]->[$j]-1;
            } else {
                $basislist->[$i+1]->[$j] = $basislist->[$i+1]->[$j]+1;
            }
        }
    }
    $difflist->[$i] = vsub($basislist->[$i+1],$basislist->[$i]);
    #print @{$difflist->[$i]},"\n";
}
$anglelist = [];
$switchs = [];
for ($i=0; $i< @{$difflist}-1; $i++) {
    $anglelist->[$i] = vdot($difflist->[$i+1],$difflist->[$i])/(vmag($difflist->[$i+1]) * vmag($difflist->[$i]));
    print $anglelist->[$i],"\n";
    if ($anglelist->[$i] < $costhreshold) {
        print "$i is basis change point\n";
        $switches->[@{$switches}] = $i;
    }
}
#print @{$switches};
# average appropriate differences for the suggested basis change definitions
if (defined($switches)) {
    print "There are ",scalar(@{$switches}+1)," suggested variables.\n";
    $suggestedvecs = [];
    for ($i=0; $i <= @{$switches}; $i++) {
        if ($i == 0) {
            $startj = 0;
        } else {
            $startj = $switches->[$i-1] + 1;
        }
        if ($i == @{$switches}) {
            $endj = @{$difflist} - 1;
        } else {
            $endj = $switches->[$i];
        }
        $tempvec = [];
        for ($k=0; $k < @{$difflist->[0]}; $k++) {
            $tempvec->[$k] = 0.0;
        }
        for ($j=$startj; $j <= $endj; $j++) {
            for ($k=0; $k < @{$difflist->[0]}; $k++) {
                $tempvec->[$k] += $difflist->[$j]->[$k]/($endj-$startj+1);
            }
        }
        for ($k=0; $k <= @{$difflist->[0]}; $k++) {
            $suggestedvecs->[$i]->[$k] = $tempvec->[$k];
            print $tempvec->[$k]," ";
        }
        print "\n";
    }
}
