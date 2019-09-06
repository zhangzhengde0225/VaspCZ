#!/usr/bin/env perl

use FindBin qw($Bin);
use lib "$Bin";
use Math::Matrix;

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
        $returnvec->[$i1] = $vec1->[$i1]+$vec2->[$i1];
    }
    return $returnvec;
}

sub vsub {
    my $vec1 = shift;
    my $vec2 = shift;

    my $i1;
    my $returnvec = [];
    for ($i1=0; $i1<@{$vec1}; $i1++) {
        $returnvec->[$i1] = $vec1->[$i1]-$vec2->[$i1];
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
    for ($i1=0; $i1<@{$matrix}; $i1++) {
        $value = 0;
        for ($j1=0;$j1<@{$matrix->[0]}; $j1++) {
            $value += $matrix->[$i1]->[$j1]*$vec->[$j1];
        }
        $returnvec->[$i1] = $value;
    }
    return $returnvec;
}

@args = @ARGV;
$filename = @args[0];

if (!defined($filename)) {
    $filename = "PPinput.txt";
}

open PPinp , "<",$filename;
$basismode = 0;
$cellmode = 0;
$basisvars = [];
$cellvars = []; # These two variables are used to store the information from the input file
$basisvarnames = [];
$cellvarnames = [];
$basisvarindex = -1;
$cellvarindex = -1;
while (<PPinp>) {
    if(lc(substr($_,0,6)) eq "basis:") {
        $basismode = 1;
        $cellmode = 0;
    } elsif (lc(substr($_,0,5)) eq "cell:") {
        $basismode = 0;
        $cellmode = 1;
    } else {
        if (/:/) { #colons follow variable names
            s/://g;	s/\s//g;  s/\n//g;
            if ($basismode) {
                $basisvarindex += 1;
                $basisvarnames->[$basisvarindex] = $_;
            } elsif ($cellmode) {
                $cellvarindex += 1;
                $cellvarnames->[$cellvarindex] = $_;
            }
        } else {
            s/\n//g;
            @linesplit = split(" ",$_);
            if ($basismode) {
                $len = scalar(@{$basisvars->[$basisvarindex]});
                if (@linesplit == 1) {
                    $basisvars->[$basisvarindex]->[$len] = $linesplit[0];
                } else {
                    $basisvars->[$basisvarindex]->[$len] = 0;
                }
            } elsif ($cellmode) {
                if (@linesplit == 6) {
                    for ($i = 0; $i < 6; $i++) {
                        $cellvars->[$cellvarindex]->[$i] = $linesplit[$i];
                    }
                } else {
                    print "Error while reading cell variable <",$cellvarnames->[$cellvarindex],">\n";
                }
            }
        }
    }
}
close(PPinp);

open output, ">ssnebproj.dat";

opendir(DIR,".") or die "couldn't open . ($!)\n";
@list=readdir(DIR);
closedir(DIR);
@directories = grep {-d && /^[0-9][0-9]$/i} @list;
@directories = sort {$a<=>$b} @directories;

# Get the baseline information from the first directory
$initbasis = [];
$finalbasis = [];
$t1i = [];
$t2i = [];
$t3i = [];
$t1f = [];
$t2f = [];
$t3f = [];
if (-e "$directories[0]/POSCAR") {
    open POSCAR, "<$directories[0]/POSCAR";
    scalar(<POSCAR>); 
    $_ = <POSCAR>; s/\n//; s/\s//g;
    $scale = $_;
    $_ = <POSCAR>; s/\n//;
    @linesplit = split(" ",$_);
    $t1i->[0] = $scale*$linesplit[0]; 
    $t1i->[1] = 0; 
    $t1i->[2] = 0;
    $_ = <POSCAR>; s/\n//;
    @linesplit = split(" ",$_);
    $t2i->[0] = $scale*$linesplit[0];
    $t2i->[1] = $scale*$linesplit[1];
    $t2i->[2] = 0;
    $_ = <POSCAR>; s/\n//;
    @linesplit = split(" ",$_);
    $t3i->[0] = $scale*$linesplit[0];
    $t3i->[1] = $scale*$linesplit[1];
    $t3i->[2] = $scale*$linesplit[2];

    $_ = <POSCAR>; s/\n//;
        @linesplit = split(" ",$_);
        $_ = $linesplit[0];
        s/\D//g;
        if ($_ == 0) {
            $_ = <POSCAR>; s/\n//;
            @linesplit = split(" ",$_);
        }
        $numatoms = 0;
        for ($i = 0; $i<@linesplit; $i++) {
            $numatoms += $linesplit[$i];
        }
    scalar(<POSCAR>);
    for ($i = 0; $i<$numatoms;$i++) {
        $_ = <POSCAR>; s/\n//;
        @linesplit = split(" ",$_);
        $initbasis->[$i*3] = $linesplit[0];
        $initbasis->[$i*3+1] = $linesplit[1];
        $initbasis->[$i*3+2] = $linesplit[2];
        for ($j=0; $j<@{$basisvars};$j++) {
            if (!defined($basisvars->[$j]->[3*$i])) {
                $basisvars->[$j]->[3*$i] = 0;
            }
            if (!defined($basisvars->[$j]->[3*$i+1])) {
                $basisvars->[$j]->[3*$i+1] = 0;
            }
            if (!defined($basisvars->[$j]->[3*$i+2])) {
                $basisvars->[$j]->[3*$i+2] = 0;
            }
        }
    }
    close(POSCAR);
} else {
    print "POSCAR file for the initial image not found, execution ending.";
    die;
}

# Get the baseline information from the last directory
if (-e "$directories[-1]/POSCAR") {
    open POSCAR, "<$directories[-1]/POSCAR";
    scalar(<POSCAR>); 
    $_ = <POSCAR>; s/\n//; s/\s//g;
    $scale = $_;
    $_ = <POSCAR>; s/\n//;
    @linesplit = split(" ",$_);
    $t1f->[0] = $scale*$linesplit[0]; 
    $t1f->[1] = 0; 
    $t1f->[2] = 0;
    $_ = <POSCAR>; s/\n//;
    @linesplit = split(" ",$_);
    $t2f->[0] = $scale*$linesplit[0];
    $t2f->[1] = $scale*$linesplit[1];
    $t2f->[2] = 0;
    $_ = <POSCAR>; s/\n//;
    @linesplit = split(" ",$_);
    $t3f->[0] = $scale*$linesplit[0];
    $t3f->[1] = $scale*$linesplit[1];
    $t3f->[2] = $scale*$linesplit[2];

    $_ = <POSCAR>; s/\n//; s/\D//g;
    @linesplit = split(" ",$_);
    if (@linesplit == 0) {
        $_ = <POSCAR>; s/\n//; s/\D//g;
    }
    scalar(<POSCAR>);
    for ($i = 0; $i<$numatoms;$i++) {
        $_ = <POSCAR>; s/\n//;
        @linesplit = split(" ",$_);
        $finalbasis->[$i*3] = $linesplit[0];
        $finalbasis->[$i*3+1] = $linesplit[1];
        $finalbasis->[$i*3+2] = $linesplit[2];
    }
    close(POSCAR);
} else {
    print "POSCAR file for the final image not found, execution ending.";
    die;
}

print output "image,";
for ($i=0; $i<=$basisvarindex; $i++) {
    print output "$basisvarnames->[$i],";
}
if ($cellvarindex >= 0) {
    $Ftot = []; # Holds the deformation for the entire transformation
    $Ftot->[0] = $t1f->[0]/$t1i->[0]; # xx
    $Ftot->[1] = $t2f->[1]/$t2i->[1]; # yy
    $Ftot->[2] = $t3f->[2]/$t3i->[2]; # zz
    $Ftot->[3] = ($t2f->[0]-$Ftot->[0]*$t2i->[0])/$t2i->[1]; #xy
    $Ftot->[4] = ($t3f->[1]-$Ftot->[1]*$t3i->[1])/$t3i->[2]; # yz
    $Ftot->[5] = ($t3f->[0]-$Ftot->[0]*$t3i->[0]-$Ftot->[3]*$t3i->[1])/$t3i->[2];
    print output "$cellvarnames->[0]";
    if ($cellvarindex > 0) {
        print output ",";
    }
    $cellvarmatrix = (new Math::Matrix($cellvars->[0]))->transpose();
    for ($i=1; $i<=$cellvarindex; $i++) {
        print output "$cellvarnames->[$i]";
        if ($cellvarindex > $i) {
            print output ",";
        }
        $cellvarmatrix = $cellvarmatrix->concat((new Math::Matrix($cellvars->[$i]))->transpose());
    }
    $cellvarmatrix = $cellvarmatrix->pinvert(); # Used to project the atom basis onto this basis set
    $tcellbegin = $cellvarmatrix->multiply((new Math::Matrix([1,1,1,0,0,0]))->transpose());
    $tcellend = $cellvarmatrix->multiply((new Math::Matrix($Ftot))->transpose());
}
print output "\n";

$basislist = [];
$Flist = [];
for ($i=0;$i<@directories;$i++) {
    if (-e "$directories[$i]/CONTCAR") {
        open POSCAR, "<$directories[$i]/CONTCAR";
    } elsif (-e "$directories[$i]/POSCAR") {
        open POSCAR, "<$directories[$i]/POSCAR";
    } else {
        die "No valid POSCAR or CONTCAR found in $directories[$i]\n";
    }

    scalar(<POSCAR>); 
    $_ = <POSCAR>; s/\n//; s/\s//g;
    $scale = $_;
    $_ = <POSCAR>; s/\n//;
    @linesplit = split(" ",$_);
    $t1->[0] = $scale*$linesplit[0];
    $t1->[1] = 0;
    $t1->[2] = 0;
    $_ = <POSCAR>; s/\n//;
    @linesplit = split(" ",$_);
    $t2->[0] = $scale*$linesplit[0];
    $t2->[1] = $scale*$linesplit[1]; $t2->[2] = 0;
    $_ = <POSCAR>; s/\n//;
    @linesplit = split(" ",$_);
    $t3->[0] = $scale*$linesplit[0];
    $t3->[1] = $scale*$linesplit[1];
    $t3->[2] = $scale*$linesplit[2];

    $_ = <POSCAR>; s/\n//; s/\D//g;
    @linesplit = split(" ",$_);
    if (@linesplit == 0) {
        $_ = <POSCAR>; s/\n//; s/\D//g;
    }
    scalar(<POSCAR>);
    for ($j = 0; $j<$numatoms;$j++) {
        $_ = <POSCAR>; s/\n//;
        @linesplit = split(" ",$_);
        $basislist->[$i]->[$j*3] = $linesplit[0];
        $basislist->[$i]->[$j*3+1] = $linesplit[1];
        $basislist->[$i]->[$j*3+2] = $linesplit[2];
    }

    $Flist->[$i]->[0] = $t1->[0]/$t1i->[0]; # xx
    $Flist->[$i]->[1] = $t2->[1]/$t2i->[1]; # yy
    $Flist->[$i]->[2] = $t3->[2]/$t3i->[2]; # zz
    $Flist->[$i]->[3] = ($t2->[0]-$Flist->[$i]->[0]*$t2i->[0])/$t2i->[1]; # xy
    $Flist->[$i]->[4] = ($t3->[1]-$Flist->[$i]->[1]*$t3i->[1])/$t3i->[2]; # yz
    $Flist->[$i]->[5] = ($t3->[0]-$Flist->[$i]->[0]*$t3i->[0]-$Flist->[$i]->[3]*$t3i->[1])/$t3i->[2]; #xz

    close(POSCAR);
}

if ($basisvarindex >= 0) {
    # Make sure there are no wrap around artifacts in the projects
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
    }

    $basisvarmatrix = (new Math::Matrix($basisvars->[0]))->transpose();
    for ($i=1;$i<=$basisvarindex;$i++) {
        $basisvarmatrix = $basisvarmatrix->concat((new Math::Matrix($basisvars->[$i]))->transpose());
    }
    $basisvarmatrix = $basisvarmatrix->pinvert(); # Used to project the atom basis onto this basis set
    $tbasisbegin = $basisvarmatrix->multiply((new Math::Matrix($basislist->[0]))->transpose());
    $tbasisend = $basisvarmatrix->multiply((new Math::Matrix($basislist->[@{$basislist}-1]))->transpose());

    $basisprojlist = [];
    for ($i=0; $i < @{$basislist}; $i++) {
        $tbasisvals = $basisvarmatrix->multiply((new Math::Matrix($basislist->[$i]))->transpose());
        for ($j = 0; $j <= $basisvarindex; $j++) {
            $value = ($tbasisvals->[$j][0]-$tbasisbegin->[$j][0])/($tbasisend->[$j][0]-$tbasisbegin->[$j][0]);
            #print output "$value,";
            $basisprojlist->[$i]->[$j] = $value;
        }
    }
}

if ($cellvarindex >= 0) {
    $cellprojlist = [];
    for ($i=0; $i < @{$Flist}; $i++) {
        $tcellvals = $cellvarmatrix->multiply((new Math::Matrix($Flist->[$i]))->transpose());
        for ($j = 0; $j <= $cellvarindex; $j++) {
            $value = ($tcellvals->[$j][0]-$tcellbegin->[$j][0])/($tcellend->[$j][0]-$tcellbegin->[$j][0]);
            $cellprojlist->[$i]->[$j] = $value;
        }
    }
}

for ($i = 0; $i < @{$basisprojlist}; $i++) {
    print output "$i,";
    for ($j=0; $j <= $basisvarindex; $j++) {
        print output $basisprojlist->[$i]->[$j],",";
    }
    for ($j=0; $j <= $cellvarindex; $j++) {
        print output $cellprojlist->[$i]->[$j];
        if ($j<$cellvarindex) {
            print output ",";
        } else {
            print output "\n";
        }
    }
}

close(output);
