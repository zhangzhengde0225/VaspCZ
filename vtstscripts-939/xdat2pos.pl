#!/usr/bin/env perl
#;-*- Perl -*-

use Cwd;
use FindBin qw($Bin);
use lib "$Bin";

@ARGV>1 || die "usage:\nxdat2pos.pl <0> <start step> <end step>\nxdat2pos.pl <1> <certain step number>\n";
if(@ARGV[0] == 0 && @ARGV < 3) { die "usage: xdat2pos.pl <0> <start step> <end step>\n"; }

# Get the length of the XDATCAR file. This is bad but works for now !!!
open TMP, "XDATCAR";
$nl = 0;
while(<TMP>){ $nl++; }
close TMP;

if(@ARGV[0] == 0){
    $NumStep = @ARGV[2] - @ARGV[1] + 1;
    bunchpos($NumStep); 
    makemovie($NumStep);
}
if(@ARGV[0] == 1) {
    $step = @ARGV[1];
    onepos($step);
}

#make a movie of the selected steps
sub makemovie {
    my ($NumStep) = @_;
    my ($scount, $poscar, $i, $j, $confile, $xyzfile, $movie, $amy);
    $movie = "movie.xyz";
    if(-e $movie) { unlink $movie; }
    for($i=0; $i<$NumStep; $i++){
        $scount = $ARGV[1] + $i;
        $poscar = "POSCAR".$scount.".out";
        $confile = $poscar.".con";
        $xyzfile = $poscar.".xyz";
        $j = `$Bin/pos2con.pl $poscar; $Bin/con2xyz.pl $confile`;
        if($i == 0) {
            $amy = `cat $xyzfile > $movie`;
        } else {
            $amy=`cat $xyzfile >> $movie`;
        }
        unlink $confile, $xyzfile;
    }
}

# Subroutine to get a series of steps from XDATCAR
sub bunchpos{
    my $NumStep = shift;
    my $scount = 0;
    # Set selected steps
    for($scount=0; $scount<$NumStep; $scount++) {
        $Step[$scount] = @ARGV[1] + $scount;
        onepos($Step[$scount]);
    }
}
# End of subroutine bunchpos

# Subroutine to get one step from XDATCAR
sub onepos {
    my ($step) = @_;
    # Get information about, number and types of atoms, box lengths ets.
    # from the POSCAR file and write out to the output POSCAR file
    open POS, "POSCAR" or die " NO POSCAR IN THIS DIRECTORY \n";
    $outfile = "POSCAR".$step.".out";
    open(OUT,">$outfile");

    $nhead=8;
    for($i=0; $i<$nhead; $i++) {
        $pos = <POS>;
#        print "line: ",$i," ",$pos;
        print OUT $pos;
        chomp($pos);
        $pos =~ s/^\s+//g;
        @pos = split /\s+/,$pos;
        if($i == 0){
            @elements = split /\s+/, $pos;
            $nel = @elements;
        }
        if($i == 5) {
            if($pos[0] =~ /^\d+$/) {
#                print "VASP4 POSCAR detected\n";
                @not[0..$nel-1] = @pos[0..$nel-1];
                    while($not[$k] != undef) {
                        $natoms += $not[$k++];
                    }   # Calculate the number of atoms
            } else {
                $atomtypeflag = 1; #check for vasp5 style POSCAR
                $nhead++;
#                print "VASP5 POSCAR detected\n";
                @elements = split /\s+/, $pos;
                $nel = @elements;
            }
        }
        if($i == 6){
            if($atomtypeflag == 1){
                @not[0..$nel-1] = @pos[0..$nel-1];
                while($not[$k] != undef) {
                    $natoms += $not[$k++];
                }   # Calculate the number of atoms
            }
        }
    }

#    print "Natoms: ",$natoms,"\n";

    $pos = undef;
    for($i=0; $i<$natoms; $i++) { $pos .= <POS>; }
#    print $pos;
    close POS;
    @pos = split /\n/, $pos;

  # check XDATCAR format
    open XDAT, "XDATCAR" or die " NO XDATCAR IN THIS DIRECTORY \n";
    for($i=0; $i<6; $i++) {
        $in = <XDAT>;}
    chomp($in);
    $in =~ s/^\s+//g;
    @in = split /\s+/,$in;
    $vasp5xdatcarflag = 0;
    if($in[0] =~ /^\d+$/) {
#        print "VASP4 XDATCAR detected\n";
    } else {
        $vasp5xdatcarflag = 1;
#        print "VASP5 XDATCAR detected\n";
    }
    close XDAT;

  # Get the right step from the XDATCAR
    open XDAT, "XDATCAR" or die " NO XDATCAR IN THIS DIRECTORY \n";
    $a = $step;
    if($vasp5xdatcarflag == 1) {
#        $st = 5 + ($natoms+1)*($a-1);   #GH: xdatcar format change
        $st = 8 + ($natoms+1)*($a-1);
    } else {
        $st = 7 + ($natoms+1)*($a-1);
    }
    $fn = $st + $natoms;
#    print "start: ",$st," end: ",$fn,"\n";

    die "THE SEARCH IS OUT OF BOUNDS\n" if($st > $nl || $fn > $nl);

    for($i=0; $i<$st-1; $i++) {
        $in = <XDAT>;
    }
    $in = <XDAT>;
    $j = 0; 
    
#    print "pos: \n";
#    print join("\n",@pos)."\n";

#    print "xdat: \n";
    for($i=$st+1; $i<=$fn; $i++) {
        $in = <XDAT>;
#        print $in;
        chomp($in);
        $in=~s/^\s+//g;
        @in=split /\s+/,$in;
        $p = $pos[$j];
        $j++ ;
        chomp($p);
        $p =~ s/^\s+//g;
        @p = split /\s+/,$p;
        printf OUT "%15.10f %15.10f %15.10f %5s %5s %5s \n",$in[0],$in[1],$in[2],$p[3],$p[4],$p[5];
#       $IN = <STDIN>;
    }
    close OUT;
    close XDAT;
    # End of subroutine onepos
}

