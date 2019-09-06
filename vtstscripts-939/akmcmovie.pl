#!/usr/bin/env perl
#;-*- Perl -*-

use FindBin qw($Bin);
use lib "$Bin";
use Vasp;
#####################################################################################
# Make movies for Adaptive kinetic Monte Carlo Simulation 1.0                       #
# by Lijun Xu and Graeme Henkelman (Univ. of Texas at Austin)                       #
# (last modifed by Lijun on July04, 2006)                                           #
#####################################################################################
if ( $/ ne "\n" )
{
    print "Warning: input record separator(" . $/
      . ") is not the default (new line). Reset it.\n";
    $/ = "\n";    #now, chomp will remove the last "\n" in a string.
}

# Note: $Bin will be a global variable and thus do not use it for other variables   #
#$dummy=substr($Bin, -1, 1);
#print "$dummy\n";
#if($dummy eq "\/") {chop($Bin);}
#print "The \$Bin is ".$Bin."\n";
print "Script directory is " . $Bin . "\n";

$argnum = @ARGV;
( $argnum > 0 )
  || die
  " USAGE: akmcmovie.pl filename(akmc.dat or st.dat?) <1 or 0 (in quadruple?)> <n (saddle-state interpolation)> <1 or 0 (both good and bad processes?)>  Note: 0 is the default.\n";

$ConfigDefault    = "POSCAR";
$ConfigDefaultCon = "POSCAR.con";
$ConfigDefaultXYZ = "POSCAR.xyz";
$min              = "min";
$min1             = "min1";
$min2             = "min2";
$dynmat           = "dynmat";
$quench           = "mins";
$ststatus         = "status";
$numprst          = "numprst";
$EventTable       = "EventTable";
chomp( $stdir = `pwd` );
$stdir =~ s/^\s+//;
@stdir = split( /\//, $stdir );

$filename          = $ARGV[0];
$quadruple         = 0;
$all               = 0;
$NumInsertedImages = 0;
if ( $argnum > 1 )
{
    $quadruple = $ARGV[1];
}
if ( $argnum > 2 )
{
    $NumInsertedImages = $ARGV[2];
    if ( $NumInsertedImages < 0 ) { $NumInsertedImages = 0; }
}
if ( $argnum > 3 )
{
    $all = $ARGV[3];
}
open( INPUT, "<$filename" ) or die "Can't open input file: $filename\n";
$i = 0;
while ( $line = <INPUT> )
{
    $line =~ s/^\s+//;
    @line = split( /\s+/, $line );
    if ( $i == 0 )
    {

        #    if(lc($line[0]) eq "akmcstep") {
        if ( lc( $line[0] ) eq "step" )
        {
            $filetype = "akmc_step";
            print "Make a movie of akmc steps\n";
        }
        else
        {
            $filetype = "akmc_state";
            print "Make a movie of mechanisms in a state\n";
        }
        $i++;
        last;
    }
    print
      "wrong format: $filename must be in the format of either index.dat or st.dat\n";
}
close INPUT;
if ( $filetype eq "akmc_step" )
{
    $akmcdir = $stdir[-1];
    $movies  = "watch_" . $akmcdir . ".xyz";
    if ( -e "$movies" ) { system "rm $movies"; }

    #print "************************************\n";
    print "Read index file that records each KMC step ......";
    $index = ReadIndex($filename);
    print "done\n";

    #print "************************************\n";
    $totalsteps = ( keys %$index ) + 1;
    $index->{$totalsteps} = [ ( "2bcreated", "unknown" ) ];
    $prest                = "";
    @vistedSts            = ();
    sub numerically { $a <=> $b }
    for $step ( sort numerically keys %$index )
    {

        if ( lc( $index->{$step}[1] ) eq "new" )
        {
            $curst = $index->{$step}[2];
        }
        elsif ( lc( $index->{$step}[1] ) eq "repeat" )
        {
            $curst = $index->{$step}[2];
        }
        elsif ( lc( $index->{$step}[0] ) eq "2bcreated" )
        {
            $curst = ".";
        }
        else
        {
            die "unrecognizable st quality in $filename.\n";
        }
        print "step: $step curst: $curst\n";

        #print "$step: "."@{$index->{$step}}\n";
        $visted = 0;
        for $item (@vistedSts)
        {
            if ( $curst eq $item )
            {

                #already visited
                $visted = 1;
            }
        }
        if ($quadruple)
        {
            if ( !$visited )
            {
                $savedBC = $ENV{'VTST_BC'};
                $ENV{'VTST_BC'} = "none";
                $dummy1 =
                  `cd $curst; $Bin/quad.pl POSCAR out; $Bin/pos2con.pl out; $Bin/quad_con.pl out.con out4.con;`;
                $ENV{'VTST_BC'} = $savedBC;
                $dummy1 =
                  `cd $curst; $Bin/pos2con.pl out4.con; $Bin/pos2con.pl out4; $Bin/pos2con.pl out4.con; $Bin/pos2con.pl out4; mv out4.con POSCAR.con; $Bin/con2xyz.pl POSCAR.con`;
                push @vistedSts, $curst;
            }
            if ( $NumInsertedImages && ( $prest ne "" ) )
            {    # not the first step
                $target = spline_2sts(
                    $prest . "/out4",
                    $curst . "/out4",
                    $NumInsertedImages
                );
                system "cat $target >> $movies; rm $target";
            }
        }
        else
        {
            if ( !$visited )
            {
                $dummy1 =
                  `cd $curst; $Bin/pos2con.pl POSCAR; $Bin/con2xyz.pl POSCAR.con`;
                push @vistedSts, $curst;
            }
            if ( $NumInsertedImages && ( $prest ne "" ) )
            {    # not the first step
                $target = spline_2sts(
                    $prest . "/POSCAR",
                    $curst . "/POSCAR",
                    $NumInsertedImages
                );
                system "cat $target >> $movies; rm $target";
            }
        }
        $target = $curst . "/POSCAR.xyz";
        system "cat $target >> $movies";
        $prest = $curst;
    }
    system "gzip $movies";
    for $item (@vistedSts)
    {
        $dummy1 = `cd $item; rm -f out out.con out4 POSCAR.con POSCAR.xyz`;
    }
}
elsif ( $filetype eq "akmc_state" )
{
    $stdir        = $stdir[-1];
    $akmcdir      = $stdir[-2];
    $mechanisms   = $akmcdir . "_" . $stdir . ".tar";
    $ihafhgewt34s = $akmcdir . "_movie_of_" . $stdir;
    print "************************************\n";

    #print "Read the current st data ......";
    $curdir = ".";
    $st = ReadSt( $curdir, $filename );

    #print "done\n";
    #print "************************************\n";
    if ($quadruple)
    {
        $savedBC = $ENV{'VTST_BC'};
        $ENV{'VTST_BC'} = "none";
        $dummy2 =
          `$Bin/quad.pl POSCAR out; $Bin/pos2con.pl out; $Bin/quad_con.pl out.con out4.con;`;
        $ENV{'VTST_BC'} = $savedBC;
        $dummy2 =
          `$Bin/pos2con.pl out4.con; $Bin/pos2con.pl out4; $Bin/pos2con.pl out4.con; $Bin/pos2con.pl out4; mv out4.con POSCAR.con; $Bin/con2xyz.pl POSCAR.con`;
    }
    else
    {
        $dummy2 = `$Bin/pos2con.pl POSCAR; $Bin/con2xyz.pl POSCAR.con;`;
    }
    if ( -e "$ihafhgewt34s" ) { system "rm -r $ihafhgewt34s"; }
    system "mkdir $ihafhgewt34s";

    #print "st(numprst) ".$st->{$numprst}[0]."\n";
    for ( $i = 1 ; $i <= $st->{$numprst}[0] ; $i++ )
    {
        $prdir = DirName( "pr", $i );
        $output = $ihafhgewt34s . "/" . $prdir . ".xyz";

        #print "output: $output\n";
        #print "quality: $st->{$prdir}[1]\n";
        if ( lc( $st->{$prdir}[1] ) eq "good" )
        {
            print "working on process $prdir, a good saddle.\n";
            $dummy4 = `cat POSCAR.xyz >> $output`;
            if ($quadruple)
            {
                $savedBC = $ENV{'VTST_BC'};
                $ENV{'VTST_BC'} = "none";
                $dummy3 =
                  `cd $prdir; $Bin/quad.pl POSCAR_sp out; $Bin/pos2con.pl out; $Bin/quad_con.pl out.con out4.con;`;
                $dummy3 =
                  `cd $prdir; cd $st->{$prdir}[3]; $Bin/quad.pl POSCAR out; $Bin/pos2con.pl out; $Bin/quad_con.pl out.con out4.con;`;
                $ENV{'VTST_BC'} = $savedBC;
                $dummy3 =
                  `cd $prdir; $Bin/pos2con.pl out4.con; $Bin/pos2con.pl out4; $Bin/pos2con.pl out4.con; $Bin/pos2con.pl out4; $Bin/con2xyz.pl out4.con; mv out4.xyz POSCAR_sp.xyz`;
                $saddle = $prdir . "/POSCAR_sp.xyz";
                if ($NumInsertedImages)
                {
                    $target = spline_2sts( "out4", $prdir . "/out4",
                        $NumInsertedImages );
                    system "cat $target >> $output; rm $target";
                }
                system "cat $saddle >> $output";
                $dummy3 =
                  `cd $prdir; cd $st->{$prdir}[3];$Bin/pos2con.pl out4.con; $Bin/pos2con.pl out4; $Bin/pos2con.pl out4.con; $Bin/pos2con.pl out4; $Bin/con2xyz.pl out4.con; mv out4.xyz POSCAR.xyz`;
                $finalst = $prdir . "/" . $st->{$prdir}[3];
                if ($NumInsertedImages)
                {
                    $target = spline_2sts(
                        $prdir . "/out4",
                        $finalst . "/out4",
                        $NumInsertedImages
                    );
                    system "cat $target >> $output; rm $target";
                }
                $finalst .= "/POSCAR.xyz";
                system "cat $finalst >> $output";
            }
            else
            {
                $dummy3 =
                  `cd $prdir; $Bin/pos2con.pl POSCAR_sp; $Bin/con2xyz.pl POSCAR_sp.con; cd $st->{$prdir}[3]; $Bin/pos2con.pl POSCAR; $Bin/con2xyz.pl POSCAR.con`;
                if ($NumInsertedImages)
                {
                    $target = spline_2sts( "POSCAR", $prdir . "/POSCAR_sp",
                        $NumInsertedImages );
                    system "cat $target >> $output; rm $target";
                }
                $saddle = $prdir . "/POSCAR_sp.xyz";
                system "cat $saddle >> $output";
                $finalst = $prdir . "/" . $st->{$prdir}[3];
                if ($NumInsertedImages)
                {
                    $target = spline_2sts(
                        $prdir . "/POSCAR_sp",
                        $finalst . "/POSCAR",
                        $NumInsertedImages
                    );
                    system "cat $target >> $output; rm $target";
                }
                $finalst .= "/POSCAR.xyz";
                system "cat $finalst >> $output";
            }

            #print "$saddle $finalst\n";
            $dummy4 =
              `cd $prdir; rm -f out out.con out4 POSCAR_sp.xyz POSCAR_sp.con; cd $st->{$prdir}[3]; rm -f out out.con out4 POSCAR.xyz POSCAR.con`;
        }
        elsif (lc( $st->{$prdir}[1] ) eq "bad"
            && $all
            && ( -e "$prdir/mins/min1/final" )
            && ( -e "$prdir/mins/min2/final" ) )
        {
            print "working on process $prdir, a bad saddle\n";
            if ($quadruple)
            {
                $savedBC = $ENV{'VTST_BC'};
                $ENV{'VTST_BC'} = "none";
                $dummy5 =
                  `cd $prdir; $Bin/quad.pl POSCAR_sp out; $Bin/pos2con.pl out; $Bin/quad_con.pl out.con out4.con;`;
                $dummy5 =
                  `cd $prdir; cd mins/min1;$Bin/quad.pl POSCAR out; $Bin/pos2con.pl out; $Bin/quad_con.pl out.con out4.con;cd ../min2;$Bin/quad.pl POSCAR out; $Bin/pos2con.pl out; $Bin/quad_con.pl out.con out4.con;`;
                $ENV{'VTST_BC'} = $savedBC;
                $dummy5 =
                  `cd $prdir; $Bin/pos2con.pl out4.con; $Bin/pos2con.pl out4; $Bin/pos2con.pl out4.con; $Bin/pos2con.pl out4; $Bin/con2xyz.pl out4.con; mv out4.xyz POSCAR_sp.xyz`;
                $dummy5 =
                  `cd $prdir; cd mins/min1;$Bin/pos2con.pl out4.con; $Bin/pos2con.pl out4; $Bin/pos2con.pl out4.con; $Bin/pos2con.pl out4; $Bin/con2xyz.pl out4.con; mv out4.xyz POSCAR.xyz; cd ../min2;$Bin/pos2con.pl out4.con; $Bin/pos2con.pl out4; $Bin/pos2con.pl out4.con; $Bin/pos2con.pl out4; $Bin/con2xyz.pl out4.con; mv out4.xyz POSCAR.xyz`;
                $saddle  = $prdir . "/POSCAR_sp.xyz";
                $initst  = $prdir . "/mins/min1/POSCAR.xyz";
                $finalst = $prdir . "/mins/min2/POSCAR.xyz";
                system "cat $initst >> $output";

                if ($NumInsertedImages)
                {
                    $target = spline_2sts(
                        $prdir . "/mins/min1/out4",
                        $prdir . "/out4",
                        $NumInsertedImages
                    );
                    system "cat $target >> $output; rm $target";
                }
                system "cat $saddle >> $output";
                if ($NumInsertedImages)
                {
                    $target = spline_2sts(
                        $prdir . "/out4",
                        $prdir . "/mins/min2/out4",
                        $NumInsertedImages
                    );
                    system "cat $target >> $output; rm $target";
                }
                system "cat $finalst >> $output";
            }
            else
            {
                $dummy5 =
                  `cd $prdir; $Bin/pos2con.pl POSCAR_sp; $Bin/con2xyz.pl POSCAR_sp.con; cd mins/min1; $Bin/pos2con.pl POSCAR; $Bin/con2xyz.pl POSCAR.con; cd ../min2; $Bin/pos2con.pl POSCAR; $Bin/con2xyz.pl POSCAR.con;`;
                $saddle  = $prdir . "/POSCAR_sp.xyz";
                $initst  = $prdir . "/mins/min1/POSCAR.xyz";
                $finalst = $prdir . "/mins/min2/POSCAR.xyz";
                system "cat $initst >> $output";
                if ($NumInsertedImages)
                {
                    $target = spline_2sts(
                        $prdir . "/mins/min1/POSCAR",
                        $prdir . "/POSCAR_sp",
                        $NumInsertedImages
                    );
                    system "cat $target >> $output; rm $target";
                }
                system "cat $saddle >> $output";
                if ($NumInsertedImages)
                {
                    $target = spline_2sts(
                        $prdir . "/POSCAR_sp",
                        $prdir . "/mins/min2/POSCAR",
                        $NumInsertedImages
                    );
                    system "cat $target >> $output; rm $target";
                }
                system "cat $finalst >> $output";
            }

            #print "$saddle $finalst\n";
            $dummy6 =
              `cd $prdir; rm -f out out.con out4 out4.con POSCAR_sp.xyz POSCAR_sp.con; cd mins/min1; rm -f out out.con out4 out4.con POSCAR.xyz POSCAR.con; cd ../min2; rm -f out out.con out4 out4.con POSCAR.xyz POSCAR.con`;
        }
    }
    system "rm -f out out.con out4 out4.con POSCAR.con POSCAR.xyz";
    system
      "tar cvf $mechanisms $ihafhgewt34s/ > /dev/null; gzip $mechanisms; rm -r $ihafhgewt34s";
}
else
{
    print "unrecognizable filetype\n";
}

###############################################################################
# FUNCTIONS                                                                   #
###############################################################################
# ------------------------------------
# Directory name
# ------------------------------------
sub DirName
{
    my $pre     = shift;
    my $num     = shift;
    my $snum    = "";
    my $dirname = "";
    $snum = sprintf "%04d", $num;
    $dirname = $pre . $snum;
    return $dirname;
}

# -----------------------------------------------
# Read index [%index=ReadIndex($IndexFile)]
# -----------------------------------------------
sub ReadIndex
{
    my $indexfile = shift;
    my $line      = "";
    my @line      = ();
    my %index     = ();
    my $j         = 0;
    my $step      = 0;
    if ( -e $indexfile )
    {
        open( INDEX, "<$indexfile" );
    }
    else
    {
        return \%index;
    }
    while ( $line = <INDEX> )
    {
        if ( $step > 0 )
        {
            $line =~ s/^\s+//;
            @line = split( /\s+/, $line );
            if ( $line[0] ne "" )
            {
                $index{ $line[0] } = [ () ];
                for ( $j = 1 ; $j < @line ; $j++ )
                {
                    push @{ $index{ $line[0] } }, lc( $line[$j] );
                }
            }
        }
        $step++;
    }
    return \%index;
}

# -----------------------------------------------------------------------
# Read st information from st.dat [%st=ReadSt($curdir,$StFile)]
# -----------------------------------------------------------------------
sub ReadSt
{
    my ( $curdir, $StFile ) = @_;
    my $stfilename = $curdir . "/" . $StFile;
    my %st         = ();
    my $line       = "";
    my @line       = ();
    my ( $i, $j, $numprst );
    %st = ();
    if ( -e $stfilename )
    {
        open( ST, "<$stfilename" );
    }
    else
    {
        $st{"status"}          = [ ("NOST") ];
        $st{"numprst"}         = [ (0) ];
        $st{"NumSearchesLeft"} = [ ("unknown") ];
        return ( \%st );
    }
    $i = 0;
    while ( $line = <ST> )
    {
        $i++;
        $line =~ s/^\s+//;
        if ( $line eq "\n" ) { $i--; next; }    # skip empty lines
        @line = split( /\s+/, $line );
        if ( $i == 1 )
        {    # the first line must be st and st energy etc.
            $st{$curdir} = [ () ];
            for ( $j = 1 ; $j < @line ; $j++ ) {
                push @{ $st{$curdir} }, $line[$j];
            }
        }
        elsif ( $i == 2 )
        {
            $st{"status"} = [ () ];
            for ( $j = 1 ; $j < @line ; $j++ )
            {
                push @{ $st{"status"} }, lc( $line[$j] );
            }
        }
        elsif ( $i == 3 )
        {
            $st{"dynmat"} = [ () ];
            for ( $j = 1 ; $j < @line ; $j++ )
            {
                push @{ $st{"dynmat"} }, lc( $line[$j] );
            }
        }
        elsif ( $i == 4 )
        {
            $st{"NumSearchesLeft"} = [ () ];
            for ( $j = 1 ; $j < @line ; $j++ )
            {
                push @{ $st{"NumSearchesLeft"} }, lc( $line[$j] );
            }
        }
        elsif ( $i == 5 )
        {
            $st{"numprst"} = [ () ];
            for ( $j = 1 ; $j < @line ; $j++ )
            {
                push @{ $st{"numprst"} }, lc( $line[$j] );
            }
        }
        else
        {
            if ( $i == 6 && lc( $line[0] ) eq "process" )
            {    #an optional description line
                $st{"prinfolist"} = [ () ];
                for ( $j = 0 ; $j < @line ; $j++ )
                {
                    push @{ $st{"prinfolist"} }, lc( $line[$j] );
                }
                $i--;
                next;
            }
            $numprst = $i - 5;
            if ( !( lc( $line[0] ) eq DirName( "pr", $numprst ) ) )
            {
                die "in $StFile, st00xx has a wrong format.\n";
            }
            $st{ lc( $line[0] ) } = [ () ];
            for ( $j = 1 ; $j < @line ; $j++ )
            {
                push @{ $st{ lc( $line[0] ) } }, lc( $line[$j] );
            }
            if ( $st{ lc( $line[0] ) }[2] =~
                /^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$/ )
            {
                $st{ lc( $line[0] ) }[2] =
                  $st{ lc( $line[0] ) }[2] +
                  $st{$curdir}[0];    # convert barriers into energy.
            }
        }
    }
    if ( $st{"numprst"}[0] != $numprst )
    {
        die "The numbers of st00xx states conflict in $stfilename\n";
    }
    close ST;
    if ( exists( $st{""} ) ) { delete $st{""} }
    return ( \%st );
}

# -----------------------------------------------------------------------
# interploate images between two poscars: return a connected xyz file
# -----------------------------------------------------------------------
sub spline_2sts
{
    my ( $POSCAR1, $POSCAR2, $n_images ) = @_;
    my ( $coordinates1, $basis, $lattice, $num_atoms, $total_atoms,
        $selectiveflag, $selective, $description );
    my (
        $coordinates2, $basis2,         $lattice2,   $num_atoms2,
        $total_atoms2, $selectiveflag2, $selective2, $description2
    );
    my ( $interpolatedXYZ, $i, $j, $k, $l, $dr, $filename, $header, $dummy );
    (
        $coordinates1, $basis, $lattice, $num_atoms, $total_atoms,
        $selectiveflag, $selective, $description
    ) = read_poscar($POSCAR1);
    (
        $coordinates2, $basis2,         $lattice2,   $num_atoms2,
        $total_atoms2, $selectiveflag2, $selective2, $description2
    ) = read_poscar($POSCAR2);
    if ( $total_atoms2 == $total_atoms && @$num_atoms == @$num_atoms2 )
    {

        for ( $i = 0 ; $i < @$num_atoms ; $i++ )
        {
            if ( $num_atoms->[$i] != $num_atoms2->[$i] )
            {
                die
                  "In spline_2sts: $POSCAR1 and $POSCAR2 do not have the same number of atoms for each type.\n";
            }
        }
    }
    else
    {
        die
          "In spline_2sts: $POSCAR1 and $POSCAR2 do not have the same number of atoms.\n";
    }
    $dr = linear_division( $coordinates1, $coordinates2, $basis, $lattice,
        $total_atoms, $n_images );
    $header = "U6e8Go52j2H67aU";
    if ( -e "$header" )
    {
        $i = 1;
        while ($i)
        {
            $i++;
            $j = int( 100 * rand() );
            $header .= "$j";
            if ( !( -e "$header" ) ) { $i = 0; }
        }
    }
    for ( $i = 0 ; $i < $n_images ; $i++ )
    {
        $coordinates1 = vsum( $coordinates1, $dr, $total_atoms );
        set_bc( $coordinates1, $total_atoms );
        $filename->[$i] = $header . "$i";
        write_poscar(
            $coordinates1, $basis,       $lattice,
            $num_atoms,    $total_atoms, $selectiveflag,
            $selective,    $description, $filename->[$i]
        );
    }
    $interpolatedXYZ = "AcMilan";
    if ( -e "$interpolatedXYZ" ) { system("rm $interpolatedXYZ"); }
    for ( $i = 0 ; $i < $n_images ; $i++ )
    {
        $j = $filename->[$i];
        $k = $j . ".con";
        $l = $j . ".xyz";
        $dummy =
          `$Bin/pos2con.pl $j; $Bin/con2xyz.pl $k; cat $l >>$interpolatedXYZ; rm $j $k $l`;
    }
    return $interpolatedXYZ;
}

sub linear_division
{
    my ( $r1, $r2, $basis, $lattice, $total_atoms, $n ) = @_;
    my ( $differ, $dr, $i, $j );
    $differ = pbc_difference( $r2, $r1, $total_atoms );
    dirkar( $differ, $basis, $lattice, $total_atoms );
    for ( $i = 0 ; $i < $total_atoms ; $i++ )
    {
        for ( $j = 0 ; $j < 3 ; $j++ )
        {
            $dr->[$i][$j] = $differ->[$i][$j] / ( $n + 1 );
        }
    }
    kardir( $dr, $basis, $lattice, $total_atoms );
    return $dr;
}

sub spline_division
{
    my ( $r1, $r2, $basis, $lattice, $total_atoms, $n ) = @_;
    my ( $differ, $dr, $i, $j );
    print
      "Warning: you are trying to use the spline interpolation in akmcmovie.pl.\n";
    print
      "Unfortunately, it is not coded yet, so it actually does the linear interpolation.\n";
    $differ = pbc_difference( $r2, $r1, $total_atoms );
    dirkar( $differ, $basis, $lattice, $total_atoms );
    for ( $i = 0 ; $i < $total_atoms ; $i++ )
    {

        for ( $j = 0 ; $j < 3 ; $j++ )
        {
            $dr->[$i][$j] = $differ->[$i][$j] / ( $n + 1 );
        }
    }
    kardir( $dr, $basis, $lattice, $total_atoms );
    return $dr;
}
