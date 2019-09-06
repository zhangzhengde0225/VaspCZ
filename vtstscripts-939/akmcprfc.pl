#!/usr/bin/env perl
# Spits out the number of force calls a given akmc process used.

if(@ARGV == 0)
{
    print "\nakmcprfc.pl Usage:\n    akmcprfc.pl [akmc process directory] <--verbose> <--rezip>\n\n";
    exit 0;
}

$rezip = 0;
$verbose = 0;

@dirs = ();
for($i = 0; $i < @ARGV; $i++)
{
    if(-e $ARGV[$i])
    {
        push(@dirs, $ARGV[$i]);
    }
    elsif($ARGV[$i] eq "--rezip")
    {
        $rezip = 1;
    }
    elsif($ARGV[$i] eq "--verbose")
    {
        $verbose = 1;
    }
}

for $procdir(@dirs)
{
    $totalfc = 0;
    if($verbose)
    {
        print "process directory $procdir:\n";
    }
    $fc += outiter("$procdir");
    $fc += outiter("$procdir/final");
    if(-e "$procdir/inter")
    {
        $fc += outiter("$procdir/inter");
    }
    $i = 0;
    while(-e "$procdir/inter$i")
    {
        $fc += outiter("$procdir/inter$i");
        $i += 1;
    }    
    if($verbose)
    {
        print "    dimer force calls: $fc\n";
    }
    $totalfc += $fc;
    $fc = outiter("$procdir/mins/min1/final");
    if(-e "$procdir/mins/min1/inter")
    {
        $fc += outiter("$procdir/mins/min1/inter");
    }
    $i = 0;
    while(-e "$procdir/mins/min1/inter$i")
    {
        $fc += outiter("$procdir/mins/min1/inter$i");
        $i += 1;
    }    
    $fc += outiter("$procdir/mins/min2/final");
    if(-e "$procdir/mins/min2/inter")
    {
        $fc += outiter("$procdir/mins/min2/inter");
    }
    $i = 0;
    while(-e "$procdir/mins/min2/inter$i")
    {
        $fc += outiter("$procdir/mins/min2/inter$i");
        $i += 1;
    }    
    if($verbose)
    {
        print "    minimization force calls: $fc\n";
    }
    $totalfc += $fc;
    if($verbose)
    {
        print "    total force calls: $totalfc\n";
    }
    else
    {
        print "$totalfc\n";
    }
}



sub outiter
{
    my($outdir) = @_;
    $iterstr = "";
    if (-e "$outdir/OUTCAR")
    {
        $iterstr = `tac $outdir/OUTCAR | grep Iteration | head -n 1`;
    }
    elsif (-e "$outdir/OUTCAR.gz")
    {
        system("gunzip", "$outdir/OUTCAR.gz");
        $iterstr = `tac $outdir/OUTCAR | grep Iteration | head -n 1`;
        if($rezip)
        {
            system("gzip", "$outdir/OUTCAR");
        }
    }
    elsif (-e "$outdir/OUTCAR.bz2")
    {
        system("bunzip2", "$outdir/OUTCAR.bz2");
        $iterstr = `tac $outdir/OUTCAR | grep Iteration | head -n 1`;
        if($rezip)
        {
            system("bzip2", "$outdir/OUTCAR");
        }
    }
    if($iterstr =~ m/(\d+)/) {
        return $1;
    }        
    else
    {
        return 0;
    }
}
