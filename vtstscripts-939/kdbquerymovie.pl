#!/usr/bin/env perl
#;-*- Perl -*-

#===================================================================================================

    use strict;
    use FindBin qw($Bin);
    use lib "$Bin";
    use File::Path;
    use File::Copy;
    use Storable qw(dclone);
    use Vasp;
    use kdbutil;

    # The number of frames to interpolate between each minimum and the saddle.
    use constant nInterpolate              => 16.0;

#===================================================================================================

    # Get the description.
    my $Desc = "";
    foreach(@ARGV)
    {
        $Desc .= $_ . " ";
    }

    # Get the matching directories.
    my @allDirs = getProcessDirs(kdbHome(), $Desc);

    # (Re)create the ./kdbmovie directory.
    rmtree("kdbmovie");
    mkdir("kdbmovie");
    
    # Loop over each matching directory...
    for(my $i = 0; $i < @allDirs; $i++)
    {
        # Get the current directory.
        my $dir = $allDirs[$i];
        
        # Open the .xyz files.
        open(MIN1, "<$dir/min1.xyz");
        open(MIN2, "<$dir/min2.xyz");
        open(SADDLE, "<$dir/saddle.xyz");
        open(INFO, ">kdbmovie/.info_" . $i);
        print INFO "$dir\n";
        close INFO;

        # Load the data in.
        my @type = ();
        my @min1 = ();
        my @min2 = ();
        my @saddle = ();
        my $numAtoms = <MIN1>;
        chomp $numAtoms;
        my $blank = <MIN1>;
        $blank = <MIN2>;
        $blank = <MIN2>;
        $blank = <SADDLE>;
        $blank = <SADDLE>;
        for(my $j = 0; $j < $numAtoms; $j++)
        {
            my $line;
            my @data;
            $line = <MIN1>;
            chomp $line;
            @data = split(" ", $line);
            $type[$j] = $data[0];
            $min1[$j] = [$data[1], $data[2], $data[3]];
            $line = <MIN2>;
            chomp $line;
            @data = split(" ", $line);
            $min2[$j] = [$data[1], $data[2], $data[3]];
            $line = <SADDLE>;
            chomp $line;
            @data = split(" ", $line);
            $saddle[$j] = [$data[1], $data[2], $data[3]];
        }
        
        # Close the input files.
        close MIN1;
        close MIN2;
        close SADDLE;

        ######################DEBUG
        print "$i\n";
        
        # Open the output file 
        open (MOVIE, ">kdbmovie/movie-" . $i . ".xyz");
        
        # The current frame number.
        my $frame = 0;

        # Start by looping over the interpolation from min1 to saddle:
        for(my $j = 0; $j < nInterpolate; $j++)
        {
            # Print the frame header.
            print MOVIE $numAtoms . "\n";
            print MOVIE "Frame $frame" . "\n";
            
            # Loop over and place the atoms for this frame.
            for(my $k = 0; $k < $numAtoms; $k++)
            {
                my $xDiff = $saddle[$k][0] - $min1[$k][0];
                my $yDiff = $saddle[$k][1] - $min1[$k][1];
                my $zDiff = $saddle[$k][2] - $min1[$k][2];
                my $xStep = $j * ($xDiff / nInterpolate);
                my $yStep = $j * ($yDiff / nInterpolate);
                my $zStep = $j * ($zDiff / nInterpolate);
                my $xSpot = $min1[$k][0] + $xStep;
                my $ySpot = $min1[$k][1] + $yStep;
                my $zSpot = $min1[$k][2] + $zStep;
                print MOVIE $type[$k] . " ";
                print MOVIE $xSpot . " ";
                print MOVIE $ySpot . " ";
                print MOVIE $zSpot . "\n";
            }
            
            # Increment the frame counter.
            $frame++;
        }
        
        # Now from saddle to min2.
        for(my $j = 0; $j < nInterpolate; $j++)
        {
            # Print the frame header.
            print MOVIE $numAtoms . "\n";
            print MOVIE "Frame $frame" . "\n";
            
            # Loop over and place the atoms for this frame.
            for(my $k = 0; $k < $numAtoms; $k++)
            {
                my $xDiff = $min2[$k][0] - $saddle[$k][0];
                my $yDiff = $min2[$k][1] - $saddle[$k][1];
                my $zDiff = $min2[$k][2] - $saddle[$k][2];
                my $xStep = $j * ($xDiff / nInterpolate);
                my $yStep = $j * ($yDiff / nInterpolate);
                my $zStep = $j * ($zDiff / nInterpolate);
                my $xSpot = $saddle[$k][0] + $xStep;
                my $ySpot = $saddle[$k][1] + $yStep;
                my $zSpot = $saddle[$k][2] + $zStep;
                print MOVIE $type[$k] . " ";
                print MOVIE $xSpot . " ";
                print MOVIE $ySpot . " ";
                print MOVIE $zSpot . "\n";
            }
            
            # Increment the frame counter.
            $frame++;
        }
        
        # Close the output.
        close MOVIE;

    # End loop over directories.
    }
    
#===================================================================================================




    









