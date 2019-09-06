#!/usr/bin/env perl

#=========================================================================================================

    # vkdbaddpr.pl
    # Adds a directory to the 'database' directory tree given a pr000N process directory from an 
    # akmc.pl run.

#=========================================================================================================

    use strict;
    use FindBin qw($Bin);
    use lib "$Bin";
    use kdbutil;
    
#=========================================================================================================

    # Loop over each process directory given.
    for(my $dirIndex = 0; $dirIndex < $#ARGV + 1; $dirIndex ++)
    {

        # Flag that tells us this process is okay to add to the kdb.
        my $goodToGo = 1;
        
        # The current directory.
        my $dir = $ARGV[$dirIndex];
        print "$dir... ";
        
        # Get the status and quality from st.dat in case we need to give a status message because
        #  of missing files.
        my $status = stdatStatus("$dir/../st.dat", substr($dir, -6));
        my $quality = stdatQuality("$dir/../st.dat", substr($dir, -6));
       
        # Check to see if the appropriate files exist in the appropriate places.
        if(!-e "$dir/final/CENTCAR") 
        {
            print "$status $quality.\n";
            $goodToGo = 0;
        }
        if(!-e "$dir/mins/min1/final/CONTCAR" && $goodToGo == 1) 
        {
            print "$status $quality.\n";
            $goodToGo = 0;
        }
        if(!-e "$dir/mins/min2/final/CONTCAR" && $goodToGo == 1) 
        {
            print "$status $quality.\n";
            $goodToGo = 0;
        }
        my $mode_file_name = "";
        if(!-e "$dir/final/NEWMODECAR" && $goodToGo == 1) 
        {
            if(!-e "$dir/final/MODECAR") 
            {
                print "$status $quality.\n";
                $goodToGo = 0;
            }
            else
            {
                $mode_file_name = "$dir/final/MODECAR";
            }
        }
        else
        {
            $mode_file_name = "$dir/final/NEWMODECAR";
        }
        
    	my $saddle_file_name = "$dir/final/CENTCAR";
        my $min1_file_name = "$dir/mins/min1/final/CONTCAR";
        my $min2_file_name = "$dir/mins/min2/final/CONTCAR";

        if($goodToGo == 1)
        {
            system("kdbaddpr.pl $min1_file_name $min2_file_name $saddle_file_name $mode_file_name");       
        }

    }


















