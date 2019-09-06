#!/usr/bin/env perl
#;-*- Perl -*-

# Usage instructions.
if (@ARGV < 1)
{
    print "\nmodemag.pl returns the magnitude of a modecar.\n\nusage:\n\n    modemag.pl MODECAR\n\n";
    exit;
}

# Open the modecars.
open(MODECAR, $ARGV[0]);

# Load the data.
@modedata = ();
while ($line = <MODECAR>)
{
    chomp($line);
    @bits = split(/ /, $line);
    for($i = 0; $i < @bits; $i++)
    {
        if ($bits[$i] =~ /^-?\d*\./)
        {
            push(@modedata1, $bits[$i]);
        }
    }
}

# Calculate the magnitude.
$sum = 0;
for($i = 0; $i < @modedata1; $i++)
{
    $sum += $modedata1[$i] * $modedata1[$i];
}
$mag = sqrt($sum);

print "$mag\n";








