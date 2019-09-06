#!/usr/bin/env perl
#;-*- Perl -*-

# Usage instructions.
if (@ARGV < 2)
{
    print "\ndotmode.pl returns the overlap between two modecars, between -1 and 1.\n\nusage:\n\n    dotmode.pl MODECAR1 MODECAR2\n\n";
    exit;
}

# Open the modecars.
open(MODECAR1, $ARGV[0]);
open(MODECAR2, $ARGV[1]);

# Load the data.
@modedata1 = ();
@modedata2 = ();
while ($line = <MODECAR1>)
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
while ($line = <MODECAR2>)
{
    chomp($line);
    @bits = split(/ /, $line);
    for($i = 0; $i < @bits; $i++)
    {
        if ($bits[$i] =~ /^-?\d*\./)
        {
            push(@modedata2, $bits[$i]);
        }
    }
}

if(not(@modedata1 eq @modedata2))
{
    print "\nThe modecars are not of equivalent length, exiting.\n\n";
    exit;
}

# Normalize.
$sum = 0;
for($i = 0; $i < @modedata1; $i++)
{
    $sum += $modedata1[$i] * $modedata1[$i];
}
$mag = sqrt($sum);
for($i = 0; $i < @modedata1; $i++)
{
    $modedata1[$i] /= $mag;
}
$sum = 0;
for($i = 0; $i < @modedata2; $i++)
{
    $sum += $modedata2[$i] * $modedata2[$i];
}
$mag = sqrt($sum);
for($i = 0; $i < @modedata2; $i++)
{
    $modedata2[$i] /= $mag;
}

# Calculate the dot product.
$sum = 0;
for($i = 0; $i < @modedata1; $i++)
{
    $sum += $modedata1[$i] * $modedata2[$i];
}

print "$sum\n";









