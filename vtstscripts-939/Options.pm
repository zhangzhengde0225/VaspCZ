#!/usr/bin/env perl

package Options;

use strict;
use warnings;

use Getopt::Std;

our($versionString,$usageString,$docString);

sub usage {
    my $errorString = $_[0];
    my $output;
    $output .= "Error: $errorString\n" if defined($errorString);
    if (not defined($errorString)) {
        $output .= "$versionString\n" if defined($versionString);
        $output .= "VTST Scripts: Version CVS HEAD\n" if not defined($versionString);
    }
    $output .= "usage: $0 $usageString\n" if defined($usageString);
    if (not defined($errorString)) {
        $output .= "$docString\n" if defined($docString);
    }
    print $output;
    exit 1;
}

sub parseOptions {
    my $optionString = $_[0];
    my %options;
    getopts($optionString,\%options) or &usage();
    &usage() if $options{h};
    return %options;
}

1;
