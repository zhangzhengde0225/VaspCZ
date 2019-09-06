#!/usr/bin/env perl
#;-*- Perl -*-

# This program fits two or matrices together using a linear least squares fit.
# The order of the fit is the first argument to the program.

use FindBin qw($Bin);
use lib "$Bin";
use Vasp;
use Math::Approx;

sub poly {
    my($n,$x) = @_;
    return $x ** $n;
}


@args=@ARGV;
(@args >= 5) || die "usage: dymfit.pl <order> <displacement 1> <matrix 1> <displacement 2> <matrix 2> ...\n";

$order = $args[0];
@args = @args[1..@args-1];

for ($i=0; $i<@args; $i=$i+2) {
    $j = 0;

    open (IN,$args[$i+1]) || die ("Can't open file $!");
    while (<IN>) {
        $_ =~ s/^\s+//;
        @temp = split(/\s+/,$_);
        $kk = @temp;
        for ($k=0; $k<@temp; $k++) {
            $matrix->[$i/2][$j][$k] = $temp[$k]*$args[$i];
        }
        $j++;
    }
    close(IN);
}
$jj = $j;

#print "Dimensions: $jj $kk\n";

for ($j=0; $j<$jj; $j++) {
    for ($k=0; $k<$kk; $k++) {
        for ($i=0; $i<@args; $i=$i+2) {
            $x{$args[$i]} = $matrix->[$i/2][$j][$k];
        }
        $x{0} = 0;
        $fit = new Math::Approx (\&poly, $order, %x);
        # $fit->print;

        printf  "%18.10f", ${$fit->{'A'}}[1]."  ";
    }
    print "\n";
}





