#                              -*- Mode: Perl -*- 
# $Basename: Approx.pm $
# $Revision: 1.1 $
# Author          : Ulrich Pfeifer
# Created On      : Wed Oct 25 10:50:38 1995
# Last Modified By: Ulrich Pfeifer
# Last Modified On: Tue Mar  3 15:14:30 1998
# Language        : Perl
# 
# (C) Copyright 1995,1998, Ulrich Pfeifer
# 

=head1 NAME

Math::Approx

=head1 METHODS

=head2 new (constructor)

    Math::Approx->new(\&poly, $degree, %data);

If the first argument to the constructor is a CODE reference, this is
used as the function to iterate over the data. Such a function must
take two arguments: The I<degree> and the I<x> value.

For interpolation with plain polynomials I<poly> can be defined as:

        sub poly {
            my($n,$x) = @_;
            return $x ** $n;
        }

If the first argument in the constructor is a FALSE value instead of a
CODE reference, then the above plain polynomial I<poly> is used as the
iterator function.

The second argument is the maximum degree which should be used for
interpolation. Degrees start with B<0>.

The rest of the arguments are treated as pairs of B<x> and B<y>
samples which should be approximated.

The constructor returns a Math::Approx reference.

=head2 approx

	$approximation->approx(17);

The method returns the approximated  B<y> value for the B<x> value
given as argument.

=head2 fit

 	$approximation->fit;

Returns the medim square error for the data points.

=head2 plot

 	$approximation->plot("tmp/app");

Prints all data pairs and the corresponding approximation pairs into
the filename given as argument. The output is suitable for usage with
gnuplot(1).


=head2 print

 	$approximation->print;

Prints information about the approximation on I<STDOUT>

=head1 EXAMPLE

        use Math::Approx;
        
        sub poly {
            my($n,$x) = @_;
            return $x ** $n;
        }
        
        for (1..20) {
            $x{$_} = sin($_/10)*cos($_/30)+0.3*rand;
        }
        
        $a = new Math::Approx (\&poly, 5, %x);
        $a->print;
        $a->plot("math-approx-demo.out");
        print "Fit: ", $a->fit, "\n";

=head1 SEE ALSO

gnuplot(1).

=head1 AUTHOR

Ulrich Pfeifer E<lt>F<pfeifer@wait.de>E<gt>

=cut


package Math::Approx;
use Math::Matrix;
use 5.004;
use strict;
use vars qw($VERSION);

# $Format: "$VERSION = sprintf '%5.3f', $ProjectVersion$;"$
$VERSION = sprintf '%5.3f', 0.2;

sub new {
    my $type = shift;
    my $func = shift || sub {my($n,$x)=@_; $x**$n;};
    my $degr = shift;
    die "Math::Approx->new : invalid degree [$degr]" unless $degr;

    my %data = @_;
    die "Math::Approx->new : empty data set" unless %data;

    my $self = {};
    my @m;
    my @x;

    $self->{'F'} = $func;       # function of two arguments
    $self->{'N'} = $degr;
    $self->{'D'} = \%data;

    for my $x (keys %data) {
        my $row = [];
        for my $n (0 .. $degr) {
            push @{$row}, &{$func}($n, $x);
        }
        push @x, $data{$x};
        push(@m, $row);
    }
    my $I = new Math::Matrix(@m);            # $I->print("Initial\n");
    my $T = $I->transpose->multiply($I);     # $T->print("Quadratic\n");
    my $x = new Math::Matrix(\@x);           # $x->print("X\n");
    my $tx = $I->transpose->
        multiply($x->transpose);             # $tx->print("TX\n");
    my $eq = $T->concat($tx);                # $eq->print("EQN\n");
    my $s = $eq->solve;                      # $s->print("SOLUTION\n");
                                             # $T->multiply($s)->print("TEST\n");
    $self->{'A'} = $s->transpose->[0];
    bless $self, $type;
}

sub approx {
    my $self = shift;
    my $x    = shift;
    my $func = $self->{'F'};
    my $degr = $self->{'N'};
    my $result;

    for my $n (0 .. $degr) {
        $result += &{$func}($n, $x) * $self->{'A'}->[$n];
    }

    $result;
}

sub fit {
    my $self = shift;
    my $result;
    my $n;

    for my $key (keys %{$self->{'D'}}) {
        $result += ($self->{'D'}->{$key}-$self->approx($key))**2;
        #print STDERR "## $result\n";
        $n++;
    }

    $result/$n;
}

sub print {
    my $self = shift;
    
    printf "Function: %s\n", $self->{'F'};
    printf "Degree:   %d\n", $self->{'N'};
    print  "Koeff: ", join(' ', @{$self->{'A'}}), "\n";
    print  "Fit: ",    $self->fit, "\n";
    print  "Data:\n";
    print  "     X          Y     Approximation\n";
    for my $key (sort {$a <=> $b} keys %{$self->{'D'}}) {
        printf("%15.10f %15.10f %15.10f\n", $key, $self->{'D'}->{$key}, 
               $self->approx($key));
    }
}

sub plot {
    my $self = shift;
    my $file = shift;
    open(OUT, ">$file") || die "Could not open $file: $!\n";
    
    print OUT "\n#data\n";
    for my $key (sort {$a <=> $b} keys %{$self->{'D'}}) {
        print OUT $key, ' ', $self->{'D'}->{$key}, "\n";
    }

    print OUT "\n#Approximation\n";
    for my $key (sort {$a <=> $b} keys %{$self->{'D'}}) {
        print OUT $key, ' ', $self->approx($key), "\n";
    }
    close(OUT);
    1;
}

1;
