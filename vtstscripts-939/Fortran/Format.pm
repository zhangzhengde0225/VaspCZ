package Fortran::Format;

use strict;
use warnings;

our $VERSION = '0.90';
use Data::Dumper;
our $DEBUG = 0;
use Carp;

=head1 NAME

Fortran::Format - Read and write data according to a standard Fortran 77 FORMAT 

=head1 SYNOPSYS

    use Fortran::Format;

    my $f = Fortran::Format->new("2('N: ',I4,2X)");
    print $f->write(1 .. 10);
    # prints the following:
    # N:    1  N:    2
    # N:    3  N:    4
    # N:    5  N:    6
    # N:    7  N:    8
    # N:    9  N:   10
    
    # if you don't want to save the format object, 
    # just chain the calls:
    Fortran::Format->new("2('N: ',I4,2X)")->write(1 .. 10);

=head1 DESCRIPTION

This is a Perl implementation of the Fortran 77 formatted input/output
facility.  One possible use is for producing input files for old Fortran
programs, making sure that their column-oriented records are rigorously
correct. Fortran formats may also have some advantages over C<printf> in some
cases: it is very easy to output an array, reusing the format as needed; and
the syntax for repeated columns is more concise. Unlike C<printf>, for good or
ill, Fortran-formatted fields B<never> exceed their desired width. For
example, compare

    printf "%3d", 12345;                            # prints "12345"
    print Fortran::Format->new("I3")->write(12345); # prints "***"

This implementation was written in pure Perl, with portability and correctness
in mind. It implements the full ANSI standard for Fortran 77 Formats (or at
least it should). It was not written with speed in mind, so if you need to
process millions of records it may not be what you need.

=head1 FORMATS

What follows is a very brief summary of Fortran formats. For a rigorous
description, see the ANSI standard. A format consists of a list of "edit
descriptors" or sublists of edit descriptors. Edit descriptors are separated by
commas, but the comma may be omitted if there's no ambiguity. Spaces and case 
are ignored, except within strings, so 'i 1 2' is the same as 'I12'.

=head2 Repeatable edit descriptors

The following edit descriptors may be repeated if they are preceded
by a number; for example, '3I4' is the same as 'I4,I4,I4' or 'I4I4I4' or
'I4,2I4'. Lists can be nested by using parentheses, so '2(I2I3)' is the same
as 'I2I3I2I3'. Most descriptors include a width I<w>. If the width is larger
than needed, the output is right-justified. If the width is not large enough,
the entire field is filled with asterisks.

=over

=item II<w>

=item II<w>.I<m>

An integer with width I<w>, and optionally a minimum number of digits
I<m> (adding zeroes on the left if needed).

=item FI<w>.I<d>

An fixed precision floating-point number with width I<w>,
and I<d> digits after the decimal point.

=item EI<w>.I<d>

=item EI<w>.I<d>EI<e>

=item DI<w>.I<d>

A number in exponential notation with width I<w>,
I<d> digits after the decimal point, and optionally I<e> digits after
the exponent.

=item GI<w>.I<d>

=item GI<w>.I<d>EI<e>

For values between 0.1 and 10^I<d>, format like I<F>. For values outside that
range, format like I<E>.

=item FI<w>

Treat the variable as Boolean and output either I<T> or I<F> in a field of 
width I<w>.

=item A

=item AI<w>

Insert a string variable. If the width is not specified, it outputs the 
entire string. If the width is smaller than the string, the string is 
truncated (instead of filling with asterisks).

=back

=head2 Non-repeatable edit descriptors

Most of the following descriptors don't output anything but act as control
strings. "Non-repeatable" descriptors can be repeated only by including them
in a repeated list within parentheses.

=over

=item 'I<string>'

Insert I<string> as is. Quotes may be escaped by doubling them; for example,
I<'Joe''s'> produces I<Joe's>.

=item I<n>HI<string>...

Insert The next I<n> characters after the H as is.

=item TI<c>

=item TLI<c>

=item TRI<c>

Move to position I<c> of the current record (T), or I<c> characters to the left
(TL), or I<c> characters to the right (TR).

=item I<n>X

Move I<n> characters to the right.

=item /

Move to the begining of the next record (the next line).

=item :

Stop producing output immediately if there are no more variables left to format.

=item S

=item SP

=item SS

Control whether the plus sign is included for positive numbers. Include it for
SP, do not include it for SS, and use the default (do not include) for S.

=item I<k>P

Scaling factor for output in exponential notation. By default, a number such
as 1.23 would be written as 0.123E+01. When a scaling factor I<k> is given,
the decimal point is shifted I<k> positions to the left and the exponent 
is decreased by I<k> orders of magnitude. With 1P the output would be 1.23E+00.

=back

=head1 METHODS

=over

=cut

=item new

    my $format = Fortran::Format->new($format_string);

Create a new format object. The string is parsed and compiled when the
object is constructed. Croaks if there is a syntax error.

=cut

# Fortran::Format->new($format_string)
# constructs and compiles a new format object
sub new {
    my $class = shift;
    $class = ref $class || $class;
    my $self = bless { 
        format => shift, 
        writer => Fortran::Format::Writer->new,
    }, $class;
    eval {
        $self->parse;
    };
    if ($@) {
        chomp $@;
        croak "Fortran::Format parse error: $@; pos=$self->{current_pos}\n",
            "$self->{format}\n", " " x $self->{current_pos}, "^\ncalled";
    }
    $self;
}

=item format

    my $format_string = $format->format;

Returns the format string used by the object.

=cut

# my $format_string = $format->format()
# returns the format string 
sub format {
    my $self = shift;
    $self->{format};
}

sub writer {
    my $self = shift;
    $self->{writer};
}

# $format->parse()
# tokenizes, parses, and compiles the format string 
sub parse {
    my $self = shift;
    my $s = $self->format;
    my $toks = $self->tokenize;
    print "$s\n" if $DEBUG;

    my $tree = Fortran::Format::RootList->build($self, 
        repeat => 1, writer => $self->writer);
    $self->{tree} = $tree;
    print Dumper $tree if $DEBUG;
}

=item write

    $output = $format->write(@data);

Formats the data. This is equivalent to the Fortran C<WRITE> statement,
except that it just returns the formatted string. It does not write 
directly to a file. Data items may be either scalar or array references 
(which can be nested).

For matrices (multidimensional arrays), the contents are formatted in
column-major order, same as in Fortran. For example,

    my $a = [[1,2],[3,4]];
    Fortran::Format->new('4I4')->write($a);

will print

    1   3   2   4

or 

    Fortran::Format->new('2I4')->write($a);

will print

    1   3
    2   4
    
This is effectively equivalent to transposing the matrix before
printing it in the row-major order that would be expected by
most non-Fortran programmers. This kludge is necessary 
to ensure that the output can be read properly by a Fortran 
program.

B<Note>: this is incompatible with Fortran::Format 0.5x, which
simply flattened the nested arrays, producing the output in row-major
order. Also note that the behavior is undefined if the nested array
is not rectangular. For example, [[1],[2,3]] will give strange results.

=cut

# my $output = $format->write(@data)
# executes the format and returns the output string
sub write {
    my ($self, @data) = @_;
    my $output;
    my $writer = $self->writer;
    $writer->begin;
    @data = _flatten(@data);
    while (@data) {
        my $data_count = @data;
        $self->{tree}->write(\@data);
        $writer->end_line;
        if (@data and @data == $data_count) { # make sure some data was used
            croak "infinite format scan for edit descriptor on writing";
        }
    }
    $writer->output;
}

# takes a list and "flattens" it by turning array references into list items
# example: flatten(1,[2,3],[4,[5,6[7]],8],9) returns (1,2,3,4,5,6,7,8,9)
sub _simple_flatten {
    my (@in) = @_;
    my @out;
    for my $item (@in) {
        if (ref $item eq 'ARRAY') {
            push @out, _simple_flatten(@$item);
        } else {
            push @out, $item;
        }
    }
    @out;
}

sub _flatten {
    my (@in) = @_;
    my @out;
    for my $item (@in) {
        if (ref $item eq 'ARRAY') {
            push @out, _colum_flatten($item);
        } else {
            push @out, $item;
        }
    }
    @out;
}

sub _transpose {
    my ($data, $offs, $size, @dims) = @_;
    unless (@dims) { return $data->[$offs] }
    my $n = pop @dims;
    my @ret;
    for my $i (0 .. $n-1) {
        push @ret, _transpose($data, $offs + $i*$size, $size*$n, @dims);
    }
    @ret;
}

sub _colum_flatten {
    my ($in) = @_;
    my @temp = _simple_flatten(@$in);
    my @dims;
    for (my $p = $in; ref $p; $p = $p->[0]) {
        push @dims, scalar @$p;
    }
    _transpose(\@temp, 0, 1, @dims);
}

=item read

    my (@results) = $format->read($fh, @input_list);

Read data from the filehandle $fh using the format ($fh can also be a string
instead of a filehandle). The input list is a list of array sizes: 1 for 
simple scalars, I<n> for simple arrays, and an array reference of dimensions
(such as [3,3]) for multidimensional arrays. For example,

    my ($i, $matrix, $j) = $format->read($fh, 1, [3,3], 2) 

will read one scalar, followed by a 3x3 matrix, followed by an array with size
two. B<Note>: this method should be called in list context!

The input list is needed because Fortran formats are reused automatically for
subsequent lines until all the variables are read.

Matrices are read in column-major order. See C<write> for details.

When reading, it is also possible to specify the length of the
resulting string variables by appending
"AI<length>". For example,

    my $s = $format->read($fh, '1A40')

will read the data I<into> a 40-character long string variable (this is
regardless of the field width specified in the format string itself). The
string will be padded with trailing spaces if needed to ensure that it is
exactly 40 characters long. This attempts to emulate Fortran's peculiar string
length semantics. It is needed if you want to read a string, write it back,
and be sure that you get the exact same output that you would get with 
Fortran.

For example,
    
    my $in = 'hello world';
    my $a5  = Fortran::Format->new('A5');
    my $a20 = Fortran::Format->new('A20');

    my ($s) = $a5->read($in, '1A10');

    print $a20->write($s);
    # prints "          hello     "

Notice that 1) C<$s> was padded with five space, to a length of ten
characters; 2) the output is right-justified to a total width of 20
characters.

Now, if we do this instead: 

    my ($s) = $a5->read($in, '1A3');
    print $a20->write($s);
    # prints "                 llo"

Five character are read from the left of the string ("hello"), but only the 
rightmost three are copied to the 3-character-long variable ("llo").

=cut

# possible way of specifying string length:
# my ($i, $matrix, $j) = $format->read($fh, 'A40' [3,A40], 1) 
# my ($i, $matrix, $j) = $format->read($fh, 'A40' '3A40', 1) 

# READ INTERFACE
# --> my ($i, $arr_ref, $j) = $format->read($fh, 1, [3,3], 1) 
sub read {
    my ($self, $input, @input_list) = @_;

    unless (wantarray) {
        croak "Fortran::Format->read should be called in list context";
    }

    $self->writer->begin(input_list => \@input_list);
    my $fh;
    if (ref $input) {
        $fh = $input;
    } else {
        open $fh, '<', \$input;
    }
    while ($self->writer->want_more) {
        $self->writer->begin_line;
        my $line = <$fh>; 
        chomp $line;
        $self->{writer}{input_line} = $line; # XXX
        $self->{tree}->read;  # read format once
        unless ($self->writer->read_something) {
            croak "infinite format scan for edit descriptor on reading";
        }
    }

    $self->writer->input_data;
}


# $format->tokenize()
# separate a string into tokens, which are stored internally by the object
# This version works for Hollerith strings
sub tokenize {
    my $self = shift;
    my $s = $self->format;

    my @chars = split '', $s;
    my $state = 0;
    my @toks;
    my ($tok, $len, $char);
    my $pos = 0;
    my $tok_pos  = $self->{current_pos} = 0;
    while (defined ($char = shift @chars) and ++$pos) {
        if ($state == 0) {
            $tok_pos = $pos - 1;
            $tok = uc $char;
            $state = 1, next if $char eq "'";       # begin string
            $state = 3, next if $char =~ /\d/;      # number
            $state = 5, next if $char =~ /[+-]/;    # sign
            next if $char eq ' ';                   # skip space
            next if $char eq ',';                   # skip comma
            push @toks, {tok => $tok, pos => $tok_pos};
        } elsif ($state == 1) {
            $tok .= $char;                          # string contents
            $state = 2, next if $char eq "'";       # quote
        } elsif ($state == 2) {
            $state = 1, next if $char eq "'";       # escaped quote
            push @toks, {tok => $tok, pos => $tok_pos};
            $state = 0, redo;                       # end of string
        } elsif ($state == 3) {
            $len = $tok, $state = 4, $tok = '', 
                next if uc $char eq 'H';            # begin H-string
            $tok .= $char, next if $char =~ /\d/;   # more digits
            next if $char eq ' ';                   # skip space
            push @toks, {tok => $tok, pos => $tok_pos};
            $state = 0, redo;                       # end of number
        } elsif ($state == 4) {
            if ($len-- == 0) {
                push @toks, {tok => "'$tok'", 
                    pos => $tok_pos};               # end of H-string
                $state = 0;
                redo;
            }
            $tok .= $char;                          # string contents
        } elsif ($state == 5) {
            $tok .= $char, next if $char =~ /\d/;   # more digits
            next if $char eq ' ';                   # skip space
            push @toks, {tok => $tok, pos => $tok_pos};
            $state = 0, redo;                       # end of number
        }
    }
    if ($state == 2 or $state == 3 or $state == 5) {
        push @toks, {tok => $tok, pos => $tok_pos};
    } elsif ($state == 1 or $state == 4) {
        $self->{current_pos} = length $self->format;
        die "unfinished string\n"; 
    }

    @toks = map { 
        if    ($_->{tok} eq '/') { $_->{tok} = "SLASH" }
        elsif ($_->{tok} eq ':') { $_->{tok} = "COLON" }
        $_ 
    } @toks;

    print Dumper \@toks if $DEBUG;
    $self->{toks} = \@toks;
}

sub get_tok {
    my ($self, $patt) = @_; 
    my $tok;
    if (! defined $patt || defined $self->peek_tok($patt)) {
        $tok = shift @{$self->{toks}};
        my $pos = $tok->{pos};
        $self->{current_pos} = $pos if $pos;
        $tok = $tok->{tok};
        print "  <$tok:$pos>\n" if $DEBUG and defined $tok;
        $self->{current_tok} = $tok;
    }
    $tok;
}

sub current_tok { $_[0]->{current_tok} }

sub peek_tok {
    my ($self, $patt) = @_;
    my $tok = $self->{toks}[0]{tok};
    defined $tok && $tok =~ /$patt/ ? $tok : undef;
}

package Fortran::Format::InputItem;

sub new {
    my ($class, %opts)  = @_;
    $class = ref $class || $class;

    my $dims = $opts{dimensions};
    $dims = [$dims] unless ref $dims;

    my $size = 1;
    my @idims;
    {no warnings; @idims = map { int } @$dims; }
    $size *= $_ for (@idims);

    my $last_dim = $dims->[-1];
    my $string_length;
    if ($last_dim =~ /^\d+A(\d+)$/) {
        $string_length = $1;
    }

    my $self = bless { 
        dimensions => \@idims, 
        size => $size,
        data => [],
        string_length => $string_length,
    }, $class;
    $self;
}

sub push_data {
    my ($self, $val) = @_;
    if ($self->{string_length}) {
        if (length $val > $self->{string_length}) {
            $val = substr $val, length($val) - $self->{string_length}; 
        } else {
            $val = sprintf "%-$self->{string_length}s", $val;
        }
    }
    push @{$self->{data}}, $val;
}

sub contents {
    my ($self) = @_;
    my $data = $self->{data};
    return undef if @$data < $self->{size};
    #use Data::Dumper; print "CONTENTS DATA:\n", Dumper $data;
    my $ret;
    if (@$data == 1) { # flatten scalars
        $ret = $data->[0];
    } else {
        $ret = _fill_array($data, 0, 1, @{$self->{dimensions}});
    }
    #print "CONTENTS RET:\n", Dumper $ret;
    $ret;
}

sub _fill_array {
    my ($data, $offs, $size, @dims) = @_;
    unless (@dims) { return $data->[$offs] }
    my $n = shift @dims;
    my @ret;
    for my $i (0 .. $n-1) {
        push @ret, _fill_array($data, $offs + $i*$size, $n*$size, @dims);
    }
    \@ret;
}

package Fortran::Format::Writer;

our $DEBUG = 0;

sub new {
    my $class = shift;
    $class = ref $class || $class;
    my $self = bless { }, $class;
}

sub begin {
    my ($self, %pars) = @_;
    $self->plus('');
    $self->bz(0);
    $self->scale(0);
    $self->reuse(0);
    $self->begin_line;
    $self->{input_data} = [];
    $self->{output}     = '';
    if ($pars{input_list}) {
        $self->{input_list} = [ map {
            Fortran::Format::InputItem->new(dimensions => $_)
        } @{$pars{input_list}} ];
        # XXX
    }
    #use Data::Dumper; print Dumper $self;
}

sub begin_line {
    my ($self) = @_;
    $self->{position} = 0;
    $self->{current_line}   = '';
    $self->{read_count} = 0;
}

sub end_line {
    my ($self) = @_;
    $self->{output} .= $self->{current_line} . "\n";
    $self->begin_line;
}

sub output {
    my ($self) = @_;
    $self->{output};
}

sub write {
    my ($self, $s) = @_;
    my $line = $self->{current_line};
    my $pos = $self->{position};

    if ($pos > length $line) { # need to pad with spaces
        $line .= " " x ($pos - length $line);
    }
    substr $line, $pos, length $s, $s;
    $self->{position} += length $s;
    $self->{current_line} = $line;
}

sub read {
    my ($self, $width) = @_;
    my $s = $self->{input_line};    
    no warnings;
    $s = substr($s, $self->{position}, $width);
    $s = sprintf "%-*s", $width, $s;
    print "extracted '$s'\n" if $DEBUG;
    $self->position(relative => $width);
    $s;
}

sub put {
    my ($self, $val) = @_;
    my $input = $self->{input_list}[0];
    print "putting '$val'\n" if $DEBUG;
    $self->{read_count}++;
    $input->push_data($val);
    my $ret = $input->contents;
    if (defined $ret) {
        print "full\n" if $DEBUG;    
        push @{$self->{input_data}}, $ret;
        shift @{$self->{input_list}};
    } else {
        print "not full yet\n" if $DEBUG;    
    }
}

sub input_data {
    my ($self) = @_;
    #use Data::Dumper; print "HI:\n", Dumper $self->{input_data};
    @{$self->{input_data}};
}

sub want_more {
    my ($self) = @_;
    scalar @{$self->{input_list}};
}

sub read_something {
    my ($self) = @_;
    $self->{read_count};
}

sub position {
    my ($self, $relative, $n) = @_;
    use Carp; confess unless @_ == 3;
    if ($relative eq 'relative') {
        $self->{position} += $n;
    } else {
        $self->{position}  = $n;
    }
    $self->{position} = 0 if $self->{position} < 0;
}

sub plus {
    my $self = shift;
    if (@_) { $self->{plus} = shift } else { $self->{plus} }
}

sub bz {
    my $self = shift;
    if (@_) { $self->{bz} = shift } else { $self->{bz} }
}

sub scale {
    my $self = shift;
    if (@_) { $self->{scale} = shift } else { $self->{scale} }
}

sub reuse {
    my $self = shift;
    if (@_) { $self->{reuse} = shift } else { $self->{reuse} }
}

package Fortran::Format::Node;

sub build {
    my $class = shift;
    my $tokenizer = shift;
    $class = ref $class || $class;
    my $self = bless { repeat => 1, @_ }, $class;
    $self->parse($tokenizer);
    $self;
}

sub new {
    my $class = shift;
    $class = ref $class || $class;
    my $self = bless { @_ }, $class;
}

sub writer {
    my $self = shift;
    $self->{writer};
}

sub write {
    my ($self, $data, $start) = @_;
    for (1 .. $self->{repeat}) {
        my $ret = $self->write_once($data, $start);
        return undef unless defined $ret; # ran out of data ?
        if (length $ret) {
            $self->writer->write($ret);
        }
    }
}

sub read {
    my ($self, $start) = @_;
    for (1 .. $self->{repeat}) {
        my $ret = $self->read_once($start);
        return undef unless defined $ret; 
    }
    1;
}

sub parse {} # do nothing

package Fortran::Format::Edit::Quote;

our @ISA = "Fortran::Format::Node";

sub parse {
    my ($self, $tokenizer) = @_;
    my $s = $tokenizer->current_tok;
    chop $s;
    substr $s, 0, 1, '';
    $self->{quoted_string} = $s;
}

sub write_once {
    my ($self, $data) = @_;
    return $self->{quoted_string};
}

package Fortran::Format::Edit::I;

our @ISA = "Fortran::Format::Node";

sub parse {
    my ($self, $tokenizer) = @_;
    my $tok = $tokenizer->get_tok('^\d+$') or die "expected \\d after I\n";
    $self->{width} = $tok;
    if ($tokenizer->get_tok('\.')) {
        $tok = $tokenizer->get_tok('^\d+$');
        defined $tok or die "expected \\d after I\\d.\n";
        $self->{min} = $tok;
    }
}

sub write_once {
    my ($self, $data) = @_;
    return undef unless @$data;
    my $i = int(shift @$data); 
    my $s = abs $i;
    if ($self->{min} and $self->{min} > length $s) { # add leading zeroes?
        my $zeroes = $self->{min}   - length $s;
        $s = "0" x $zeroes . $s;
    }
    if ($i < 0) {  # add negative sign?
        $s = "-$s";
    } else {
        $s = $self->writer->plus . $s;
    }
    if (defined $self->{min} and $self->{min} == 0 and $s == 0) {
        $s = ''; # zero with zero with must be output as blank
    }
    $s = sprintf "%$self->{width}s", $s; # right-justify
    if (length $s > $self->{width}) { # too wide?
        $s = "*" x $self->{width};
    }
    $s;
}

sub read_once {
    my ($self) = @_;
    return undef unless $self->writer->want_more;
    my $s = $self->writer->read($self->{width});
    if ($s =~ /^ *-?[\d ]+$/) {
        $s =~ s/^ +//;
        if ($self->writer->bz) {
            $s =~ s/ /0/g;
        } else {
            $s =~ s/ //g;
        }
        no warnings;
        my $i = int($s);
        #print "I parsed '$i'\n";
        $self->writer->put($i);
    } else {
        die "invalid integer '$s'\n";
    }
    1;
}

package Fortran::Format::Edit::F;

our @ISA = "Fortran::Format::Node";

sub parse {
    my ($self, $tokenizer) = @_;
    my $tok = $tokenizer->get_tok('^\d+$') or die "expected \\d after F\n";
    $self->{width} = $tok;
    $tokenizer->get_tok('^\.$') or die "expected . after F\\d\n";
    $tok = $tokenizer->get_tok('^\d+$'); 
    defined $tok or die "expected \\d after F\\d.\n";
    $self->{precision} = $tok;
}

sub write_once {
    my ($self, $data) = @_;
    return undef unless @$data;
    my $f = shift @$data; 
    $f *= 10 ** ($self->writer->scale);
    my $s = sprintf "%.$self->{precision}f", abs $f;
    if ($f < 0.0 and $s =~ /[1-9]/) { 
        # must only include negative sign for non-zero output
        $s = "-$s";
    } else {
        $s = $self->writer->plus . $s;
    }
    if ($self->{precision} == 0) {
        $s .= '.'; # must include decimal point even for Fn.0
    }
    $s = sprintf "%$self->{width}s", $s; # right-justify

    # Remove optional zero if width is too big by one
    $s =~ s/^([+-]?)0.(\d)/$1.$2/ if length $s == $self->{width} + 1; 
    if (length $s > $self->{width}) { # too wide?
        $s = "*" x $self->{width};
    }

    $s;
}

sub read_once {
    my ($self) = @_;
    return undef unless $self->writer->want_more;
    my $s = $self->writer->read($self->{width});
    my $f;

    if ($s =~ /^ *-?(?:[\d ]*\.?[\d ]*)$/ and $s =~ /\d/) {
        $s =~ s/^ +//; # remove leading spaces
        if ($self->writer->bz) {
            $s =~ s/ /0/g;
        } else {
            $s =~ s/ //g;
        }
        unless ($s =~ /\./) {
            substr $s, length($s) - $self->{precision}, 0, '.';
        }
        #no warnings;
        $f = $s / 10**($self->writer->scale);
        #print "F parsed '$i'\n";
        #$self->writer->put($i);
    } elsif ($s =~ /^[ .]*$/) {
        $f = 0;
    } else {
        die "invalid F number'$s'\n";
    }
    $self->writer->put($f);
    1;
}

package Fortran::Format::Edit::D;

our @ISA = "Fortran::Format::Node";

sub parse {
    my ($self, $tokenizer) = @_;
    $self->{width} = $tokenizer->get_tok('^\d+$') 
        or die "expected \\d after [DE]\n";
    $tokenizer->get_tok('\.') or die "expected . after [DE]\\d\n";
    my $tok = $tokenizer->get_tok('^\d+$'); 
    defined $tok or die "expected \\d after [DE]\\d.\n";
    $self->{precision} = $tok;
}

sub write_once {
    my ($self, $data) = @_;
    return undef unless @$data;
    my $s;  # working string

    my $d = shift @$data; 

    # shorthand
    my $scale     = $self->writer->scale;
    my $width     = $self->{width};
    my $precision = $self->{precision};
    my $exp_width = $self->{exp_width} || 0;

    # get exponent
    my $spf = sprintf "%.3E", $d;
    my ($exp) = $spf =~ /E(.*)/g;  # maybe floor log10 abs is faster?

    # normalize to "apparent" magnitude
    my $dnorm = abs $d * 10**($scale - $exp - 1);

    # validate scale factor range (from standard, 13.5.9.2.2)
    if ($scale <= 0 and -$precision < $scale
        or $scale > 0 and ($precision + 2) > $scale) {

        # apply scale factor
        $exp += -$scale + 1 if ($d != 0.0);
        $precision += -$scale + 1 if ($scale > 0);

        if ( !$exp_width ) { # calculate default exp. width
            $exp_width = (abs $exp > 99) ? 3 : 2;
        }

        # format exponent
        my $exp_s = sprintf "%+0*d", $exp_width + 1, $exp;
        if ($self->{exp_width} or $exp_width != 3) { # add optional E
            $exp_s = $self->exp_char . "$exp_s";
        }

        # proceed if exponent didn't overflow
        if (length $exp_s <= $exp_width + 2) {
            # format string (at last!)
            $s = sprintf "%.${precision}f$exp_s", $dnorm;

            # add sign if needed
            if ($d < 0.0 and $s =~ /[1-9]/) { 
                # must only include negative sign for non-zero output
                $s = "-$s";
            } else {
                $s = $self->writer->plus . $s;
            }

            # must include decimal point even for Fn.0
            $s =~ s/(\d)(E?[+-])/$1.$2/ unless ($s =~ /\./);

            # right-justify
            $s = sprintf "%${width}s", $s; 

            # Remove optional zero if width is too big by one
            $s =~ s/^([+-]?)0.(\d)/$1.$2/ if length $s == $width + 1; 

            # make sure final result did not overflow
            $s = undef if length $s > $width;
        }
    } 
    $s || "*" x $width;
}

sub exp_char { "D" }


package Fortran::Format::Edit::E;

our @ISA = "Fortran::Format::Edit::D";

sub parse {
    my ($self, $tokenizer) = @_;
    $self->SUPER::parse($tokenizer); # mostly similar to D
    if ($tokenizer->get_tok('^E$')) {
        $self->{exp_width} = $tokenizer->get_tok('^\d+$') 
            or die "expected \\d after E\\d.\\dE\n";
    }
}

sub exp_char { "E" }


package Fortran::Format::Edit::G;

our @ISA = "Fortran::Format::Edit::E";

sub write_once {
    my ($self, $data) = @_;
    return undef unless @$data;
    my $s;  # working string

    my $d = $data->[0]; # just peek to decide who'll handle the formatting

    # shorthand
    my $scale     = $self->writer->scale;
    my $width     = $self->{width};
    my $precision = $self->{precision};
    my $exp_width = $self->{exp_width} || 0;

    # get exponent
    my $spf = sprintf "%.3E", $d;
    my ($exp) = $spf =~ /E(.*)/g;  # maybe floor log10 abs is faster?
    
    if ($exp < -1 or $exp >= $precision) {
        # format as E
        $s = $self->SUPER::write_once($data);
    } else {
        my $right_margin = $exp_width ? $exp_width + 2 : 4;

        $self->{width} -= $right_margin;
        $self->{precision} = $precision - $exp - 1;
        $s  = $self->Fortran::Format::Edit::F::write_once($data);
        $s .= " " x $right_margin;
        $self->{precision} = $precision;
        $self->{width} = $width;
    }
    $s || "*" x $width;
}

package Fortran::Format::Edit::L;

our @ISA = "Fortran::Format::Node";

sub parse {
    my ($self, $tokenizer) = @_;
    $self->{width} = $tokenizer->get_tok('^\d+$') 
        or die "expected \\d after F\n";
}

sub write_once {
    my ($self, $data) = @_;
    return undef unless @$data;
    my $l = shift @$data; 
    sprintf "%$self->{width}s", $l ? 'T' : 'F';
}

sub read_once {
    my ($self) = @_;
    return undef unless $self->writer->want_more;
    my $s = $self->writer->read($self->{width});
    my $b;

    if ($s =~ /^ *\.?[tT]/) {
        $b = 1;
    } elsif ($s =~ /^ *\.?[fF]/) {
        $b = 0;
    } else {
        die "invalid F format '$s'\n";
    }
    $self->writer->put($b);
    1;
}

package Fortran::Format::Edit::X;

our @ISA = "Fortran::Format::Node";

sub write_once {
    my ($self, $data) = @_;
    $self->writer->position( relative => 1 );
    "";
}

package Fortran::Format::Edit::SLASH;

our @ISA = "Fortran::Format::Node";

sub write_once {
    my ($self, $data) = @_;
    $self->writer->end_line;
    "";
}

package Fortran::Format::Edit::COLON;

our @ISA = "Fortran::Format::Node";

sub write_once {
    my ($self, $data) = @_;
    return undef unless @$data;
    "";
}


package Fortran::Format::Edit::A;

our @ISA = "Fortran::Format::Node";

sub parse {
    my ($self, $tokenizer) = @_;
    $self->{width} = $tokenizer->get_tok('^\d+$');
}

sub write_once {
    my ($self, $data) = @_;
    return undef unless @$data;
    my $datum = shift @$data; 
    my $s;
    if (defined $self->{width}) {
        if (length $datum > $self->{width}) { # truncate
            $s = substr $datum, 0, $self->{width};
        } else { # justify
            $s = sprintf "%$self->{width}s", $datum;
        }
    } else { # use as is
        $s = $datum; 
    }
    $s;
}

sub read_once {
    my ($self) = @_;
    return undef unless $self->writer->want_more;
    my $s = $self->writer->read($self->{width});
    $self->writer->put($s);
    1;
}


package Fortran::Format::Edit::S;

our @ISA = "Fortran::Format::Node";

sub parse {
    my ($self, $tokenizer) = @_;
    $self->{plus} = ''; # default is no plus
    if (my $tok = $tokenizer->get_tok('^[SP]$')) {
        $self->{plus} = '+' if $tok eq 'P';
    }
}

sub write_once {
    my ($self) = @_;
    $self->writer->plus($self->{plus});
    '';
}

package Fortran::Format::Edit::B;

our @ISA = "Fortran::Format::Node";

sub parse {
    my ($self, $tokenizer) = @_;
    my $tok = $tokenizer->get_tok('^[NZ]$')
        or die "expected [NZ] after B\n";
    $self->{bz} = $tok eq 'Z' ? 1 : 0;
}

sub read_once {
    my ($self) = @_;
    $self->writer->bz($self->{bz});
    1;
}

package Fortran::Format::Edit::P;

our @ISA = "Fortran::Format::Node";

sub write_once {
    my ($self) = @_;
    $self->writer->scale($self->{scale});
    '';
}

sub read_once {
    write_once(@_);
}

package Fortran::Format::Edit::T;

our @ISA = "Fortran::Format::Node";

sub parse {
    my ($self, $tokenizer) = @_;
    if ($tokenizer->get_tok('^R$')) {
        my $tok = $tokenizer->get_tok('^\d+$') 
            or die "expected \\d after TR\n";
        $self->{delta} = $tok;
    } elsif ($tokenizer->get_tok('^L$')) {
        my $tok = $tokenizer->get_tok('^\d+$') 
            or die "expected \\d after TL\n";
        $self->{delta} = -$tok;
    } elsif (my $tok = $tokenizer->get_tok('^\d+$')) {
        $self->{position} = $tok;
    } else {
        die "expected \\d after T\n";
    }
}

sub write_once {
    my ($self) = @_;
    if ($self->{position}) { # absolute position (T)
        $self->writer->position( absolute => $self->{position} - 1 ); # Fortran is 1-based
    } elsif ($self->{delta}) { # relative position (TR, TL)
        $self->writer->position( relative => $self->{delta} );
    }
    '';
}


package Fortran::Format::List;

our @ISA = "Fortran::Format::Node";

sub nodes {
    my ($self) = @_;
    @{$self->{nodes}}
}

sub parse {
    my ($self, $tokenizer) = @_;
    $self->{nodes} = my $nodes = [];
    my $repeat = 1;

    while (defined (my $tok = $tokenizer->get_tok)) {
        if ($tok =~ /^[+-]?\d+$/) {
            # should check that next token is repeatable and $tok > 0
            if ($tokenizer->get_tok('P')) {  # scale factor
                push @$nodes, Fortran::Format::Edit::P->build($tokenizer,
                    writer => $self->writer, scale => $tok );
            } elsif ($tokenizer->peek_tok('^[IFEDGLAX(]$')) {
                if ($tok =~ /^[+-]/ or $tok == 0) {
                    die "repeat count should be unsigned and non-zero\n";
                } else {
                    $repeat = $tok;
                }
            } else {
                die "number not followed by repeatable token\n";
            }
        } elsif ($tok eq '(') {
            push @$nodes, $self->{last_list} = Fortran::Format::List->build(
                $tokenizer, 
                repeat => $repeat, 
                writer => $self->writer
            );
            $repeat = 1;
        } elsif ($tok eq ')') {
            return; # end of list
        } elsif ($tok =~ /^'/) {
            push @$nodes, Fortran::Format::Edit::Quote->build($tokenizer,
                writer => $self->writer);
        } elsif ($tok =~ /^[IFEDGLAX]$/i) { # repeatable tokens
            # NOTE: X is technically not a repeatable token; the
            # "repeat" count is suposedly mandatory, but at least g77, ifc, 
            # and pgf77 don't really care (and neither do most programmers)
            push @$nodes, "Fortran::Format::Edit::$tok"->build(
                $tokenizer,
                writer => $self->writer,
                repeat => $repeat, 
            );
            $repeat = 1;
        } elsif ($tok =~ /^([STB]|SLASH|COLON)$/) { # non-repeatable tokens
            push @$nodes, "Fortran::Format::Edit::$tok"->build(
                $tokenizer,
                writer => $self->writer
            );
        } else {
            die "invalid or unimplemented token: $tok\n";
        }
    }
}

sub write_once {
    my ($self, $data, $start) = @_;

    my $started;
    for my $node ($self->nodes) {
        next if $start and !$started and $node != $start;
        $started = 1;

        my $ret = $node->write($data);
        return undef unless defined $ret; # ran out of data ?
        if (length $ret) {
            $self->{writer}->write($ret);
        }
    }
    ''; # this function does not produce new text
}

sub read_once {
    my ($self, $start) = @_;

    my $started;
    for my $node ($self->nodes) {
        next if $start and !$started and $node != $start;
        $started = 1;

        my $ret = $node->read;
        return undef unless defined $ret;
    }
    1;
}

package Fortran::Format::RootList;

our @ISA = "Fortran::Format::List";

sub write {
    my ($self, $data) = @_;
    if ($self->writer->reuse() and $self->{last_list}) {
        $self->SUPER::write($data, $self->{last_list});
    } else {
        $self->SUPER::write($data);
    }
    $self->writer->reuse(1);
    ''; # this function does not produce new text
}

sub read {
    my ($self) = @_;
    if ($self->writer->reuse() and $self->{last_list}) {
        $self->SUPER::read($self->{last_list});
    } else {
        $self->SUPER::read;
    }
    $self->writer->reuse(1);
    ''; # this function does not produce new text
}





1;


=back

=head1 VERSION

0.90

=head1 SEE ALSO

The Fortran format specification: 
L<http://www.fortran.com/fortran/F77_std/rjcnf0001-sh-13.html>

=head1 AUTHOR

Ivan Tubert-Brohman E<lt>itub@cpan.orgE<gt>

=head1 COPYRIGHT

Copyright (c) 2005 Ivan Tubert-Brohman. All rights reserved. This program is
free software; you can redistribute it and/or modify it under the same terms
as Perl itself.

=cut

