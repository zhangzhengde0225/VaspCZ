#  Copyright (c) 1996, 1997 by Steffen Beyer. All rights reserved.
#  Copyright (c) 1999 by Rodolphe Ortalo. All rights reserved.
#  Copyright (c) 2001,2002 by Jonathan Leto. All rights reserved.
#  This package is free software; you can redistribute it and/or
#  modify it under the same terms as Perl itself.

package Math::MatrixReal;

use strict;
use Carp;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);
require Exporter;

@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(min max);
%EXPORT_TAGS = (all => [@EXPORT_OK]);
$VERSION = '2.01';

use overload
     'neg' => '_negate',
       '~' => '_transpose',
    'bool' => '_boolean',
       '!' => '_not_boolean',
      '""' => '_stringify',
     'abs' => '_norm',
       '+' => '_add',
       '-' => '_subtract',
       '*' => '_multiply',
      '**' => '_exponent',
      '+=' => '_assign_add',
      '-=' => '_assign_subtract',
      '*=' => '_assign_multiply',
     '**=' => '_assign_exponent',
      '==' => '_equal',
      '!=' => '_not_equal',
       '<' => '_less_than',
      '<=' => '_less_than_or_equal',
       '>' => '_greater_than',
      '>=' => '_greater_than_or_equal',
      'eq' => '_equal',
      'ne' => '_not_equal',
      'lt' => '_less_than',
      'le' => '_less_than_or_equal',
      'gt' => '_greater_than',
      'ge' => '_greater_than_or_equal',
       '=' => '_clone',
     'exp' => '_exp',
     'sin' => '_sin',
     'cos' => '_cos',
'fallback' =>   undef;

sub new
{
    croak "Usage: \$new_matrix = Math::MatrixReal->new(\$rows,\$columns);"
      if (@_ != 3);

    my ($proto,$rows,$cols) =  @_;
    my $class = ref($proto) || $proto || 'Math::MatrixReal';
    my($i,$j,$this);

    croak "Math::MatrixReal::new(): number of rows must be integer > 0"
      unless ($rows > 0 and  $rows == int($rows) );

    croak "Math::MatrixReal::new(): number of columns must be integer > 0"
      unless ($cols > 0 and $cols == int($cols) );

    $this = [ [ ], $rows, $cols ];

    # Creates first empty row
    my $empty = [ ];
    $#$empty = $cols - 1; # Lengthens the array
    for (my $j = 0; $j < $cols; $j++)
    {
        $empty->[$j] = 0.0;
    }
    $this->[0][0] = $empty;
    # Creates other rows (by copying)
    for (my $i = 1; $i < $rows; $i++)
    {
        my $arow = [ ];
        @$arow = @$empty;
        $this->[0][$i] = $arow;
    }
    bless($this, $class);
    return($this);
}
sub new_diag {
    croak "Usage: \$new_matrix = Math::MatrixReal->new_diag( [ 1, 2, 3] );"
        unless (@_ == 2 );
    my ($proto,$diag) = @_;
    my $class = ref($proto) || $proto || 'Math::MatrixReal';
    my $matrix;
    my $n = scalar(@$diag);

    croak "Math::MatrixReal::new_diag(): Third argument must be an arrayref" 
        unless (ref($diag) eq "ARRAY");

    $matrix = Math::MatrixReal->new($n,$n);
    $matrix = $matrix->each_diag(  sub { shift @$diag } );
    return $matrix;
}

sub new_from_string
{
    croak "Usage: \$new_matrix = Math::MatrixReal->new_from_string(\$string);"
      if (@_ != 2);

    my ($proto,$string)  = @_;
    my $class  = ref($proto) || $proto || 'Math::MatrixReal';
    my($line,$values);
    my($rows,$cols);
    my($row,$col);
    my($warn,$this);

    $warn = $rows = $cols = 0;

    $values = [ ];
    while ($string =~ m!^\s*
  \[ \s+ ( (?: [+-]? \d+ (?: \. \d* )? (?: E [+-]? \d+ )? \s+ )+ ) \] \s*? \n
    !ix)
    {
        $line = $1;
        $string = $';
        $values->[$rows] = [ ];
        @{$values->[$rows]} = split(' ', $line);
        $col = @{$values->[$rows]};
        if ($col != $cols)
        {
            unless ($cols == 0) { $warn = 1; }
            if ($col > $cols) { $cols = $col; }
        }
        $rows++;
    }
    if ($string !~ m/^\s*$/)
    {
        #croak
        print "Math::MatrixReal::new_from_string(): syntax error in input string";
        print "String is\n$string\n---\n";
        croak;
    }
    if ($rows == 0)
    {
        croak "Math::MatrixReal::new_from_string(): empty input string";
    }
    if ($warn)
    {
        warn "Math::MatrixReal::new_from_string(): missing elements will be set to zero!\n";
    }
    $this = Math::MatrixReal::new($class,$rows,$cols);
    for ( $row = 0; $row < $rows; $row++ )
    {
        for ( $col = 0; $col < @{$values->[$row]}; $col++ )
        {
            $this->[0][$row][$col] = $values->[$row][$col];
        }
    }
    return($this);
}
# from Math::MatrixReal::Ext1 (msouth@fulcrum.org)
sub new_from_cols {
    my $this = shift;
    my $extra_args = ( @_ > 1 && ref($_[-1]) eq 'HASH' ) ? pop : {};
    $extra_args->{_type} = 'column';

    $this->_new_from_rows_or_cols(@_, $extra_args );
}
# from Math::MatrixReal::Ext1 (msouth@fulcrum.org)
sub new_from_columns {
    my $this = shift;
    $this->new_from_cols(@_);
}
# from Math::MatrixReal::Ext1 (msouth@fulcrum.org)
sub new_from_rows {
    my $this = shift;
    my $extra_args = ( @_ > 1 && ref($_[-1]) eq 'HASH' ) ? pop : {};
    $extra_args->{_type} = 'row';

    $this->_new_from_rows_or_cols(@_, $extra_args );
}

# from Math::MatrixReal::Ext1 (msouth@fulcrum.org)
sub _new_from_rows_or_cols {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $ref_to_vectors = shift;

    # these additional args are internal at the moment, 
    # but in the future the user could pass e.g. {pad=>1} to
    # request padding
    my $args = pop;
    my $vector_type = $args->{_type};
    die "Internal ".__PACKAGE__." error" unless $vector_type =~ /^(row|column)$/;

    # step back one frame because this private method is 
    # not how the user called it
    my $caller_subname = (caller(1))[3];

    # note--this die() could be inconvenient if someone had something
    # really fancy that knew how to be dereffed as an array 
    # (can you do that with a tied scalar?), but I'm not putting
    # the rest of the world through an eval--they can just
    # deref and pass a reference themselves.  If that ever happens
    # we can add an arg to skip this check
    croak "$caller_subname: need a reference to an array of ${vector_type}s" unless ref($ref_to_vectors) eq 'ARRAY';
    my @vectors = @{$ref_to_vectors};

    my $matrix;

    my $other_type = {row=>'column', column=>'row'}->{$vector_type};

    my %matrix_dim = (
        $vector_type => scalar( @vectors ), 
        $other_type  => 0,  # we will correct this in a bit
    );

    # row and column indices are one based
    my $current_vector_count = 1; 
    foreach my $current_vector (@vectors) {
        # dimension is one-based, so we're
        # starting with one here and incrementing
        # as we go.  The other dimension is fixed (for now, until
        # we add the 'pad' option), and gets set later
        my $ref = ref( $current_vector ) ;

        if ( $ref eq '' ) {
            # we hope this is a properly formatted Math::MatrixReal string,
            # but if not we just let the Math::MatrixReal die() do it's
            # thing
            $current_vector = $class->new_from_string( $current_vector );
        }
        elsif ( $ref eq 'ARRAY' ) {
            my @array = @$current_vector;
            croak "$caller_subname: one $vector_type you gave me was a ref to an array with no elements" unless @array ;
            # we need to create the right kind of string based on whether
            # they said they were sending us rows or columns:
            if ($vector_type eq 'row') {
                $current_vector = $class->new_from_string( '[ '. join( " ", @array) ." ]\n" );
            }
            else {
                $current_vector = $class->new_from_string( '[ '. join( " ]\n[ ", @array) ." ]\n" );
            }
        }
        elsif ( $ref ne 'HASH' and $current_vector->isa('Math::MatrixReal') ) {
            # it's already a Math::MatrixReal something.
            # we don't need to do anything, it will all
            # work out
        }
        else {
            # we have no idea, error time!
            croak "$caller_subname: I only know how to deal with array refs, strings, and things that inherit from Math::MatrixReal\n";
        }

        # starting now we know $current_vector isa Math::MatrixReal thingy
        my @vector_dims = $current_vector->dim;

        #die unless the appropriate dimension is 1
        croak "$caller_subname: I don't accept $other_type vectors"
            unless ($vector_dims[ $vector_type eq 'row' ? 0 : 1 ] == 1) ;

        # the other dimension is the length of our vector
        my $length =  $vector_dims[ $vector_type eq 'row' ? 1 : 0 ];

        # set the "other" dimension to the length of this
        # vector the first time through
        $matrix_dim{$other_type} ||= $length;

        # die unless length of this vector matches the first length
        croak "$caller_subname: one $vector_type has [$length] elements and another one had [$matrix_dim{$other_type}]--all of the ${vector_type}s passed in must have the same dimension"
              unless ($length == $matrix_dim{$other_type}) ;

        # create the matrix the first time through
        $matrix ||= $class->new($matrix_dim{row}, $matrix_dim{column});

        # step along the vector assigning the value of each element
        # to the correct place in the matrix we're building
        foreach my $element_index ( 1..$length ){
            # args for vector assignment:
            # initialize both to one and reset the correct
            # one below
            my ($v_r, $v_c) = (1,1);

            # args for matrix assignment
            my ($row_index, $column_index, $value);

            if ($vector_type eq 'row') {
                $row_index           = $current_vector_count;
                $v_c = $column_index = $element_index;
            }
            else {
                $v_r = $row_index    = $element_index;
                $column_index        = $current_vector_count;
            }
            $value = $current_vector->element($v_r, $v_c);
            $matrix->assign($row_index, $column_index, $value);
        }
        $current_vector_count ++ ;
    }
    return $matrix;
}

sub shadow
{
    croak "Usage: \$new_matrix = \$some_matrix->shadow();"
      if (@_ != 1);

    my($matrix) = @_;
    my($temp);

    $temp = $matrix->new($matrix->[1],$matrix->[2]);
    return($temp);
}


sub copy
{
    croak "Usage: \$matrix1->copy(\$matrix2);"
      if (@_ != 2);

    my($matrix1,$matrix2) = @_;
    my($rows1,$cols1) = ($matrix1->[1],$matrix1->[2]);
    my($rows2,$cols2) = ($matrix2->[1],$matrix2->[2]);
    my($i,$j);

    croak "Math::MatrixReal::copy(): matrix size mismatch"
      unless (($rows1 == $rows2) && ($cols1 == $cols2));

    for ( $i = 0; $i < $rows1; $i++ )
    {
        my $r1 = []; # New array ref
        my $r2 = $matrix2->[0][$i];
        @$r1 = @$r2; # Copy whole array directly
        $matrix1->[0][$i] = $r1;
    }
    if (defined $matrix2->[3]) # is an LR decomposition matrix!
    {
        $matrix1->[3] = $matrix2->[3]; # $sign
        $matrix1->[4] = $matrix2->[4]; # $perm_row
        $matrix1->[5] = $matrix2->[5]; # $perm_col
    }
}

sub clone
{
    croak "Usage: \$twin_matrix = \$some_matrix->clone();"
      if (@_ != 1);

    my($matrix) = @_;
    my($temp);

    $temp = $matrix->new($matrix->[1],$matrix->[2]);
    $temp->copy($matrix);
    return($temp);
}

## trace() : return the sum of the diagonal elements
sub trace {
    croak "Usage: \$trace = \$matrix->trace();" if (@_ != 1);

    my $matrix = shift;
    my($rows,$cols) = ($matrix->[1],$matrix->[2]);
    my $trace = 0;
    my($j);

    croak "Math::MatrixReal::trace(): matrix is not quadratic"
      unless ($rows == $cols);


    for ( $j = 0; $j < $cols; $j++ )
    {
        $trace += $matrix->[0][$j][$j];
    }
    return($trace);
}

## return the minor corresponding to $r and $c
## eliminate row $r and col $c , and return the $r-1 by $c-1 matrix
sub minor {
    croak "Usage: \$minor = \$matrix->minor(\$r,\$c);" unless (@_ == 3);
    my ($matrix,$r,$c) = @_;
    my ($rows,$cols) = $matrix->dim();

    croak "Math::MatrixReal::minor(): \$matrix must be at least 2x2"
        unless ($rows > 1 and $cols > 1);
    croak "Math::MatrixReal::minor(): $r and $c must be positive" 
        unless ($r > 0 and $c > 0 );
    croak "Math::MatrixReal::minor(): matrix has no $r,$c element" 
        unless ($r <= $rows and $c <= $cols );

    my $minor = new Math::MatrixReal($rows-1,$cols-1);
    my ($i,$j) = (0,0);

    ## assign() might have been easier, but this should be faster
    ## the element can be in any of 4 regions compared to the eliminated
    ## row and col:
    ## above and to the left, above and to the right
    ## below and to the left, below and to the right

    for(; $i < $rows; $i++){
        for(;$j < $rows; $j++ ){
            if( $i >= $r && $j >= $c ){
                $minor->[0][$i-1][$j-1] = $matrix->[0][$i][$j];
            } elsif ( $i >= $r && $j < $c ){
                $minor->[0][$i-1][$j] = $matrix->[0][$i][$j];
            } elsif ( $i < $r && $j < $c ){
                $minor->[0][$i][$j] = $matrix->[0][$i][$j];
            } elsif ( $i < $r && $j >= $c ){
                $minor->[0][$i][$j-1] = $matrix->[0][$i][$j];
            } else {
                croak "Very bad things";
            }
        }
        $j = 0;
    }
    return ($minor);
}
sub swap_col {
    croak "Usage: \$matrix->swap_col(\$col1,\$col2); " unless (@_ == 3);
    my ($matrix,$col1,$col2) = @_;
    my ($rows,$cols) = $matrix->dim();
    my (@temp);

    croak "Math::MatrixReal::swap_col(): col index is not valid"
        unless ( $col1 <= $cols && $col2 <= $cols &&
             $col1 == int($col1) &&
             $col2 == int($col2) );
    $col1--;$col2--;
    for(my $i=0;$i < $rows;$i++){
        $temp[$i] = $matrix->[0][$i][$col1];
        $matrix->[0][$i][$col1] = $matrix->[0][$i][$col2];
        $matrix->[0][$i][$col2] =  $temp[$i];
    }
}
sub swap_row {
    croak "Usage: \$matrix->swap_row(\$row1,\$row2); " unless (@_ == 3);
    my ($matrix,$row1,$row2) = @_;
    my ($rows,$cols) = $matrix->dim();
    my (@temp);

    croak "Math::MatrixReal::swap_row(): row index is not valid"
        unless ( $row1 <= $rows && $row2 <= $rows && 
             $row1 == int($row1) && 
             $row2 == int($row2) ); 
    $row1--;$row2--;
    for(my $j=0;$j < $cols;$j++){
        $temp[$j] = $matrix->[0][$row1][$j];
        $matrix->[0][$row1][$j] = $matrix->[0][$row2][$j];
        $matrix->[0][$row2][$j] =  $temp[$j];
    }
}
# TODO: docs
# no worky
sub assign_row {
    croak "Usage: \$matrix->assign_row(\$row,\$row_vec);"  unless (@_ == 3);
    my ($matrix,$row,$row_vec) = @_;
    my ($rows1,$cols1) = $matrix->dim();
    my ($rows2,$cols2) = $row_vec->dim();
    
    croak "Math::MatrixReal::assign_row(): row mismatch" if ($rows1 != $rows2);
    croak "Math::MatrixReal::assign_row(): not a row vector" unless( $cols2 == 1);

    @{$matrix->[0][--$row]} = @{$row_vec->[0][0]};

}
# returns the number of zeroes in a row
sub _count_zeroes_row {
    my ($matrix) = @_;
    my ($rows,$cols) = $matrix->dim();
    my $count = 0;
    croak "_count_zeroes_row(): only 1 row, buddy" unless ($rows == 1);

    for(my $i=0;$i < $cols;$i++){
        $count++ unless $matrix->[0][0][$i];
    }
    return $count;
}
## divide a row by it's largest abs() element
sub _normalize_row {
    my ($matrix) = @_;
    my ($rows,$cols) = $matrix->dim();
    my $new_row = Math::MatrixReal->new(1,$cols);

    my $big = abs($matrix->[0][0][0]);
    for(my $j=0;$j < $cols; $j++ ){
        $big = $big < abs($matrix->[0][0][$j])
            ? abs($matrix->[0][0][$j]) : $big;
    }
    next unless $big;
    # now $big is biggest element in row
    for(my $j = 0;$j < $cols; $j++ ){
        $new_row->[0][0][$j]  = $matrix->[0][0][$j] / $big;
    }
    return $new_row;

}
#### doesn't work
=crap
sub row_echelon {
    my ($matrix) = @_;
    my ($rows,$cols) = $matrix->dim();
    my $big;
    my $tmprow;

    for(my $i = 0;$i < $rows; $i++ ){
        $big = abs($matrix->[0][$i][0]);
        for(my $j=0;$j < $cols; $j++ ){
            $big = $big < abs($matrix->[0][$i][$j]) 
                ? abs($matrix->[0][$i][$j]) : $big;
        }
        next unless $big;
        # now $big is biggest element in row
                for(my $j = 0;$j < $cols; $j++ ){
                        $matrix->[0][$i][$j]  /= $big;
                }
        # now all elements are between [-1,1] and 
        # dependence can be easily checked
    }

    for(my $i = 0;$i < $rows; $i++ ){
        for(my $k = 0; $k < $rows; $k++ ){
            next if ($i == $k);
            print "i,k = $i,$k\n";
            $tmprow = $matrix->row($i+1) + $matrix->row($k+1);
            print "tmprow= $tmprow\n";
            print "tmprow has " . $tmprow->_count_zeroes_row . " zeroes.\n";
            print "i+1 has " . $matrix->row($i+1)->_count_zeroes_row . " zeroes\n";
            print "k+1 has " . $matrix->row($k+1)->_count_zeroes_row . " zeroes\n";

            $tmprow = $tmprow->_normalize_row;
            if( $tmprow->norm_sum == 0 ){
                for(my $j=0;$j<$cols;$j++){
                    $matrix->[0][$k][$j] = 0;
                }
            } elsif ( $tmprow->_count_zeroes_row > $matrix->row($i+1)->_count_zeroes_row ){
                print "tmprow has more zeroes, replacing row " . ($i+1) . "\n";
                ## assign row
                $matrix->assign_row($i+1,$tmprow);
            } elsif ( $tmprow->_count_zeroes_row > $matrix->row($k+1)->_count_zeroes_row ){
                print "tmprow has more zeroes, replacing row " . ($k+1) . "\n";
                $matrix->assign_row($k+1,$tmprow);
            }
                    
        }
    }
    return $matrix;
}
=cut

sub cofactor {
    my ($matrix) = @_;
    my ($rows,$cols) = $matrix->dim();

    croak "Math::MatrixReal::cofactor(): Matrix is not quadratic" 
        unless ($rows == $cols);

    # black magic ahead
    my $cofactor = $matrix->each( sub { my($v,$i,$j) = @_;
        ($i+$j) % 2 == 0 ? $matrix->minor($i,$j)->det()
        : -1*$matrix->minor($i,$j)->det(); 
        } );
    return ($cofactor);
}

sub adjoint {
    my ($matrix) = @_;
    return ~($matrix->cofactor);
}
sub row
{
    croak "Usage: \$row_vector = \$matrix->row(\$row);"
      if (@_ != 2);

    my($matrix,$row) = @_;
    my($rows,$cols) = ($matrix->[1],$matrix->[2]);
    my($temp);
    my($j);

    croak "Math::MatrixReal::row(): row index out of range"
      if (($row < 1) || ($row > $rows));

    $row--;
    $temp = $matrix->new(1,$cols);
    for ( $j = 0; $j < $cols; $j++ )
    {
        $temp->[0][0][$j] = $matrix->[0][$row][$j];
    }
    return($temp);
}

sub column
{
    croak "Usage: \$column_vector = \$matrix->column(\$column);"
      if (@_ != 2);

    my($matrix,$col) = @_;
    my($rows,$cols) = ($matrix->[1],$matrix->[2]);
    my($temp);
    my($i);

    croak "Math::MatrixReal::column(): column index out of range"
      if (($col < 1) || ($col > $cols));

    $col--;
    $temp = $matrix->new($rows,1);
    for ( $i = 0; $i < $rows; $i++ )
    {
        $temp->[0][$i][0] = $matrix->[0][$i][$col];
    }
    return($temp);
}

sub _undo_LR
{
    croak "Usage: \$matrix->_undo_LR();"
      if (@_ != 1);

    my($this) = @_;

    undef $this->[3];
    undef $this->[4];
    undef $this->[5];
}

sub zero
{
    croak "Usage: \$matrix->zero();"
      if (@_ != 1);

    my($this) = @_;
    my($rows,$cols) = ($this->[1],$this->[2]);
    my($i,$j);

    $this->_undo_LR();

    # Zero first row
    for (my $j = 0; $j < $cols; $j++ )
    {
        $this->[0][0][$j] = 0.0;
    }
    # Then propagate to other rows
    for (my $i = 0; $i < $rows; $i++)
    {
        @{$this->[0][$i]} = @{$this->[0][0]};
    }
}

sub one
{
    croak "Usage: \$matrix->one();"
      if (@_ != 1);

    my($this) = @_;
    my($rows,$cols) = ($this->[1],$this->[2]);
    my($i,$j);

# No need for this: done by the 'zero()'
#    $this->_undo_LR();
    $this->zero(); # We rely on zero() efficiency
    for (my $i = 0; $i < $rows; $i++ )
    {
        $this->[0][$i][$i] = 1.0;
    }
}

sub assign
{
    croak "Usage: \$matrix->assign(\$row,\$column,\$value);"
      if (@_ != 4);

    my($this,$row,$col,$value) = @_;
    my($rows,$cols) = ($this->[1],$this->[2]);

    croak "Math::MatrixReal::assign(): row index out of range"
      if (($row < 1) || ($row > $rows));

    croak "Math::MatrixReal::assign(): column index out of range"
      if (($col < 1) || ($col > $cols));

    $this->_undo_LR();

    $this->[0][--$row][--$col] = $value;
}

sub element
{
    croak "Usage: \$value = \$matrix->element(\$row,\$column);"
      if (@_ != 3);

    my($this,$row,$col) = @_;
    my($rows,$cols) = ($this->[1],$this->[2]);

    croak "Math::MatrixReal::element(): row index out of range"
      if (($row < 1) || ($row > $rows));

    croak "Math::MatrixReal::element(): column index out of range"
      if (($col < 1) || ($col > $cols));

    return( $this->[0][--$row][--$col] );
}

sub dim  #  returns dimensions of a matrix
{
    croak "Usage: (\$rows,\$columns) = \$matrix->dim();"
      if (@_ != 1);

    my($matrix) = @_;

    return( $matrix->[1], $matrix->[2] );
}

sub norm_one  #  maximum of sums of each column
{
    croak "Usage: \$norm_one = \$matrix->norm_one();"
      if (@_ != 1);

    my($this) = @_;
    my($rows,$cols) = ($this->[1],$this->[2]);

    my $max = 0.0;
    for (my $j = 0; $j < $cols; $j++)
    {
        my $sum = 0.0;
        for (my $i = 0; $i < $rows; $i++)
        {
            $sum += abs( $this->[0][$i][$j] );
        }
        $max = $sum if ($sum > $max);
    }
    return($max);
}
## sum of absolute value of every element
sub norm_sum {
    croak "Usage: \$norm_sum = \$matrix->norm_sum();"
        unless (@_ == 1);
    my ($matrix) = @_;
    my $norm = 0;
    $matrix->each( sub { $norm+=abs(shift); } );
    return $norm;
}
sub norm_frobenius {
    my ($m) = @_;
    my ($r,$c) = $m->dim;
    my $s=0;

    $m->each( sub { $s+=abs(shift)**2 } );
    return sqrt($s);
}

# Vector Norm 
sub norm_p {
    my ($v,$p) = @_;
    # sanity check on $p
    croak "Math::MatrixReal:norm_p: argument must be a row or column vector"
        unless ( $v->is_row_vector || $v->is_col_vector );
    croak "Math::MatrixReal::norm_p: $p must be >= 1" 
        unless ($p =~ m/Inf(inity)?/i || $p >= 1);

    if( $p =~ m/^(Inf|Infinity)$/i ){
        my $max = $v->element(1,1);
        $v->each ( sub { my $x=abs(shift); $max = $x if( $x > $max ); } );
        return $max;
    }

    my $s=0;
    $v->each( sub { $s+= (abs(shift))**$p; } );
    return $s ** (1/$p);

}
sub norm_max  #  maximum of sums of each row
{
    croak "Usage: \$norm_max = \$matrix->norm_max();"
      if (@_ != 1);

    my($this) = @_;
    my($rows,$cols) = ($this->[1],$this->[2]);

    my $max = 0.0;
    for (my $i = 0; $i < $rows; $i++)
    {
        my $sum = 0.0;
        for (my $j = 0; $j < $cols; $j++)
        {
            $sum += abs( $this->[0][$i][$j] );
        }
        $max = $sum if ($sum > $max);
    }
    return($max);
}

sub negate
{
    croak "Usage: \$matrix1->negate(\$matrix2);"
      if (@_ != 2);

    my($matrix1,$matrix2) = @_;
    my($rows1,$cols1) = ($matrix1->[1],$matrix1->[2]);
    my($rows2,$cols2) = ($matrix2->[1],$matrix2->[2]);

    croak "Math::MatrixReal::negate(): matrix size mismatch"
      unless (($rows1 == $rows2) && ($cols1 == $cols2));

    $matrix1->_undo_LR();

    for (my $i = 0; $i < $rows1; $i++ )
    {
        for (my $j = 0; $j < $cols1; $j++ )
        {
            $matrix1->[0][$i][$j] = -($matrix2->[0][$i][$j]);
        }
    }
}
## each(): evaluate a coderef on each element and return a new matrix
## of said 
sub each {
    croak "Usage: \$new_matrix = \$matrix->each( \&sub );" unless (@_ == 2 );
    my($matrix,$function) = @_;
    my($rows,$cols) = ($matrix->[1],$matrix->[2]);
    my($new_matrix) = $matrix->clone();

    croak "Math::MatrixReal::each(): argument is not a sub reference" unless ref($function);
    $new_matrix->_undo_LR();

    for (my $i = 0; $i < $rows; $i++ ) {
        for (my $j = 0; $j < $cols; $j++ ) {
            no strict 'refs';
            # $i,$j are 1-based as of 1.7
            $new_matrix->[0][$i][$j] = &{ $function }($matrix->[0][$i][$j],$i+1,$j+1) ;
        }
    }
    return ($new_matrix);
}
## each_diag(): same as each() but only diag elements
sub each_diag { 
    croak "Usage: \$new_matrix = \$matrix->each_diag( \&sub );" unless (@_ == 2 );
    my($matrix,$function) = @_;
    my($rows,$cols) = ($matrix->[1],$matrix->[2]);
    my($new_matrix) = $matrix->clone();

    croak "Math::MatrixReal::each(): argument is not a sub reference" unless ref($function);
    croak "Matrix is not quadratic" unless ($rows == $cols);

    $new_matrix->_undo_LR();

    for (my $i = 0; $i < $rows; $i++ ) {
        for (my $j = 0; $j < $cols; $j++ ) {
            next unless ($i == $j);
            no strict 'refs';
            # $i,$j are 1-based as of 1.7
            $new_matrix->[0][$i][$j] = &{ $function }($matrix->[0][$i][$j],$i+1,$j+1) ;
        }
    }
    return ($new_matrix);
}

## Make computing the inverse more user friendly
sub inverse {
    croak "Usage: \$inverse = \$matrix->inverse();" unless (@_ == 1);
    my ($matrix) = @_;
    return $matrix->decompose_LR->invert_LR;
}

sub transpose {

    croak "Usage: \$matrix1->transpose(\$matrix2);"
      if (@_ != 2);

    my($matrix1,$matrix2) = @_;
    my($rows1,$cols1) = ($matrix1->[1],$matrix1->[2]);
    my($rows2,$cols2) = ($matrix2->[1],$matrix2->[2]);

    croak "Math::MatrixReal::transpose(): matrix size mismatch"
      unless (($rows1 == $cols2) && ($cols1 == $rows2));

    $matrix1->_undo_LR();

    if ($rows1 == $cols1)
    {
        # more complicated to make in-place possible!

        for (my $i = 0; $i < $rows1; $i++)
        {
            for (my $j = ($i + 1); $j < $cols1; $j++)
            {
                my $swap              = $matrix2->[0][$i][$j];
                $matrix1->[0][$i][$j] = $matrix2->[0][$j][$i];
                $matrix1->[0][$j][$i] = $swap;
            }
            $matrix1->[0][$i][$i] = $matrix2->[0][$i][$i];
        }
    }
    else # ($rows1 != $cols1)
    {
        for (my $i = 0; $i < $rows1; $i++)
        {
            for (my $j = 0; $j < $cols1; $j++)
            {
                $matrix1->[0][$i][$j] = $matrix2->[0][$j][$i];
            }
        }
    }
}

sub add
{
    croak "Usage: \$matrix1->add(\$matrix2,\$matrix3);"
      if (@_ != 3);

    my($matrix1,$matrix2,$matrix3) = @_;
    my($rows1,$cols1) = ($matrix1->[1],$matrix1->[2]);
    my($rows2,$cols2) = ($matrix2->[1],$matrix2->[2]);
    my($rows3,$cols3) = ($matrix3->[1],$matrix3->[2]);
    my($i,$j);

    croak "Math::MatrixReal::add(): matrix size mismatch"
      unless (($rows1 == $rows2) && ($rows1 == $rows3) &&
              ($cols1 == $cols2) && ($cols1 == $cols3));

    $matrix1->_undo_LR();

    for ( $i = 0; $i < $rows1; $i++ )
    {
        for ( $j = 0; $j < $cols1; $j++ )
        {
            $matrix1->[0][$i][$j] =
            $matrix2->[0][$i][$j] + $matrix3->[0][$i][$j];
        }
    }
}

sub subtract
{
    croak "Usage: \$matrix1->subtract(\$matrix2,\$matrix3);"
      if (@_ != 3);

    my($matrix1,$matrix2,$matrix3) = @_;
    my($rows1,$cols1) = ($matrix1->[1],$matrix1->[2]);
    my($rows2,$cols2) = ($matrix2->[1],$matrix2->[2]);
    my($rows3,$cols3) = ($matrix3->[1],$matrix3->[2]);
    my($i,$j);

    croak "Math::MatrixReal::subtract(): matrix size mismatch"
      unless (($rows1 == $rows2) && ($rows1 == $rows3) &&
              ($cols1 == $cols2) && ($cols1 == $cols3));

    $matrix1->_undo_LR();

    for ( $i = 0; $i < $rows1; $i++ )
    {
        for ( $j = 0; $j < $cols1; $j++ )
        {
            $matrix1->[0][$i][$j] =
            $matrix2->[0][$i][$j] - $matrix3->[0][$i][$j];
        }
    }
}

sub multiply_scalar
{
    croak "Usage: \$matrix1->multiply_scalar(\$matrix2,\$scalar);"
      if (@_ != 3);

    my($matrix1,$matrix2,$scalar) = @_;
    my($rows1,$cols1) = ($matrix1->[1],$matrix1->[2]);
    my($rows2,$cols2) = ($matrix2->[1],$matrix2->[2]);
    my($i,$j);

    croak "Math::MatrixReal::multiply_scalar(): matrix size mismatch"
      unless (($rows1 == $rows2) && ($cols1 == $cols2));

    $matrix1->_undo_LR();

    for ( $i = 0; $i < $rows1; $i++ )
    {
        for ( $j = 0; $j < $cols1; $j++ )
        {
            $matrix1->[0][$i][$j] = $matrix2->[0][$i][$j] * $scalar;
        }
    }
}

sub multiply
{
    croak "Usage: \$product_matrix = \$matrix1->multiply(\$matrix2);"
      if (@_ != 2);

    my($matrix1,$matrix2) = @_;
    my($rows1,$cols1) = ($matrix1->[1],$matrix1->[2]);
    my($rows2,$cols2) = ($matrix2->[1],$matrix2->[2]);
    my($temp);

    croak "Math::MatrixReal::multiply(): matrix size mismatch"
      unless ($cols1 == $rows2);

    $temp = $matrix1->new($rows1,$cols2);
    for (my $i = 0; $i < $rows1; $i++ )
    {
        for (my $j = 0; $j < $cols2; $j++ )
        {
            my $sum = 0.0;
            for (my $k = 0; $k < $cols1; $k++ )
            {
                $sum += ( $matrix1->[0][$i][$k] * $matrix2->[0][$k][$j] );
            }
            $temp->[0][$i][$j] = $sum;
        }
    }
    return($temp);
}

sub exponent {       
    croak "Usage: \$matrix_exp = \$matrix1->exponent(\$integer);" if(@_ != 2 );
    my($matrix,$argument) = @_;
    my($rows,$cols) = ($matrix->[1],$matrix->[2]);
    my($name) = "'**'"; 
    my($temp) = $matrix->clone();       

    croak "Matrix is not quadratic" unless ($rows == $cols);
    croak "Exponent must be integer" unless ($argument =~ m/^[+-]?\d+$/ );

    return($matrix) if ($argument == 1);

    $temp->_undo_LR();

    # negative exponent is (A^-1)^n
    if( $argument < 0 ){ 
        my $LR = $matrix->decompose_LR();
        my $inverse = $LR->invert_LR();
        unless (defined $inverse){
            carp "Matrix has no inverse";
            return undef;
        }
        $temp = $inverse->clone();
        if( $inverse ){
            return($inverse) if ($argument == -1);
                for( 2 .. abs($argument) ){ 
                    $temp = multiply($inverse,$temp);
                    }
                return($temp);
            } else {   
            # TODO: is this the right behaviour?
                    carp "Cannot compute negative exponent, inverse does not exist";
            return undef;
        }
    # matrix to zero power is identity matrix
    } elsif( $argument == 0 ){
        $temp->one();
        return ($temp);
    }

    # if it is diagonal, just raise diagonal entries to power
    if( $matrix->is_diagonal() ){
        $temp = $temp->each_diag( sub { (shift)**$argument } );
        return ($temp);
    
    } else {
        for( 2 .. $argument ){
            $temp = multiply($matrix,$temp);
        }
        return ($temp);
    }
}

sub min
{
    croak "Usage: \$minimum = Math::MatrixReal::min(\$number1,\$number2);"
      if (@_ != 2);

    return( $_[0] < $_[1] ? $_[0] : $_[1] );
}

sub max
{
    croak "Usage: \$maximum = Math::MatrixReal::max(\$number1,\$number2);"
      if (@_ != 2);

    return( $_[0] > $_[1] ? $_[0] : $_[1] );
}

sub kleene
{
    croak "Usage: \$minimal_cost_matrix = \$cost_matrix->kleene();"
      if (@_ != 1);

    my($matrix) = @_;
    my($rows,$cols) = ($matrix->[1],$matrix->[2]);
    my($i,$j,$k,$n);
    my($temp);

    croak "Math::MatrixReal::kleene(): matrix is not quadratic"
      unless ($rows == $cols);

    $temp = $matrix->new($rows,$cols);
    $temp->copy($matrix);
    $temp->_undo_LR();
    $n = $rows;
    for ( $i = 0; $i < $n; $i++ )
    {
        $temp->[0][$i][$i] = min( $temp->[0][$i][$i] , 0 );
    }
    for ( $k = 0; $k < $n; $k++ )
    {
        for ( $i = 0; $i < $n; $i++ )
        {
            for ( $j = 0; $j < $n; $j++ )
            {
                $temp->[0][$i][$j] = min( $temp->[0][$i][$j] ,
                                        ( $temp->[0][$i][$k] +
                                          $temp->[0][$k][$j] ) );
            }
        }
    }
    return($temp);
}

sub normalize
{
    croak "Usage: (\$norm_matrix,\$norm_vector) = \$matrix->normalize(\$vector);"
      if (@_ != 2);

    my($matrix,$vector) = @_;
    my($rows,$cols) = ($matrix->[1],$matrix->[2]);
    my($norm_matrix,$norm_vector);
    my($max,$val);
    my($i,$j,$n);

    croak "Math::MatrixReal::normalize(): matrix is not quadratic"
      unless ($rows == $cols);

    $n = $rows;

    croak "Math::MatrixReal::normalize(): vector is not a column vector"
      unless ($vector->[2] == 1);

    croak "Math::MatrixReal::normalize(): matrix and vector size mismatch"
      unless ($vector->[1] == $n);

    $norm_matrix = $matrix->new($n,$n);
    $norm_vector = $vector->new($n,1);

    $norm_matrix->copy($matrix);
    $norm_vector->copy($vector);

    $norm_matrix->_undo_LR();

    for ( $i = 0; $i < $n; $i++ )
    {
        $max = abs($norm_vector->[0][$i][0]);
        for ( $j = 0; $j < $n; $j++ )
        {
            $val = abs($norm_matrix->[0][$i][$j]);
            if ($val > $max) { $max = $val; }
        }
        if ($max != 0)
        {
            $norm_vector->[0][$i][0] /= $max;
            for ( $j = 0; $j < $n; $j++ )
            {
                $norm_matrix->[0][$i][$j] /= $max;
            }
        }
    }
    return($norm_matrix,$norm_vector);
}

sub decompose_LR
{
    croak "Usage: \$LR_matrix = \$matrix->decompose_LR();"
      if (@_ != 1);

    my($matrix) = @_;
    my($rows,$cols) = ($matrix->[1],$matrix->[2]);
    my($perm_row,$perm_col);
    my($row,$col,$max);
    my($i,$j,$k,$n);
    my($sign) = 1;
    my($swap);
    my($temp);

    croak "Math::MatrixReal::decompose_LR(): matrix is not quadratic"
      unless ($rows == $cols);

    $temp = $matrix->new($rows,$cols);
    $temp->copy($matrix);
    $n = $rows;
    $perm_row = [ ];
    $perm_col = [ ];
    for ( $i = 0; $i < $n; $i++ )
    {
        $perm_row->[$i] = $i;
        $perm_col->[$i] = $i;
    }
    NONZERO:
    for ( $k = 0; $k < $n; $k++ ) # use Gauss's algorithm:
    {
        # complete pivot-search:

        $max = 0;
        for ( $i = $k; $i < $n; $i++ )
        {
            for ( $j = $k; $j < $n; $j++ )
            {
                if (($swap = abs($temp->[0][$i][$j])) > $max)
                {
                    $max = $swap;
                    $row = $i;
                    $col = $j;
                }
            }
        }
        last NONZERO if ($max == 0); # (all remaining elements are zero)
        if ($k != $row) # swap row $k and row $row:
        {
            $sign = -$sign;
            $swap             = $perm_row->[$k];
            $perm_row->[$k]   = $perm_row->[$row];
            $perm_row->[$row] = $swap;
            for ( $j = 0; $j < $n; $j++ )
            {
                # (must run from 0 since L has to be swapped too!)

                $swap                = $temp->[0][$k][$j];
                $temp->[0][$k][$j]   = $temp->[0][$row][$j];
                $temp->[0][$row][$j] = $swap;
            }
        }
        if ($k != $col) # swap column $k and column $col:
        {
            $sign = -$sign;
            $swap             = $perm_col->[$k];
            $perm_col->[$k]   = $perm_col->[$col];
            $perm_col->[$col] = $swap;
            for ( $i = 0; $i < $n; $i++ )
            {
                $swap                = $temp->[0][$i][$k];
                $temp->[0][$i][$k]   = $temp->[0][$i][$col];
                $temp->[0][$i][$col] = $swap;
            }
        }
        for ( $i = ($k + 1); $i < $n; $i++ )
        {
            # scan the remaining rows, add multiples of row $k to row $i:

            $swap = $temp->[0][$i][$k] / $temp->[0][$k][$k];
            if ($swap != 0)
            {
                # calculate a row of matrix R:

                for ( $j = ($k + 1); $j < $n; $j++ )
                {
                    $temp->[0][$i][$j] -= $temp->[0][$k][$j] * $swap;
                }

                # store matrix L in same matrix as R:

                $temp->[0][$i][$k] = $swap;
            }
        }
    }
    $temp->[3] = $sign;
    $temp->[4] = $perm_row;
    $temp->[5] = $perm_col;
    return($temp);
}

sub solve_LR
{
    croak "Usage: (\$dimension,\$x_vector,\$base_matrix) = \$LR_matrix->solve_LR(\$b_vector);"
      if (@_ != 2);

    my($LR_matrix,$b_vector) = @_;
    my($rows,$cols) = ($LR_matrix->[1],$LR_matrix->[2]);
    my($dimension,$x_vector,$base_matrix);
    my($perm_row,$perm_col);
    my($y_vector,$sum);
    my($i,$j,$k,$n);

    croak "Math::MatrixReal::solve_LR(): not an LR decomposition matrix"
      unless ((defined $LR_matrix->[3]) && ($rows == $cols));

    $n = $rows;

    croak "Math::MatrixReal::solve_LR(): vector is not a column vector"
      unless ($b_vector->[2] == 1);

    croak "Math::MatrixReal::solve_LR(): matrix and vector size mismatch"
      unless ($b_vector->[1] == $n);

    $perm_row = $LR_matrix->[4];
    $perm_col = $LR_matrix->[5];

    $x_vector    =   $b_vector->new($n,1);
    $y_vector    =   $b_vector->new($n,1);
    $base_matrix = $LR_matrix->new($n,$n);

    # calculate "x" so that LRx = b  ==>  calculate Ly = b, Rx = y:

    for ( $i = 0; $i < $n; $i++ ) # calculate $y_vector:
    {
        $sum = $b_vector->[0][($perm_row->[$i])][0];
        for ( $j = 0; $j < $i; $j++ )
        {
            $sum -= $LR_matrix->[0][$i][$j] * $y_vector->[0][$j][0];
        }
        $y_vector->[0][$i][0] = $sum;
    }

    $dimension = 0;
    for ( $i = ($n - 1); $i >= 0; $i-- ) # calculate $x_vector:
    {
        if ($LR_matrix->[0][$i][$i] == 0)
        {
            if ($y_vector->[0][$i][0] != 0)
            {
                return(); # a solution does not exist!
            }
            else
            {
                $dimension++;
                $x_vector->[0][($perm_col->[$i])][0] = 0;
            }
        }
        else
        {
            $sum = $y_vector->[0][$i][0];
            for ( $j = ($i + 1); $j < $n; $j++ )
            {
                $sum -= $LR_matrix->[0][$i][$j] *
                    $x_vector->[0][($perm_col->[$j])][0];
            }
            $x_vector->[0][($perm_col->[$i])][0] =
                $sum / $LR_matrix->[0][$i][$i];
        }
    }
    if ($dimension)
    {
        if ($dimension == $n)
        {
            $base_matrix->one();
        }
        else
        {
            for ( $k = 0; $k < $dimension; $k++ )
            {
                $base_matrix->[0][($perm_col->[($n-$k-1)])][$k] = 1;
                for ( $i = ($n-$dimension-1); $i >= 0; $i-- )
                {
                    $sum = 0;
                    for ( $j = ($i + 1); $j < $n; $j++ )
                    {
                        $sum -= $LR_matrix->[0][$i][$j] *
                            $base_matrix->[0][($perm_col->[$j])][$k];
                    }
                    $base_matrix->[0][($perm_col->[$i])][$k] =
                        $sum / $LR_matrix->[0][$i][$i];
                }
            }
        }
    }
    return( $dimension, $x_vector, $base_matrix );
}

sub invert_LR
{
    croak "Usage: \$inverse_matrix = \$LR_matrix->invert_LR();"
      if (@_ != 1);

    my($matrix) = @_;
    my($rows,$cols) = ($matrix->[1],$matrix->[2]);
    my($inv_matrix,$x_vector,$y_vector);
    my($i,$j,$n);

    croak "Math::MatrixReal::invert_LR(): not an LR decomposition matrix"
      unless ((defined $matrix->[3]) && ($rows == $cols));

    $n = $rows;
    if ($matrix->[0][$n-1][$n-1] != 0)
    {
        $inv_matrix = $matrix->new($n,$n);
        $y_vector   = $matrix->new($n,1);
        for ( $j = 0; $j < $n; $j++ )
        {
            if ($j > 0)
            {
                $y_vector->[0][$j-1][0] = 0;
            }
            $y_vector->[0][$j][0] = 1;
            if (($rows,$x_vector,$cols) = $matrix->solve_LR($y_vector))
            {
                for ( $i = 0; $i < $n; $i++ )
                {
                    $inv_matrix->[0][$i][$j] = $x_vector->[0][$i][0];
                }
            }
            else
            {
                die "Math::MatrixReal::invert_LR(): unexpected error - please inform author!\n";
            }
        }
        return($inv_matrix);
    }
    else { return(); } # matrix is not invertible!
}

sub condition
{
    # 1st matrix MUST be the inverse of 2nd matrix (or vice-versa)
    # for a meaningful result!

    # make this work when given no args

    croak "Usage: \$condition = \$matrix->condition(\$inverse_matrix);"
      if (@_ != 2);

    my($matrix1,$matrix2) = @_;
    my($rows1,$cols1) = ($matrix1->[1],$matrix1->[2]);
    my($rows2,$cols2) = ($matrix2->[1],$matrix2->[2]);

    croak "Math::MatrixReal::condition(): 1st matrix is not quadratic"
      unless ($rows1 == $cols1);

    croak "Math::MatrixReal::condition(): 2nd matrix is not quadratic"
      unless ($rows2 == $cols2);

    croak "Math::MatrixReal::condition(): matrix size mismatch"
      unless (($rows1 == $rows2) && ($cols1 == $cols2));

    return( $matrix1->norm_one() * $matrix2->norm_one() );
}

## easy to use determinant
## very fast if matrix is diagonal or triangular

sub det {
    croak "Usage: \$determinant = \$matrix->det_LR();" unless (@_ == 1);
    my ($matrix) = @_;
    my ($rows,$cols) = $matrix->dim();
    my $det = 1;

    croak "Math::MatrixReal::det(): Matrix is not quadratic"
        unless ($rows == $cols);
    
    # diagonal will match too
    if( $matrix->is_upper_triangular() ){
        $matrix->each_diag( sub { $det*=shift; } );
    } elsif ( $matrix->is_lower_triangular() ){
        $matrix->each_diag( sub { $det*=shift; } );
    } else {
        return $matrix->decompose_LR->det_LR();
    }
    return $det;
}

sub det_LR  #  determinant of LR decomposition matrix
{
    croak "Usage: \$determinant = \$LR_matrix->det_LR();"
      if (@_ != 1);

    my($matrix) = @_;
    my($rows,$cols) = ($matrix->[1],$matrix->[2]);
    my($k,$det);

    croak "Math::MatrixReal::det_LR(): not an LR decomposition matrix"
      unless ((defined $matrix->[3]) && ($rows == $cols));

    $det = $matrix->[3];
    for ( $k = 0; $k < $rows; $k++ )
    {
        $det *= $matrix->[0][$k][$k];
    }
    return($det);
}

sub rank_LR {
    return (shift)->order_LR;
}

sub order_LR  #  order of LR decomposition matrix (number of non-zero equations)
{
    croak "Usage: \$order = \$LR_matrix->order_LR();"
      if (@_ != 1);

    my($matrix) = @_;
    my($rows,$cols) = ($matrix->[1],$matrix->[2]);
    my($order);

    croak "Math::MatrixReal::order_LR(): not an LR decomposition matrix"
      unless ((defined $matrix->[3]) && ($rows == $cols));

    ZERO:
    for ( $order = ($rows - 1); $order >= 0; $order-- )
    {
        last ZERO if ($matrix->[0][$order][$order] != 0);
    }
    return(++$order);
}

sub scalar_product
{
    croak "Usage: \$scalar_product = \$vector1->scalar_product(\$vector2);"
      if (@_ != 2);

    my($vector1,$vector2) = @_;
    my($rows1,$cols1) = ($vector1->[1],$vector1->[2]);
    my($rows2,$cols2) = ($vector2->[1],$vector2->[2]);
    my($k,$sum);

    croak "Math::MatrixReal::scalar_product(): 1st vector is not a column vector"
      unless ($cols1 == 1);

    croak "Math::MatrixReal::scalar_product(): 2nd vector is not a column vector"
      unless ($cols2 == 1);

    croak "Math::MatrixReal::scalar_product(): vector size mismatch"
      unless ($rows1 == $rows2);

    $sum = 0;
    for ( $k = 0; $k < $rows1; $k++ )
    {
        $sum += $vector1->[0][$k][0] * $vector2->[0][$k][0];
    }
    return($sum);
}

sub vector_product
{
    croak "Usage: \$vector_product = \$vector1->vector_product(\$vector2);"
      if (@_ != 2);

    my($vector1,$vector2) = @_;
    my($rows1,$cols1) = ($vector1->[1],$vector1->[2]);
    my($rows2,$cols2) = ($vector2->[1],$vector2->[2]);
    my($temp);
    my($n);

    croak "Math::MatrixReal::vector_product(): 1st vector is not a column vector"
      unless ($cols1 == 1);

    croak "Math::MatrixReal::vector_product(): 2nd vector is not a column vector"
      unless ($cols2 == 1);

    croak "Math::MatrixReal::vector_product(): vector size mismatch"
      unless ($rows1 == $rows2);

    $n = $rows1;

    croak "Math::MatrixReal::vector_product(): only defined for 3 dimensions"
      unless ($n == 3);

    $temp = $vector1->new($n,1);
    $temp->[0][0][0] = $vector1->[0][1][0] * $vector2->[0][2][0] -
                       $vector1->[0][2][0] * $vector2->[0][1][0];
    $temp->[0][1][0] = $vector1->[0][2][0] * $vector2->[0][0][0] -
                       $vector1->[0][0][0] * $vector2->[0][2][0];
    $temp->[0][2][0] = $vector1->[0][0][0] * $vector2->[0][1][0] -
                       $vector1->[0][1][0] * $vector2->[0][0][0];
    return($temp);
}

sub length
{
    croak "Usage: \$length = \$vector->length();"
      if (@_ != 1);

    my($vector) = @_;
    my($rows,$cols) = ($vector->[1],$vector->[2]);
    my($k,$comp,$sum);

    croak "Math::MatrixReal::length(): vector is not a column vector"
      unless ($cols == 1);

    $sum = 0;
    for ( $k = 0; $k < $rows; $k++ )
    {
        $comp = $vector->[0][$k][0];
        $sum += $comp * $comp;
    }
    return( sqrt( $sum ) );
}

sub _init_iteration
{
    croak "Usage: \$which_norm = \$matrix->_init_iteration();"
      if (@_ != 1);

    my($matrix) = @_;
    my($rows,$cols) = ($matrix->[1],$matrix->[2]);
    my($ok,$max,$sum,$norm);
    my($i,$j,$n);

    croak "Math::MatrixReal::_init_iteration(): matrix is not quadratic"
      unless ($rows == $cols);

    $ok = 1;
    $n = $rows;
    for ( $i = 0; $i < $n; $i++ )
    {
        if ($matrix->[0][$i][$i] == 0) { $ok = 0; }
    }
    if ($ok)
    {
        $norm = 1; # norm_one
        $max = 0;
        for ( $j = 0; $j < $n; $j++ )
        {
            $sum = 0;
            for ( $i = 0; $i < $j; $i++ )
            {
                $sum += abs($matrix->[0][$i][$j]);
            }
            for ( $i = ($j + 1); $i < $n; $i++ )
            {
                $sum += abs($matrix->[0][$i][$j]);
            }
            $sum /= abs($matrix->[0][$j][$j]);
            if ($sum > $max) { $max = $sum; }
        }
        $ok = ($max < 1);
        unless ($ok)
        {
            $norm = -1; # norm_max
            $max = 0;
            for ( $i = 0; $i < $n; $i++ )
            {
                $sum = 0;
                for ( $j = 0; $j < $i; $j++ )
                {
                    $sum += abs($matrix->[0][$i][$j]);
                }
                for ( $j = ($i + 1); $j < $n; $j++ )
                {
                    $sum += abs($matrix->[0][$i][$j]);
                }
                $sum /= abs($matrix->[0][$i][$i]);
                if ($sum > $max) { $max = $sum; }
            }
            $ok = ($max < 1)
        }
    }
    if ($ok) { return($norm); }
    else     { return(0); }
}

sub solve_GSM  #  Global Step Method
{
    croak "Usage: \$xn_vector = \$matrix->solve_GSM(\$x0_vector,\$b_vector,\$epsilon);"
      if (@_ != 4);

    my($matrix,$x0_vector,$b_vector,$epsilon) = @_;
    my($rows1,$cols1) = (   $matrix->[1],   $matrix->[2]);
    my($rows2,$cols2) = ($x0_vector->[1],$x0_vector->[2]);
    my($rows3,$cols3) = ( $b_vector->[1], $b_vector->[2]);
    my($norm,$sum,$diff);
    my($xn_vector);
    my($i,$j,$n);

    croak "Math::MatrixReal::solve_GSM(): matrix is not quadratic"
      unless ($rows1 == $cols1);

    $n = $rows1;

    croak "Math::MatrixReal::solve_GSM(): 1st vector is not a column vector"
      unless ($cols2 == 1);

    croak "Math::MatrixReal::solve_GSM(): 2nd vector is not a column vector"
      unless ($cols3 == 1);

    croak "Math::MatrixReal::solve_GSM(): matrix and vector size mismatch"
      unless (($rows2 == $n) && ($rows3 == $n));

    return() unless ($norm = $matrix->_init_iteration());

    $xn_vector = $x0_vector->new($n,1);

    $diff = $epsilon + 1;
    while ($diff >= $epsilon)
    {
        for ( $i = 0; $i < $n; $i++ )
        {
            $sum = $b_vector->[0][$i][0];
            for ( $j = 0; $j < $i; $j++ )
            {
                $sum -= $matrix->[0][$i][$j] * $x0_vector->[0][$j][0];
            }
            for ( $j = ($i + 1); $j < $n; $j++ )
            {
                $sum -= $matrix->[0][$i][$j] * $x0_vector->[0][$j][0];
            }
            $xn_vector->[0][$i][0] = $sum / $matrix->[0][$i][$i];
        }
        $x0_vector->subtract($x0_vector,$xn_vector);
        if ($norm > 0) { $diff = $x0_vector->norm_one(); }
        else           { $diff = $x0_vector->norm_max(); }
        for ( $i = 0; $i < $n; $i++ )
        {
            $x0_vector->[0][$i][0] = $xn_vector->[0][$i][0];
        }
    }
    return($xn_vector);
}

sub solve_SSM  #  Single Step Method
{
    croak "Usage: \$xn_vector = \$matrix->solve_SSM(\$x0_vector,\$b_vector,\$epsilon);"
      if (@_ != 4);

    my($matrix,$x0_vector,$b_vector,$epsilon) = @_;
    my($rows1,$cols1) = (   $matrix->[1],   $matrix->[2]);
    my($rows2,$cols2) = ($x0_vector->[1],$x0_vector->[2]);
    my($rows3,$cols3) = ( $b_vector->[1], $b_vector->[2]);
    my($norm,$sum,$diff);
    my($xn_vector);
    my($i,$j,$n);

    croak "Math::MatrixReal::solve_SSM(): matrix is not quadratic"
      unless ($rows1 == $cols1);

    $n = $rows1;

    croak "Math::MatrixReal::solve_SSM(): 1st vector is not a column vector"
      unless ($cols2 == 1);

    croak "Math::MatrixReal::solve_SSM(): 2nd vector is not a column vector"
      unless ($cols3 == 1);

    croak "Math::MatrixReal::solve_SSM(): matrix and vector size mismatch"
      unless (($rows2 == $n) && ($rows3 == $n));

    return() unless ($norm = $matrix->_init_iteration());

    $xn_vector = $x0_vector->new($n,1);
    $xn_vector->copy($x0_vector);

    $diff = $epsilon + 1;
    while ($diff >= $epsilon)
    {
        for ( $i = 0; $i < $n; $i++ )
        {
            $sum = $b_vector->[0][$i][0];
            for ( $j = 0; $j < $i; $j++ )
            {
                $sum -= $matrix->[0][$i][$j] * $xn_vector->[0][$j][0];
            }
            for ( $j = ($i + 1); $j < $n; $j++ )
            {
                $sum -= $matrix->[0][$i][$j] * $xn_vector->[0][$j][0];
            }
            $xn_vector->[0][$i][0] = $sum / $matrix->[0][$i][$i];
        }
        $x0_vector->subtract($x0_vector,$xn_vector);
        if ($norm > 0) { $diff = $x0_vector->norm_one(); }
        else           { $diff = $x0_vector->norm_max(); }
        for ( $i = 0; $i < $n; $i++ )
        {
            $x0_vector->[0][$i][0] = $xn_vector->[0][$i][0];
        }
    }
    return($xn_vector);
}

sub solve_RM  #  Relaxation Method
{
    croak "Usage: \$xn_vector = \$matrix->solve_RM(\$x0_vector,\$b_vector,\$weight,\$epsilon);"
      if (@_ != 5);

    my($matrix,$x0_vector,$b_vector,$weight,$epsilon) = @_;
    my($rows1,$cols1) = (   $matrix->[1],   $matrix->[2]);
    my($rows2,$cols2) = ($x0_vector->[1],$x0_vector->[2]);
    my($rows3,$cols3) = ( $b_vector->[1], $b_vector->[2]);
    my($norm,$sum,$diff);
    my($xn_vector);
    my($i,$j,$n);

    croak "Math::MatrixReal::solve_RM(): matrix is not quadratic"
      unless ($rows1 == $cols1);

    $n = $rows1;

    croak "Math::MatrixReal::solve_RM(): 1st vector is not a column vector"
      unless ($cols2 == 1);

    croak "Math::MatrixReal::solve_RM(): 2nd vector is not a column vector"
      unless ($cols3 == 1);

    croak "Math::MatrixReal::solve_RM(): matrix and vector size mismatch"
      unless (($rows2 == $n) && ($rows3 == $n));

    return() unless ($norm = $matrix->_init_iteration());

    $xn_vector = $x0_vector->new($n,1);
    $xn_vector->copy($x0_vector);

    $diff = $epsilon + 1;
    while ($diff >= $epsilon)
    {
        for ( $i = 0; $i < $n; $i++ )
        {
            $sum = $b_vector->[0][$i][0];
            for ( $j = 0; $j < $i; $j++ )
            {
                $sum -= $matrix->[0][$i][$j] * $xn_vector->[0][$j][0];
            }
            for ( $j = ($i + 1); $j < $n; $j++ )
            {
                $sum -= $matrix->[0][$i][$j] * $xn_vector->[0][$j][0];
            }
            $xn_vector->[0][$i][0] = $weight * ( $sum / $matrix->[0][$i][$i] )
                                   + (1 - $weight) * $xn_vector->[0][$i][0];
        }
        $x0_vector->subtract($x0_vector,$xn_vector);
        if ($norm > 0) { $diff = $x0_vector->norm_one(); }
        else           { $diff = $x0_vector->norm_max(); }
        for ( $i = 0; $i < $n; $i++ )
        {
            $x0_vector->[0][$i][0] = $xn_vector->[0][$i][0];
        }
    }
    return($xn_vector);
}

# Core householder reduction routine (when eigenvector
# are wanted).
# Adapted from: Numerical Recipes, 2nd edition.
sub _householder_vectors ($)
{
    my ($Q) = @_;
    my ($rows, $cols) = ($Q->[1], $Q->[2]);
    
    # Creates tridiagonal
    # Set up tridiagonal needed elements
    my $d = []; # N Diagonal elements 0...n-1
    my $e = []; # N-1 Off-Diagonal elements 0...n-2
    
    my @p = ();
    for (my $i = ($rows-1); $i > 1; $i--)
    {
    my $scale = 0.0;
    # Computes norm of one column (below diagonal)
    for (my $k = 0; $k < $i; $k++)
    {
        $scale += abs($Q->[0][$i][$k]);
    }
    if ($scale == 0.0)
    { # skip the transformation
        $e->[$i-1] = $Q->[0][$i][$i-1];
    }
    else
    {
        my $h = 0.0;
        for (my $k = 0; $k < $i; $k++)
        { # Used scaled Q for transformation
            $Q->[0][$i][$k] /= $scale;
            # Form sigma in h
            $h += $Q->[0][$i][$k] * $Q->[0][$i][$k];
        }
        my $t1 = $Q->[0][$i][$i-1];
        my $t2 = (($t1 >= 0.0) ? -sqrt($h) : sqrt($h));
        $e->[$i-1] = $scale * $t2; # Update off-diagonals
        $h -= $t1 * $t2;
        $Q->[0][$i][$i-1] -= $t2;
        my $f = 0.0;
        for (my $j = 0; $j < $i; $j++)
        {
            $Q->[0][$j][$i] = $Q->[0][$i][$j] / $h;
            my $g = 0.0;
            for (my $k = 0; $k <= $j; $k++)
            {
                $g += $Q->[0][$j][$k] * $Q->[0][$i][$k];
            }
            for (my $k = $j+1; $k < $i; $k++)
            {
                $g += $Q->[0][$k][$j] * $Q->[0][$i][$k];
            }
            # Form elements of P
            $p[$j] = $g / $h;
            $f += $p[$j] * $Q->[0][$i][$j];
        }
        my $hh = $f / ($h + $h);
        for (my $j = 0; $j < $i; $j++)
        {
            my $t3 = $Q->[0][$i][$j];
            my $t4 = $p[$j] - $hh * $t3;
            $p[$j] = $t4;
            for (my $k = 0; $k <= $j; $k++)
            {
                $Q->[0][$j][$k] -= $t3 * $p[$k]
                + $t4 * $Q->[0][$i][$k];
            }
        }
    }
    }
    # Updates for i == 0,1
    $e->[0] = $Q->[0][1][0];    
    $d->[0] = $Q->[0][0][0]; # i==0
    $Q->[0][0][0] = 1.0;
    $d->[1] = $Q->[0][1][1]; # i==1
    $Q->[0][1][1] = 1.0;
    $Q->[0][1][0] = $Q->[0][0][1] = 0.0;
    for (my $i = 2; $i < $rows; $i++)
    {
        for (my $j = 0; $j < $i; $j++)
        {
            my $g = 0.0;
            for (my $k = 0; $k < $i; $k++)
            {
                $g += $Q->[0][$i][$k] * $Q->[0][$k][$j];
            }
            for (my $k = 0; $k < $i; $k++)
            {
                $Q->[0][$k][$j] -= $g * $Q->[0][$k][$i];
            }
        }
        $d->[$i] = $Q->[0][$i][$i];
        # Reset row and column of Q for next iteration
        $Q->[0][$i][$i] = 1.0;
        for (my $j = 0; $j < $i; $j++)
        {
            $Q->[0][$i][$j] = $Q->[0][$j][$i] = 0.0;
        }
    }
    return ($d, $e);
}

# Computes sqrt(a*a + b*b), but more carefully...
sub _pythag ($$)
{
    my ($a, $b) = @_;
    my $aa = abs($a);
    my $ab = abs($b);
    if ($aa > $ab)
    {
        # NB: Not needed!: return 0.0 if ($aa == 0.0);
        my $t = $ab / $aa;
        return ($aa * sqrt(1.0 + $t*$t));
    }
    else
    {
        return 0.0 if ($ab == 0.0);
        my $t = $aa / $ab;
        return ($ab * sqrt(1.0 + $t*$t));
    }
}

# QL algorithm with implicit shifts to determine the eigenvalues
# of a tridiagonal matrix. Internal routine.
sub _tridiagonal_QLimplicit
{
    my ($EV, $d, $e) = @_;
    my ($rows, $cols) = ($EV->[1], $EV->[2]);

    $e->[$rows-1] = 0.0;
    # Start real computation
    for (my $l = 0; $l < $rows; $l++)
    {
        my $iter = 0;
        my $m;
        OUTER:
        do {
            for ($m = $l; $m < ($rows - 1); $m++)
            {
                my $dd = abs($d->[$m]) + abs($d->[$m+1]);
                last if ((abs($e->[$m]) + $dd) == $dd);
            }
            if ($m != $l)
            {
                ## why only allow 30 iterations?
                croak("Too many iterations!") if ($iter++ >= 30);
                my $g = ($d->[$l+1] - $d->[$l])
                    / (2.0 * $e->[$l]);
                my $r = _pythag($g, 1.0);
                $g = $d->[$m] - $d->[$l]
                    + $e->[$l] / ($g + (($g >= 0.0) ? abs($r) : -abs($r)));
                my ($p,$s,$c) = (0.0, 1.0,1.0);
            for (my $i = ($m-1); $i >= $l; $i--)
            {
                my $ii = $i + 1;
                my $f = $s * $e->[$i];
                my $t = _pythag($f, $g);
                $e->[$ii] = $t;
                if ($t == 0.0)
                {
                    $d->[$ii] -= $p;
                    $e->[$m] = 0.0;
                    next OUTER;
                }
                my $b = $c * $e->[$i];
                $s = $f / $t;
                $c = $g / $t;
                $g = $d->[$ii] - $p;
                my $t2 = ($d->[$i] - $g) * $s + 2.0 * $c * $b;
                $p = $s * $t2;
                $d->[$ii] = $g + $p;
                $g = $c * $t2 - $b;
                for (my $k = 0; $k < $rows; $k++)
                {
                    my $t1 = $EV->[0][$k][$ii];
                    my $t2 = $EV->[0][$k][$i];
                    $EV->[0][$k][$ii] = $s * $t2 + $c * $t1;
                    $EV->[0][$k][$i] = $c * $t2 - $s * $t1;
                }
            }
            $d->[$l] -= $p;
            $e->[$l] = $g;
            $e->[$m] = 0.0;
            }
        } while ($m != $l);
    }
    return;
}

# Core householder reduction routine (when eigenvector
# are NOT wanted).
sub _householder_values ($)
{
    my ($Q) = @_; # NB: Q is destroyed on output...
    my ($rows, $cols) = ($Q->[1], $Q->[2]);
    
    # Creates tridiagonal
    # Set up tridiagonal needed elements
    my $d = []; # N Diagonal elements 0...n-1
    my $e = []; # N-1 Off-Diagonal elements 0...n-2
    
    my @p = ();
    for (my $i = ($rows - 1); $i > 1; $i--)
    {
        my $scale = 0.0;
        for (my $k = 0; $k < $i; $k++)
        {
            $scale += abs($Q->[0][$i][$k]);
        }
        if ($scale == 0.0)
        { # skip the transformation
            $e->[$i-1] = $Q->[0][$i][$i-1];
        }
        else
        {
            my $h = 0.0;
            for (my $k = 0; $k < $i; $k++)
            { # Used scaled Q for transformation
                $Q->[0][$i][$k] /= $scale;
                # Form sigma in h
                $h += $Q->[0][$i][$k] * $Q->[0][$i][$k];
            }
            my $t = $Q->[0][$i][$i-1];
            my $t2 = (($t >= 0.0) ? -sqrt($h) : sqrt($h));
            $e->[$i-1] = $scale * $t2; # Updates off-diagonal
            $h -= $t * $t2;
            $Q->[0][$i][$i-1] -= $t2;
            my $f = 0.0;
            for (my $j = 0; $j < $i; $j++)
            {
                my $g = 0.0;
                for (my $k = 0; $k <= $j; $k++)
                {
                    $g += $Q->[0][$j][$k] * $Q->[0][$i][$k];
                }
                for (my $k = $j+1; $k < $i; $k++)
                {
                    $g += $Q->[0][$k][$j] * $Q->[0][$i][$k];
                }
                # Form elements of P
                $p[$j] = $g / $h;
                $f += $p[$j] * $Q->[0][$i][$j];
            }
            my $hh = $f / ($h + $h);
            for (my $j = 0; $j < $i; $j++)
            {
                my $t = $Q->[0][$i][$j];
                my $g = $p[$j] - $hh * $t;
                $p[$j] = $g;
                for (my $k = 0; $k <= $j; $k++)
                {
                    $Q->[0][$j][$k] -= $t * $p[$k]
                    + $g * $Q->[0][$i][$k];
                }
            }
        }
    }
    # Updates for i==1
    $e->[0] =  $Q->[0][1][0];
    # Updates diagonal elements
    for (my $i = 0; $i < $rows; $i++)
    {
        $d->[$i] =  $Q->[0][$i][$i];
    }
    return ($d, $e);
}

# QL algorithm with implicit shifts to determine the
# eigenvalues ONLY. This is O(N^2) only...
sub _tridiagonal_QLimplicit_values
{
    my ($M, $d, $e) = @_; # NB: M is not touched...
    my ($rows, $cols) = ($M->[1], $M->[2]);

    $e->[$rows-1] = 0.0;
    # Start real computation
    for (my $l = 0; $l < $rows; $l++)
    {
        my $iter = 0;
        my $m;
        OUTER:
        do {
            for ($m = $l; $m < ($rows - 1); $m++)
            {
                my $dd = abs($d->[$m]) + abs($d->[$m+1]);
                last if ((abs($e->[$m]) + $dd) == $dd);
            }
            if ($m != $l)
            {
                croak("Too many iterations!") if ($iter++ >= 30);
                my $g = ($d->[$l+1] - $d->[$l])
                    / (2.0 * $e->[$l]);
                my $r = _pythag($g, 1.0);
                $g = $d->[$m] - $d->[$l]
                    + $e->[$l] / ($g + (($g >= 0.0) ? abs($r) : -abs($r)));
                my ($p,$s,$c) = (0.0, 1.0,1.0);
                for (my $i = ($m-1); $i >= $l; $i--)
                {
                    my $ii = $i + 1;
                    my $f = $s * $e->[$i];
                    my $t = _pythag($f, $g);
                    $e->[$ii] = $t;
                    if ($t == 0.0)
                    {
                        $d->[$ii] -= $p;
                        $e->[$m] = 0.0;
                        next OUTER;
                    }
                    my $b = $c * $e->[$i];
                    $s = $f / $t;
                    $c = $g / $t;
                    $g = $d->[$ii] - $p;
                    my $t2 = ($d->[$i] - $g) * $s + 2.0 * $c * $b;
                    $p = $s * $t2;
                    $d->[$ii] = $g + $p;
                    $g = $c * $t2 - $b;
                }
                $d->[$l] -= $p;
                $e->[$l] = $g;
                $e->[$m] = 0.0;
            }
        } while ($m != $l);
    }
    return;
}

# Householder reduction of a real, symmetric matrix A.
# Returns a tridiagonal matrix T and the orthogonal matrix
# Q effecting the transformation between A and T.
sub householder ($)
{
    my ($A) = @_;
    my ($rows, $cols) = ($A->[1], $A->[2]);

    croak "Matrix is not quadratic"
        unless ($rows = $cols);
    croak "Matrix is not symmetric"
        unless ($A->is_symmetric());

    # Copy given matrix TODO: study if we should do in-place modification
    my $Q = $A->clone();

    # Do the computation of tridiagonal elements and of
    # transformation matrix
    my ($diag, $offdiag) = $Q->_householder_vectors();

    # Creates the tridiagonal matrix
    my $T = $A->shadow();
    for (my $i = 0; $i < $rows; $i++)
    { # Set diagonal
        $T->[0][$i][$i] = $diag->[$i];
    }
    for (my $i = 0; $i < ($rows-1); $i++)
    { # Set off diagonals
        $T->[0][$i+1][$i] = $offdiag->[$i];
        $T->[0][$i][$i+1] = $offdiag->[$i];
    }
    return ($T, $Q);
}

# QL algorithm with implicit shifts to determine the eigenvalues
# and eigenvectors of a real tridiagonal matrix - or of a matrix
# previously reduced to tridiagonal form.
sub tri_diagonalize ($;$)
{
    my ($T,$Q) = @_; # Q may be 0 if the original matrix is really tridiagonal

    my ($rows, $cols) = ($T->[1], $T->[2]);

    croak "Matrix is not quadratic"
        unless ($rows = $cols);
    croak "Matrix is not tridiagonal"
        unless ($T->is_tridiagonal()); # DONE

    my $EV;
    # Obtain/Creates the todo eigenvectors matrix
    if ($Q)
    {
        $EV = $Q->clone();
    }
    else
    {
        $EV = $T->shadow();
        $EV->one();
    }
    # Allocates diagonal vector
    my $diag = [ ];
    # Initializes it with T
    for (my $i = 0; $i < $rows; $i++)
    {
        $diag->[$i] = $T->[0][$i][$i];
    }
    # Allocate temporary vector for off-diagonal elements
    my $offdiag = [ ];
    for (my $i = 1; $i < $rows; $i++)
    {
        $offdiag->[$i-1] = $T->[0][$i][$i-1];
    }

    # Calls the calculus routine
    $EV->_tridiagonal_QLimplicit($diag, $offdiag);

    # Allocate eigenvalues vector
    my $v = Math::MatrixReal->new($rows,1);
    # Fills it
    for (my $i = 0; $i < $rows; $i++)
    {
        $v->[0][$i][0] = $diag->[$i];
    }
    return ($v, $EV);
}

# Main routine for diagonalization of a real symmetric
# matrix M. Operates by transforming M into a tridiagonal
# matrix and then obtaining the eigenvalues and eigenvectors
# for that matrix (taking into account the transformation to
# tridiagonal).
sub sym_diagonalize ($)
{
    my ($M) = @_;
    my ($rows, $cols) = ($M->[1], $M->[2]);
    
    croak "Matrix is not quadratic"
        unless ($rows = $cols);
    croak "Matrix is not symmetric"
        unless ($M->is_symmetric());
    
    # Copy initial matrix
    # TODO: study if we should allow in-place modification
    my $VEC = $M->clone();

    # Do the computation of tridiagonal elements and of
    # transformation matrix
    my ($diag, $offdiag) = $VEC->_householder_vectors();
    # Calls the calculus routine for diagonalization
    $VEC->_tridiagonal_QLimplicit($diag, $offdiag);

    # Allocate eigenvalues vector
    my $val = Math::MatrixReal->new($rows,1);
    # Fills it
    for (my $i = 0; $i < $rows; $i++)
    {
        $val->[0][$i][0] = $diag->[$i];
    }
    return ($val, $VEC);
}

# Householder reduction of a real, symmetric matrix A.
# Returns a tridiagonal matrix T equivalent to A.
sub householder_tridiagonal ($)
{
    my ($A) = @_;
    my ($rows, $cols) = ($A->[1], $A->[2]);

    croak "Matrix is not quadratic"
        unless ($rows = $cols);
    croak "Matrix is not symmetric"
        unless ($A->is_symmetric());

    # Copy given matrix
    my $Q = $A->clone();

    # Do the computation of tridiagonal elements and of
    # transformation matrix
    # Q is destroyed after reduction
    my ($diag, $offdiag) = $Q->_householder_values();

    # Creates the tridiagonal matrix in Q (avoid allocation)
    my $T = $Q;
    $T->zero();
    for (my $i = 0; $i < $rows; $i++)
    { # Set diagonal
        $T->[0][$i][$i] = $diag->[$i];
    }
    for (my $i = 0; $i < ($rows-1); $i++)
    { # Set off diagonals
        $T->[0][$i+1][$i] = $offdiag->[$i];
        $T->[0][$i][$i+1] = $offdiag->[$i];
    }
    return $T;
}

# QL algorithm with implicit shifts to determine ONLY
# the eigenvalues a real tridiagonal matrix - or of a
# matrix previously reduced to tridiagonal form.
sub tri_eigenvalues ($;$)
{
    my ($T) = @_;
    my ($rows, $cols) = ($T->[1], $T->[2]);

    croak "Matrix is not quadratic"
        unless ($rows = $cols);
    croak "Matrix is not tridiagonal"
        unless ($T->is_tridiagonal() ); # DONE

    # Allocates diagonal vector
    my $diag = [ ];
    # Initializes it with T
    for (my $i = 0; $i < $rows; $i++)
    {
        $diag->[$i] = $T->[0][$i][$i];
    }
    # Allocate temporary vector for off-diagonal elements
    my $offdiag = [ ];
    for (my $i = 1; $i < $rows; $i++)
    {
        $offdiag->[$i-1] = $T->[0][$i][$i-1];
    }

    # Calls the calculus routine (T is not touched)
    $T->_tridiagonal_QLimplicit_values($diag, $offdiag);

    # Allocate eigenvalues vector
    my $v = Math::MatrixReal->new($rows,1);
    # Fills it
    for (my $i = 0; $i < $rows; $i++)
    {
        $v->[0][$i][0] = $diag->[$i];
    }
    return $v;
}

## more general routine than sym_eigenvalues
sub eigenvalues ($){
    my ($matrix) = @_;
    my ($rows,$cols) = $matrix->dim();
    my $i=0;

    #return sym_eigenvalues($matrix) if $matrix->is_symmetric();
    
    croak "Matrix is not quadratic" unless ($rows == $cols);

    if($matrix->is_upper_triangular() || $matrix->is_lower_triangular() ){
        my $l = Math::MatrixReal->new($rows,1);
        for(;$i < $rows; $i++ ){
            $l->[0][$i][0] = $matrix->[0][$i][$i];
        }
        return $l;
    }

    return sym_eigenvalues($matrix) if $matrix->is_symmetric();

    carp "Math::MatrixReal::eigenvalues(): Matrix is not symmetric or triangular";
    return undef;

}
# Main routine for diagonalization of a real symmetric
# matrix M. Operates by transforming M into a tridiagonal
# matrix and then obtaining the eigenvalues and eigenvectors
# for that matrix (taking into account the transformation to
# tridiagonal).
sub sym_eigenvalues ($)
{
    my ($M) = @_;
    my ($rows, $cols) = ($M->[1], $M->[2]);
    
    croak "Matrix is not quadratic"
        unless ($rows == $cols); 
    croak "Matrix is not symmetric"
        unless ($M->is_symmetric());

    # Copy matrix in temporary
    my $A = $M->clone();
    # Do the computation of tridiagonal elements and of
    # transformation matrix. A is destroyed
    my ($diag, $offdiag) = $A->_householder_values();
    # Calls the calculus routine for diagonalization
    # (M is not touched)
    $M->_tridiagonal_QLimplicit_values($diag, $offdiag);

    # Allocate eigenvalues vector
    my $val = Math::MatrixReal->new($rows,1);
    # Fills it
    for (my $i = 0; $i < $rows; $i++)
    {
        $val->[0][$i][0] = $diag->[$i];
    }
    return $val;
}
#TODO: docs+test
sub is_positive_definite {
    my ($matrix) = @_;
    my ($r,$c) = $matrix->dim;
    my $ev = $matrix->eigenvalues;
    my $pos = 1;
    $ev->each(sub { my $x = shift; if ($x <= 0){ $pos=0;return; } } );
    return $pos;
}

sub is_row_vector {
    my ($m) = @_;
    my ($r,$c) = $m->dim;
    return 1 if ($r == 1);
}
sub is_col_vector {
    my ($m) = @_;
    my ($r,$c) = $m->dim;
    return 1 if ($c == 1);
}

sub is_orthogonal($) {
    my ($matrix) = @_;
    ##croak "Math::MatrixReal::is_orthogonal(): Matrix is not quadratic"

    return 0 unless( $matrix->is_quadratic() );

    my $one = $matrix->shadow();
    $one->one;

    abs(~$matrix * $matrix - $one) < 1e-12 ? return 1 : return 0;

}

sub is_positive($) {
    my ($m) = @_;
    my $pos = 1;
    $m->each( sub { if( (shift) <= 0){ $pos = 0;return;} } );
    return $pos;
}
sub is_negative($) {
    my ($m) = @_;
    my $neg = 1;
    $m->each( sub { if( (shift) >= 0){ $neg = 0;return;} } );
    return $neg;
}


sub is_periodic($$) {
    my ($m,$k) = @_;
        return 0 unless $m->is_quadratic();
    abs($m**(int($k)+1) - $m) < 1e-12 ? return 1 : return 0;
}
sub is_idempotent($) {
    return (shift)->is_periodic(1);
}

# Boolean check routine to see if a matrix is
# symmetric
sub is_symmetric ($)
{
  my ($M) = @_;
  my ($rows, $cols) = ($M->[1], $M->[2]);
  # if it is not quadratic it cannot be symmetric...
  return 0 unless ($rows == $cols);
  # skip when $i=$j?
  for (my $i = 1; $i < $rows; $i++)
    {
      for (my $j = 0; $j < $i; $j++)
        {
          return 0 unless ($M->[0][$i][$j] == $M->[0][$j][$i]);
        }
    }
  return 1;
}
# Boolean check to see if matrix is tridiagonal
sub is_tridiagonal ($) {
    my ($M) = @_;
    my ($rows,$cols) = ($M->[1],$M->[2]);
    my ($i,$j) = (0,0); 
    # if it is not quadratic it cannot be tridiag
    return 0 unless ($rows == $cols);

    for(;$i < $rows; $i++ ){
        for(;$j < $cols; $j++ ){
            #print "debug: testing $i,$j = " . $M->[0][$i][$j] . "\n";
            # skip diag and diag+-1
            next if ($i == $j);
            next if ($i+1 == $j);
            next if ($i-1 == $j);
            return 0 if $M->[0][$i][$j];
        }
        $j = 0;
    }
    return 1;
}
# Boolean check to see if matrix is upper triangular
# i.e all nonzero elements are above main diagonal
sub is_upper_triangular {
    my ($M) = @_;
    my ($rows,$cols) = $M->dim();
    my ($i,$j) = (1,0);
    return 0 unless ($rows == $cols);
    
    for(;$i < $rows; $i++ ){
        for(;$j < $cols;$j++ ){
            next if ($i <= $j);
            return 0 if $M->[0][$i][$j];
        }
        $j = 0;
    }
    return 1;
}
# Boolean check to see if matrix is lower triangular
# i.e all nonzero elements are lower main diagonal
sub is_lower_triangular {
    my ($M) = @_;
    my ($rows,$cols) = $M->dim();
    my ($i,$j) = (0,1);
    return 0 unless ($rows == $cols);

    for(;$i < $rows; $i++ ){
        for(;$j < $cols;$j++ ){
            next if ($i >= $j);
            return 0 if $M->[0][$i][$j];
        }
        $j = 0;
    }
    return 1;
}

# Boolean check to see if matrix is diagonal
sub is_diagonal ($) {
    my ($M) = @_;
    my ($rows,$cols) = ($M->[1],$M->[2]);
    my ($i,$j) = (0,0);
    return 0 unless ($rows == $cols );
    for(;$i < $rows; $i++ ){
        for(;$j < $cols; $j++ ){
            # skip diag elements
            next if ($i == $j);
            return 0 if $M->[0][$i][$j];
        }
        $j = 0;
    }
    return 1;
}
sub is_quadratic ($) {
    croak "Usage: \$matrix->is_quadratic()" unless (@_ == 1);
    my ($matrix) = @_;
    $matrix->[1] == $matrix->[2] ? return 1 : return 0;
}

sub is_square($) {
    croak "Usage: \$matrix->is_square()" unless (@_ == 1);
    return (shift)->is_quadratic();
}

sub is_LR($) {
    croak "Usage: \$matrix->is_LR()" unless (@_ == 1);
    return (shift)->[3] ? 1 : 0;
}
###
sub is_normal{
    my ($matrix) = @_;
    my ($rows,$cols) = $matrix->dim;
    return 0 unless ($rows == $cols);
    
    return 1 if ( ~$matrix * $matrix - $matrix * ~$matrix < 1e-8 );
    return 0;

}
sub is_skew_symmetric{
    my ($m) = @_;
    my ($rows, $cols) = $m->dim;
    # if it is not quadratic it cannot be skew symmetric...
    return 0 unless ($rows == $cols);
    for (my $i = 1; $i < $rows; $i++) {
        for (my $j = 0; $j < $i; $j++) {
            return 0 unless ($m->[0][$i][$j] == -$m->[0][$j][$i]);
        }
    }
    return 1;

}
####
sub is_gramian{
    my ($m) = @_;
    my ($rows,$cols) = $m->dim;
    my $neg=0;
    # gramian matrix must be symmetric
    return 0 unless $m->is_symmetric;

    # must have all non-negative eigenvalues
    my $ev = $m->eigenvalues;
    $ev->each(sub { $neg++ if ((shift)<0) } );

    return $neg ? 0 : 1;
}
sub is_binary{
    my ($m) = @_;
    my ($rows, $cols) = $m->dim;
    for (my $i = 0; $i < $rows; $i++) {
        for (my $j = 0; $j < $cols; $j++) {
            return 0 unless ($m->[0][$i][$j] == 1 || $m->[0][$i][$j] == 0);
        }
    }
    return 1;

}
sub as_scilab {
    return (shift)->as_matlab;
}

sub as_matlab {
    my ($m) = shift;
    my %args = ( 
        format => "%s",
        name => "",
        semi => 0,
        @_);
    my ($row,$col) = $m->dim;
    my $s = "";
    
    if( $args{name} ){
        $s = "$args{name} = ";
    }
    $s .= "[";
    $m->each( 
        sub { my($x,$i,$j) = @_;
            $s .= sprintf(" $args{format}",$x);
            $s .= ";\n" if( $j == $col && $i != $row);
        }
        );
    $s .= "]";
    $s .= ";" if $args{semi};
    return $s;
}
#TODO: docs+test
sub as_yacas{
    my ($m) = shift;
    my %args = (
        format => "%s",
        name => "",
        semi => 0,
        @_);
    my ($row,$col) = $m->dim;
    my $s = "";

    if( $args{name} ){
        $s = "$args{name} := ";
    }
    $s .= "{";

    $m->each(
        sub { my($x,$i,$j) = @_;
            $s .= "{" if ($j == 1);
                        $s .= sprintf("$args{format}",$x);
                        $s .= "," if( $j != $col );
            $s .= "}," if ($j == $col && $i != $row);
        }
    );
    $s .= "}}";

    return $s;
}
#TODO: any ideas for a good test?
sub as_latex{
    my ($m) = shift;
    my %args = ( 
        format => "%s",
        name => "",
        align => "c",
        @_);


    my ($row,$col) = $m->dim;
    my $inside;
    my $s = <<LATEX;
\$
\\left( \\begin{array}{%COLS%}
%INSIDE%\\end{array} \\right)
\$
LATEX
    $args{align} = lc $args{align};
    if( $args{align} !~ m/^(c|l|r)$/ ){
        croak "Math::MatrixReal::as_latex(): Invalid alignment '$args{align}'";
    }

    $s =~ s/%COLS%/$args{align} x $col/em;

    if( $args{name} ){
        $s = "\$$args{name} = \$ " . $s;
    }
    $m->each( 
        sub { my($x,$i,$j)=@_;

            $x = sprintf($args{format},$x);

            # last element in each row gets a \\
            if ($j == $col && $i != $row){
                $inside .= "$x \\\\"."\n";
            # the annoying last line has neither
            } elsif( $j == $col && $i == $row){ 
                $inside .= "$x\n";
            } else {
                $inside .= "$x&";
            }
        } 
    );
    $s =~ s/%INSIDE%/$inside/gm;
    return $s;
}
#### 
sub spectral_radius {
    my ($matrix) = @_;
    my ($r,$c) = $matrix->dim;
    my $ev = $matrix->eigenvalues;
    my $radius=0;
    $ev->each(sub { my $x = shift; $radius = $x if (abs($x) > $radius); } );
    return $radius;
}


                ########################################
                #                                      #
                # define overloaded operators section: #
                #                                      #
                ########################################

sub _negate
{
    my($object,$argument,$flag) = @_;
#   my($name) = "neg"; 
    my($temp);

    $temp = $object->new($object->[1],$object->[2]);
    $temp->negate($object);
    return($temp);
}

sub _transpose
{
    my($object,$argument,$flag) = @_;
#   my($name) = "'~'"; 
    my($temp);

    $temp = $object->new($object->[2],$object->[1]);
    $temp->transpose($object);
    return($temp);

}

sub _boolean
{
    my($object,$argument,$flag) = @_;
#   my($name) = "bool"; 
    my($rows,$cols) = ($object->[1],$object->[2]);
    my($i,$j,$result);

    #TODO: ugly...
    $result = 0;
    BOOL:
    for ( $i = 0; $i < $rows; $i++ )
    {
        for ( $j = 0; $j < $cols; $j++ )
        {
            if ($object->[0][$i][$j] != 0)
            {
                $result = 1;
                last BOOL;
            }
        }
    }
    return($result);
}
#TODO: ugly copy+paste
sub _not_boolean
{
    my($object,$argument,$flag) = @_;
#   my($name) = "'!'"; 
    my($rows,$cols) = ($object->[1],$object->[2]);
    my($i,$j,$result);

    $result = 1;
    NOTBOOL:
    for ( $i = 0; $i < $rows; $i++ )
    {
        for ( $j = 0; $j < $cols; $j++ )
        {
            if ($object->[0][$i][$j] != 0)
            {
                $result = 0;
                last NOTBOOL;
            }
        }
    }
    return($result);
}

sub _stringify
{
    my($object,$argument,$flag) = @_;
#   my($name) = '""'; 
    my($rows,$cols) = ($object->[1],$object->[2]);
    my($i,$j,$s);

    $s = '';
    for ( $i = 0; $i < $rows; $i++ )
    {
        $s .= "[ ";
        for ( $j = 0; $j < $cols; $j++ )
        {
            $s .= sprintf("% #-19.12E ", $object->[0][$i][$j]);
        }
        $s .= "]\n";
    }
    return($s);
}

sub _norm
{
    my($object,$argument,$flag) = @_;
#   my($name) = "abs"; 

    return( $object->norm_one() );
}

sub _add
{
    my($object,$argument,$flag) = @_;
    my($name) = "'+'"; 
    my($temp);

    if ((defined $argument) && ref($argument) &&
        (ref($argument) !~ /^SCALAR$|^ARRAY$|^HASH$|^CODE$|^REF$/))
    {
        if (defined $flag)
        {
            $temp = $object->new($object->[1],$object->[2]);
            $temp->add($object,$argument);
            return($temp);
        }
        else
        {
            $object->add($object,$argument);
            return($object);
        }
    }
    else
    {
        croak "Math::MatrixReal $name: wrong argument type";
    }
}

sub _subtract
{
    my($object,$argument,$flag) = @_;
    my($name) = "'-'"; 
    my($temp);

    if ((defined $argument) && ref($argument) &&
        (ref($argument) !~ /^SCALAR$|^ARRAY$|^HASH$|^CODE$|^REF$/))
    {
        if (defined $flag)
        {
            $temp = $object->new($object->[1],$object->[2]);
            if ($flag) { $temp->subtract($argument,$object); }
            else       { $temp->subtract($object,$argument); }
            return($temp);
        }
        else
        {
            $object->subtract($object,$argument);
            return($object);
        }
    }
    else
    {
        croak "Math::MatrixReal $name: wrong argument type";
    }
}

sub _exponent 
{
    my($matrix,$argument,$flag) = @_;
    my($rows,$cols) = ($matrix->[1],$matrix->[2]);
    my($name) = "'**'"; 

    return $matrix->exponent( $argument );
}


sub _multiply
{
    my($object,$argument,$flag) = @_;
    my($name) = "'*'"; 
    my($temp);

    if ((defined $argument) && ref($argument) &&
        (ref($argument) !~ /^SCALAR$|^ARRAY$|^HASH$|^CODE$|^REF$/))
    {
        if ((defined $flag) && $flag)
        {
            return( multiply($argument,$object) );
        }
        else
        {
            return( multiply($object,$argument) );
        }
    }
    elsif ((defined $argument) && !(ref($argument)))
    {
        if (defined $flag)
        {
            $temp = $object->new($object->[1],$object->[2]);
            $temp->multiply_scalar($object,$argument);
            return($temp);
        }
        else
        {
            $object->multiply_scalar($object,$argument);
            return($object);
        }
    }
    else
    {
        croak "Math::MatrixReal $name: wrong argument type";
    }
}

sub _assign_add
{
    my($object,$argument,$flag) = @_;
#   my($name) = "'+='"; 

    return( &_add($object,$argument,undef) );
}

sub _assign_subtract
{
    my($object,$argument,$flag) = @_;
#   my($name) = "'-='"; 

    return( &_subtract($object,$argument,undef) );
}

sub _assign_multiply
{
    my($object,$argument,$flag) = @_;
#   my($name) = "'*='"; 

    return( &_multiply($object,$argument,undef) );
}

sub _assign_exponent {
    my($object,$arg,$flag) = @_;
    return ( &_exponent($object,$arg,undef) );
}

sub _equal
{
    my($object,$argument,$flag) = @_;
    my($name) = "'=='"; 
    my($rows,$cols) = ($object->[1],$object->[2]);
    my($i,$j,$result);

    if ((defined $argument) && ref($argument) &&
        (ref($argument) !~ /^SCALAR$|^ARRAY$|^HASH$|^CODE$|^REF$/))
    {
        $result = 1;
        EQUAL:
        for ( $i = 0; $i < $rows; $i++ )
        {
            for ( $j = 0; $j < $cols; $j++ )
            {
                if ($object->[0][$i][$j] != $argument->[0][$i][$j])
                {
                    $result = 0;
                    last EQUAL;
                }
            }
        }
        return($result);
    }
    else
    {
        croak "Math::MatrixReal $name: wrong argument type";
    }
}

sub _not_equal
{
    my($object,$argument,$flag) = @_;
    my($name) = "'!='"; 
    my($rows,$cols) = ($object->[1],$object->[2]);
    my($i,$j,$result);

    if ((defined $argument) && ref($argument) &&
        (ref($argument) !~ /^SCALAR$|^ARRAY$|^HASH$|^CODE$|^REF$/))
    {
        $result = 0;
        NOTEQUAL:
        for ( $i = 0; $i < $rows; $i++ )
        {
            for ( $j = 0; $j < $cols; $j++ )
            {
                if ($object->[0][$i][$j] != $argument->[0][$i][$j])
                {
                    $result = 1;
                    last NOTEQUAL;
                }
            }
        }
        return($result);
    }
    else
    {
        croak "Math::MatrixReal $name: wrong argument type";
    }
}

sub _less_than
{
    my($object,$argument,$flag) = @_;
    my($name) = "'<'"; 

    if ((defined $argument) && ref($argument) &&
        (ref($argument) !~ /^SCALAR$|^ARRAY$|^HASH$|^CODE$|^REF$/))
    {
        if ((defined $flag) && $flag)
        {
            return( $argument->norm_one() < $object->norm_one() );
        }
        else
        {
            return( $object->norm_one() < $argument->norm_one() );
        }
    }
    elsif ((defined $argument) && !(ref($argument)))
    {
        if ((defined $flag) && $flag)
        {
            return( abs($argument) < $object->norm_one() );
        }
        else
        {
            return( $object->norm_one() < abs($argument) );
        }
    }
    else
    {
        croak "Math::MatrixReal $name: wrong argument type";
    }
}

sub _less_than_or_equal
{
    my($object,$argument,$flag) = @_;
    my($name) = "'<='"; 

    if ((defined $argument) && ref($argument) &&
        (ref($argument) !~ /^SCALAR$|^ARRAY$|^HASH$|^CODE$|^REF$/))
    {
        if ((defined $flag) && $flag)
        {
            return( $argument->norm_one() <= $object->norm_one() );
        }
        else
        {
            return( $object->norm_one() <= $argument->norm_one() );
        }
    }
    elsif ((defined $argument) && !(ref($argument)))
    {
        if ((defined $flag) && $flag)
        {
            return( abs($argument) <= $object->norm_one() );
        }
        else
        {
            return( $object->norm_one() <= abs($argument) );
        }
    }
    else
    {
        croak "Math::MatrixReal $name: wrong argument type";
    }
}

sub _greater_than
{
    my($object,$argument,$flag) = @_;
    my($name) = "'>'"; 

    if ((defined $argument) && ref($argument) &&
        (ref($argument) !~ /^SCALAR$|^ARRAY$|^HASH$|^CODE$|^REF$/))
    {
        if ((defined $flag) && $flag)
        {
            return( $argument->norm_one() > $object->norm_one() );
        }
        else
        {
            return( $object->norm_one() > $argument->norm_one() );
        }
    }
    elsif ((defined $argument) && !(ref($argument)))
    {
        if ((defined $flag) && $flag)
        {
            return( abs($argument) > $object->norm_one() );
        }
        else
        {
            return( $object->norm_one() > abs($argument) );
        }
    }
    else
    {
        croak "Math::MatrixReal $name: wrong argument type";
    }
}

sub _greater_than_or_equal
{
    my($object,$argument,$flag) = @_;
    my($name) = "'>='"; 

    if ((defined $argument) && ref($argument) &&
        (ref($argument) !~ /^SCALAR$|^ARRAY$|^HASH$|^CODE$|^REF$/))
    {
        if ((defined $flag) && $flag)
        {
            return( $argument->norm_one() >= $object->norm_one() );
        }
        else
        {
            return( $object->norm_one() >= $argument->norm_one() );
        }
    }
    elsif ((defined $argument) && !(ref($argument)))
    {
        if ((defined $flag) && $flag)
        {
            return( abs($argument) >= $object->norm_one() );
        }
        else
        {
            return( $object->norm_one() >= abs($argument) );
        }
    }
    else
    {
        croak "Math::MatrixReal $name: wrong argument type";
    }
}
sub _exp {
    my ($matrix,$arg,$flag) = @_;
    my $new_matrix = $matrix->clone();
    my ($rows,$cols) = $matrix->dim();
    
    $new_matrix->_undo_LR();

    croak "Math::MatrixReal::exp(): Matrix is not quadratic" unless ($rows == $cols);
    croak "Math::MatrixReal::exp(): Only diagonal matrices supported" unless ( $matrix->is_diagonal() );

    $new_matrix = $matrix->each_diag( sub { exp(shift) } );
    return $new_matrix;

}
sub _cos {
    my ($matrix,$arg,$flag) = @_;
    my $new_matrix = $matrix->clone();
    my ($rows,$cols) = $matrix->dim();

    $new_matrix->_undo_LR();

    croak "Math::MatrixReal::cos(): Matrix is not quadratic" unless ($rows == $cols);
    croak "Math::MatrixReal::cos(): Only diagonal matrices supported" unless ( $matrix->is_diagonal() );

    $new_matrix = $matrix->each_diag( sub { cos(shift) } );
    return $new_matrix;
}
sub _sin {
    my ($matrix,$arg,$flag) = @_;
    my $new_matrix = $matrix->clone();
    my ($rows,$cols) = $matrix->dim();

    $new_matrix->_undo_LR();

    croak "Math::MatrixReal::sin(): Matrix is not quadratic" unless ($rows == $cols);
    croak "Math::MatrixReal::sin(): Only diagonal matrices supported" unless ( $matrix->is_diagonal() );

    $new_matrix = $matrix->each_diag( sub { sin(shift) } );
    return $new_matrix;

}

sub _clone
{
    my($object,$argument,$flag) = @_;
#   my($name) = "'='"; 
    my($temp);

    $temp = $object->new($object->[1],$object->[2]);
    $temp->copy($object);
    $temp->_undo_LR();
    return($temp);
}

sub _trace
{
    my($text,$object,$argument,$flag) = @_;

    unless (defined $object)   { $object   = 'undef'; };
    unless (defined $argument) { $argument = 'undef'; };
    unless (defined $flag)     { $flag     = 'undef'; };
    if (ref($object))   { $object   = ref($object);   }
    if (ref($argument)) { $argument = ref($argument); }
    print "$text: \$obj='$object' \$arg='$argument' \$flag='$flag'\n";
}

1;

__END__

=head1 NAME

Math::MatrixReal - Matrix of Reals

Implements the data type "matrix of reals" (and consequently also
"vector of reals").

=head1 DESCRIPTION

Implements the data type "matrix of reals", which can be used almost
like any other basic Perl type thanks to B<OPERATOR OVERLOADING>, i.e.,

  $product = $matrix1 * $matrix2;

does what you would like it to do (a matrix multiplication).

Also features many important operations and methods: matrix norm,
matrix transposition, matrix inverse, determinant of a matrix, order
and numerical condition of a matrix, scalar product of vectors, vector
product of vectors, vector length, projection of row and column vectors,
a comfortable way for reading in a matrix from a file, the keyboard or
your code, and many more.

Allows to solve linear equation systems using an efficient algorithm
known as "L-R-decomposition" and several approximative (iterative) methods.

Features an implementation of Kleene's algorithm to compute the minimal
costs for all paths in a graph with weighted edges (the "weights" being
the costs associated with each edge).

=head1 SYNOPSIS

=head2 Constructor Methods And Such

=over 4
=item *

C<use Math::MatrixReal;>

Makes the methods and overloaded operators of this module available
to your program.

=item *

C<$new_matrix = new Math::MatrixReal($rows,$columns);>

The matrix object constructor method. A new matrix of size $rows by $columns
will be created, with the value C<0.0> for all elements.

Note that this method is implicitly called by many of the other methods
in this module.

=item *

C<$new_matrix = $some_matrix-E<gt>>C<new($rows,$columns);>

Another way of calling the matrix object constructor method.

Matrix "C<$some_matrix>" is not changed by this in any way.

=item *

C<$new_matrix = $matrix-E<gt>new_from_cols( [ $column_vector|$array_ref|$string, ... ] )>

Creates a new matrix given a reference to an array of any of the following:

=over 4

=item * column vectors ( n by 1 Math::MatrixReal matrices )

=item * references to arrays

=item * strings properly formatted to create a column with Math::MatrixReal's
C<new_from_string> command

=back

You may mix and match these as you wish.  However, all must be of the
same dimension--no padding happens automatically.  Example: 

    my $matrix = Math::MatrixReal->new_from_cols( [ [1,2], [3,4] ] );
    print $matrix;

will print

    [  1.000000000000E+00  3.000000000000E+00 ]
    [  2.000000000000E+00  4.000000000000E+00 ]


=item * 

C<new_from_rows( [ $row_vector|$array_ref|$string, ... ] )>

Creates a new matrix given a reference to an array of any of the following:

=over 4

=item * row vectors ( 1 by n Math::MatrixReal matrices )

=item * references to arrays

=item * strings properly formatted to create a row with Math::MatrixReal's C<new_from_string> command

=back

You may mix and match these as you wish.  However, all must be of the
same dimension--no padding happens automatically. Example:

        my $matrix = Math::MatrixReal->new_from_rows( [ [1,2], [3,4] ] );
        print $matrix;

will print

        [  1.000000000000E+00  2.000000000000E+00 ]
        [  3.000000000000E+00  4.000000000000E+00 ]


=item *

C<$new_matrix = Math::MatrixReal-E<gt>new_diag( $array_ref );>

This method allows you to create a diagonal matrix by only specifying
the diagonal elements. Example: 

    $matrix = Math::MatrixReal->new_diag( [ 1,2,3,4 ] );
    print $matrix;

will print

    [  1.000000000000E+00  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00 ]
    [  0.000000000000E+00  2.000000000000E+00  0.000000000000E+00  0.000000000000E+00 ]
    [  0.000000000000E+00  0.000000000000E+00  3.000000000000E+00  0.000000000000E+00 ]
    [  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00  4.000000000000E+00 ]


=item *

C<$new_matrix = Math::MatrixReal-E<gt>>C<new_from_string($string);>

This method allows you to read in a matrix from a string (for
instance, from the keyboard, from a file or from your code).

The syntax is simple: each row must start with "C<[ >" and end with
"C< ]\n>" ("C<\n>" being the newline character and "C< >" a space or
tab) and contain one or more numbers, all separated from each other
by spaces or tabs.

Additional spaces or tabs can be added at will, but no comments.

Examples:

  $string = "[ 1 2 3 ]\n[ 2 2 -1 ]\n[ 1 1 1 ]\n";
  $matrix = Math::MatrixReal->new_from_string($string);
  print "$matrix";

By the way, this prints

  [  1.000000000000E+00  2.000000000000E+00  3.000000000000E+00 ]
  [  2.000000000000E+00  2.000000000000E+00 -1.000000000000E+00 ]
  [  1.000000000000E+00  1.000000000000E+00  1.000000000000E+00 ]

But you can also do this in a much more comfortable way using the
shell-like "here-document" syntax:

  $matrix = Math::MatrixReal->new_from_string(<<'MATRIX');
  [  1  0  0  0  0  0  1  ]
  [  0  1  0  0  0  0  0  ]
  [  0  0  1  0  0  0  0  ]
  [  0  0  0  1  0  0  0  ]
  [  0  0  0  0  1  0  0  ]
  [  0  0  0  0  0  1  0  ]
  [  1  0  0  0  0  0 -1  ]
  MATRIX

You can even use variables in the matrix:

  $c1 =   2  /  3;
  $c2 =  -2  /  5;
  $c3 =  26  /  9;

  $matrix = Math::MatrixReal->new_from_string(<<"MATRIX");

      [   3    2    0   ]
      [   0    3    2   ]
      [  $c1  $c2  $c3  ]

  MATRIX

(Remember that you may use spaces and tabs to format the matrix to
your taste)

Note that this method uses exactly the same representation for a
matrix as the "stringify" operator "": this means that you can convert
any matrix into a string with C<$string = "$matrix";> and read it back
in later (for instance from a file!).

Note however that you may suffer a precision loss in this process
because only 13 digits are supported in the mantissa when printed!!

If the string you supply (or someone else supplies) does not obey
the syntax mentioned above, an exception is raised, which can be
caught by "eval" as follows:

  print "Please enter your matrix (in one line): ";
  $string = <STDIN>;
  $string =~ s/\\n/\n/g;
  eval { $matrix = Math::MatrixReal->new_from_string($string); };
  if ($@)
  {
      print "$@";
      # ...
      # (error handling)
  }
  else
  {
      # continue...
  }

or as follows:

  eval { $matrix = Math::MatrixReal->new_from_string(<<"MATRIX"); };
  [   3    2    0   ]
  [   0    3    2   ]
  [  $c1  $c2  $c3  ]
  MATRIX
  if ($@)
  # ...

Actually, the method shown above for reading a matrix from the keyboard
is a little awkward, since you have to enter a lot of "\n"'s for the
newlines.

A better way is shown in this piece of code:

  while (1)
  {
      print "\nPlease enter your matrix ";
      print "(multiple lines, <ctrl-D> = done):\n";
      eval { $new_matrix =
          Math::MatrixReal->new_from_string(join('',<STDIN>)); };
      if ($@)
      {
          $@ =~ s/\s+at\b.*?$//;
          print "${@}Please try again.\n";
      }
      else { last; }
  }

Possible error messages of the "new_from_string()" method are:

  Math::MatrixReal::new_from_string(): syntax error in input string
  Math::MatrixReal::new_from_string(): empty input string

If the input string has rows with varying numbers of columns,
the following warning will be printed to STDERR:

  Math::MatrixReal::new_from_string(): missing elements will be set to zero!

If everything is okay, the method returns an object reference to the
(newly allocated) matrix containing the elements you specified.

=item *

C<$new_matrix = $some_matrix-E<gt>shadow();>

Returns an object reference to a B<NEW> but B<EMPTY> matrix
(filled with zero's) of the B<SAME SIZE> as matrix "C<$some_matrix>".

Matrix "C<$some_matrix>" is not changed by this in any way.

=item *

C<$matrix1-E<gt>copy($matrix2);>

Copies the contents of matrix "C<$matrix2>" to an B<ALREADY EXISTING>
matrix "C<$matrix1>" (which must have the same size as matrix "C<$matrix2>"!).

Matrix "C<$matrix2>" is not changed by this in any way.

=item *

C<$twin_matrix = $some_matrix-E<gt>clone();>

Returns an object reference to a B<NEW> matrix of the B<SAME SIZE> as
matrix "C<$some_matrix>". The contents of matrix "C<$some_matrix>" have
B<ALREADY BEEN COPIED> to the new matrix "C<$twin_matrix>". This
is the method that the operator "=" is overloaded to when you type
C<$a = $b>, when C<$a> and C<$b> are matrices.

Matrix "C<$some_matrix>" is not changed by this in any way.


=head2 Matrix Row, Column and Element operations

=item *

C<$row_vector = $matrix-E<gt>row($row);>

This is a projection method which returns an object reference to
a B<NEW> matrix (which in fact is a (row) vector since it has only
one row) to which row number "C<$row>" of matrix "C<$matrix>" has
already been copied.

Matrix "C<$matrix>" is not changed by this in any way.

=item *

C<$column_vector = $matrix-E<gt>column($column);>

This is a projection method which returns an object reference to
a B<NEW> matrix (which in fact is a (column) vector since it has
only one column) to which column number "C<$column>" of matrix
"C<$matrix>" has already been copied.

Matrix "C<$matrix>" is not changed by this in any way.

=item *

C<$matrix-E<gt>assign($row,$column,$value);>

Explicitly assigns a value "C<$value>" to a single element of the
matrix "C<$matrix>", located in row "C<$row>" and column "C<$column>",
thereby replacing the value previously stored there.

=item *

C<$value = $matrix-E<gt>>C<element($row,$column);>

Returns the value of a specific element of the matrix "C<$matrix>",
located in row "C<$row>" and column "C<$column>".

=item *

C<$new_matrix = $matrix-E<gt>each( \&function )>;

Creates a new matrix by evaluating a code reference on each element of the 
given matrix. The function is passed the element, the row index and the column
index, in that order. The value the function returns ( or the value of the last
executed statement ) is the value given to the corresponding element in $new_matrix.

Example:

    # add 1 to every element in the matrix
    $matrix = $matrix->each ( sub { (shift) + 1 } );


Example:

    my $cofactor = $matrix->each( sub { my(undef,$i,$j) = @_;
        ($i+$j) % 2 == 0 ? $matrix->minor($i,$j)->det()
        : -1*$matrix->minor($i,$j)->det();
        } );

This code needs some explanation. For each element of $matrix, it throws away the actual value
and stores the row and column indexes in $i and $j. Then it sets element [$i,$j] in $cofactor
to the determinant of C<$matrix-E<gt>minor($i,$j)> if it is an "even" element, or C<-1*$matrix-E<gt>minor($i,$j)>
if it is an "odd" element.

=item *

C<$new_matrix = $matrix-E<gt>each_diag( \&function )>;

Creates a new matrix by evaluating a code reference on each diagonal element of the 
given matrix. The function is passed the element, the row index and the column
index, in that order. The value the function returns ( or the value of the last
executed statement ) is the value given to the corresponding element in $new_matrix.


=item * 

C<$matrix-E<gt>swap_col( $col1, $col2 );>

This method takes two one-based column numbers and swaps the values of each element in each column.
C<$matrix-E<gt>swap_col(2,3)> would replace column 2 in $matrix with column 3, and replace column
3 with column 2. 

=item * 

C<$matrix-E<gt>swap_row( $row1, $row2 );>

This method takes two one-based row numbers and swaps the values of each element in each row.
C<$matrix-E<gt>swap_row(2,3)> would replace row 2 in $matrix with row 3, and replace row
3 with row 2. 


=head2 Matrix Operations

=item *

C<$det = $matrix-E<gt>det();>

Returns the determinant of the matrix, without going through
the rigamarole of computing a LR decomposition. This method should
be much faster than LR decomposition if the matrix is diagonal or
triangular. Otherwise, it is just a wrapper for 
C<$matrix-E<gt>decompose_LR-E<gt>det_LR>. If the determinant is zero, 
there is no inverse and vice-versa. Only quadratic matrices have 
determinants.

=item *

C<$inverse = $matrix-E<gt>inverse();>

Returns the inverse of a matrix, without going through the
rigamarole of computing a LR decomposition. If no inverse exists,
undef is returned and an error is printed via C<carp()>.
This is nothing but a wrapper for C<$matrix-E<gt>decompose_LR-E<gt>invert_LR>.

=item *

C<($rows,$columns) = $matrix-E<gt>dim();>

Returns a list of two items, representing the number of rows
and columns the given matrix "C<$matrix>" contains.

=item *

C<$norm_one = $matrix-E<gt>norm_one();>

Returns the "one"-norm of the given matrix "C<$matrix>".

The "one"-norm is defined as follows:

For each column, the sum of the absolute values of the elements in the
different rows of that column is calculated. Finally, the maximum
of these sums is returned.

Note that the "one"-norm and the "maximum"-norm are mathematically
equivalent, although for the same matrix they usually yield a different
value.

Therefore, you should only compare values that have been calculated
using the same norm!

Throughout this package, the "one"-norm is (arbitrarily) used
for all comparisons, for the sake of uniformity and comparability,
except for the iterative methods "solve_GSM()", "solve_SSM()" and
"solve_RM()" which use either norm depending on the matrix itself.

=item *

C<$norm_max = $matrix-E<gt>norm_max();>

Returns the "maximum"-norm of the given matrix $matrix.

The "maximum"-norm is defined as follows:

For each row, the sum of the absolute values of the elements in the
different columns of that row is calculated. Finally, the maximum
of these sums is returned.

Note that the "maximum"-norm and the "one"-norm are mathematically
equivalent, although for the same matrix they usually yield a different
value.

Therefore, you should only compare values that have been calculated
using the same norm!

Throughout this package, the "one"-norm is (arbitrarily) used
for all comparisons, for the sake of uniformity and comparability,
except for the iterative methods "solve_GSM()", "solve_SSM()" and
"solve_RM()" which use either norm depending on the matrix itself.

=item *

C<$norm_sum = $matrix-E<gt>norm_sum();>

This is a very simple norm which is defined as the sum of the 
absolute values of every element.

=item *

C<$p_norm> = $matrix-E<gt>norm_p($n);>

This function returns the "p-norm" of a vector. The argument $n
must be a number greater than or equal to 1 or the string "Inf".
The p-norm is defined as (sum(x_i^p))^(1/p). In words, it raised
each element to the p-th power, adds them up, and then takes the
p-th root of that number. If the string "Inf" is passed, the
"infinity-norm" is computed, which is really the limit of the 
p-norm as p goes to infinity. It is defined as the maximum element
of the vector. Also, note that the familiar Euclidean distance 
between two vectors is just a special case of a p-norm, when p is
equal to 2.

Example:
    $a = Math::MatrixReal->new_from_cols([[1,2,3]]);
    $p1   = $a->norm_p(1);
        $p2   = $a->norm_p(2);    
        $p3   = $a->norm_p(3);    
    $pinf = $a->norm_p("Inf");

    print "(1,2,3,Inf) norm:\n$p1\n$p2\n$p3\n$pinf\n";

    $i1 = $a->new_from_rows([[1,0]]);
    $i2 = $a->new_from_rows([[0,1]]);

    # this should be sqrt(2) since it is the same as the 
    # hypotenuse of a 1 by 1 right triangle

    $dist  = ($i1-$i2)->norm_p(2);
    print "Distance is $dist, which should be " . sqrt(2) . "\n";

Output:

    (1,2,3,Inf) norm:
    6
    3.74165738677394139
    3.30192724889462668
    3

    Distance is 1.41421356237309505, which should be 1.41421356237309505



=item *

C<$frob_norm> = C<$matrix-E<gt>norm_frobenius();>

This norm is similar to that of a p-norm where p is 2, except it
acts on a B<matrix>, not a vector. Each element of the matrix is 
squared, this is added up, and then a square root is taken. 

=item *

C<$matrix-E<gt>spectral_radius();>

Returns the maximum value of the absolute value of all eigenvalues.
Currently this computes B<all> eigenvalues, then sifts through them
to find the largest in absolute value. Needless to say, this is very
inefficient, and in the future an algorithm that computes only the 
largest eigenvalue may be implemented.

=item *

C<$matrix1-E<gt>transpose($matrix2);>

Calculates the transposed matrix of matrix $matrix2 and stores
the result in matrix "C<$matrix1>" (which must already exist and have
the same size as matrix "C<$matrix2>"!).

This operation can also be carried out "in-place", i.e., input and
output matrix may be identical.

Transposition is a symmetry operation: imagine you rotate the matrix
along the axis of its main diagonal (going through elements (1,1),
(2,2), (3,3) and so on) by 180 degrees.

Another way of looking at it is to say that rows and columns are
swapped. In fact the contents of element C<(i,j)> are swapped
with those of element C<(j,i)>.

Note that (especially for vectors) it makes a big difference if you
have a row vector, like this:

  [ -1 0 1 ]

or a column vector, like this:

  [ -1 ]
  [  0 ]
  [  1 ]

the one vector being the transposed of the other!

This is especially true for the matrix product of two vectors:

               [ -1 ]
  [ -1 0 1 ] * [  0 ]  =  [ 2 ] ,  whereas
               [  1 ]

                             *     [ -1  0  1 ]
  [ -1 ]                                            [  1  0 -1 ]
  [  0 ] * [ -1 0 1 ]  =  [ -1 ]   [  1  0 -1 ]  =  [  0  0  0 ]
  [  1 ]                  [  0 ]   [  0  0  0 ]     [ -1  0  1 ]
                          [  1 ]   [ -1  0  1 ]

So be careful about what you really mean!

Hint: throughout this module, whenever a vector is explicitly required
for input, a B<COLUMN> vector is expected!

=item *

C<$trace = $matrix-E<gt>trace();>

This returns the trace of the matrix, which is defined as
the sum of the diagonal elements. The matrix must be
quadratic.

=item *

C<$minor = $matrix-E<gt>minor($row,$col);>

Returns the minor matrix corresponding to $row and $col. $matrix must be quadratic.
If $matrix is n rows by n cols, the minor of $row and $col will be an (n-1) by (n-1)
matrix. The minor is defined as crossing out the row and the col specified and returning
the remaining rows and columns as a matrix. This method is used by C<cofactor()>.

=item *

C<$cofactor = $matrix-E<gt>cofactor();>

The cofactor matrix is constructed as follows:

For each element, cross out the row and column that it sits in.
Now, take the determinant of the matrix that is left in the other
rows and columns.
Multiply the determinant by (-1)^(i+j), where i is the row index,
and j is the column index. 
Replace the given element with this value.

The cofactor matrix can be used to find the inverse of the matrix. One formula for the
inverse of a matrix is the cofactor matrix transposed divided by the original
determinant of the matrix. 

The following two inverses should be exactly the same:

    my $inverse1 = $matrix->inverse;
    my $inverse2 = ~($matrix->cofactor)->each( sub { (shift)/$matrix->det() } );

Caveat: Although the cofactor matrix is simple algorithm to compute the inverse of a matrix, and
can be used with pencil and paper for small matrices, it is comically slower than 
the native C<inverse()> function. Here is a small benchmark:

    # $matrix1 is 15x15
    $det = $matrix1->det;
    timethese( 10,
        {'inverse' => sub { $matrix1->inverse(); },
          'cofactor' => sub { (~$matrix1->cofactor)->each ( sub { (shift)/$det; } ) }
        } );


    Benchmark: timing 10 iterations of LR, cofactor, inverse...
        inverse:  1 wallclock secs ( 0.56 usr +  0.00 sys =  0.56 CPU) @ 17.86/s (n=10)
    cofactor: 36 wallclock secs (36.62 usr +  0.01 sys = 36.63 CPU) @  0.27/s (n=10)

=item *

C<$adjoint = $matrix-E<gt>adjoint();>

The adjoint is just the transpose of the cofactor matrix. This method is 
just an alias for C< ~($matrix-E<gt>cofactor)>.


=head2 Arithmetic Operations

=item *

C<$matrix1-E<gt>add($matrix2,$matrix3);>

Calculates the sum of matrix "C<$matrix2>" and matrix "C<$matrix3>"
and stores the result in matrix "C<$matrix1>" (which must already exist
and have the same size as matrix "C<$matrix2>" and matrix "C<$matrix3>"!).

This operation can also be carried out "in-place", i.e., the output and
one (or both) of the input matrices may be identical.

=item *

C<$matrix1-E<gt>subtract($matrix2,$matrix3);>

Calculates the difference of matrix "C<$matrix2>" minus matrix "C<$matrix3>"
and stores the result in matrix "C<$matrix1>" (which must already exist
and have the same size as matrix "C<$matrix2>" and matrix "C<$matrix3>"!).

This operation can also be carried out "in-place", i.e., the output and
one (or both) of the input matrices may be identical.

Note that this operation is the same as
C<$matrix1-E<gt>add($matrix2,-$matrix3);>, although the latter is
a little less efficient.

=item *

C<$matrix1-E<gt>multiply_scalar($matrix2,$scalar);>

Calculates the product of matrix "C<$matrix2>" and the number "C<$scalar>"
(i.e., multiplies each element of matrix "C<$matrix2>" with the factor
"C<$scalar>") and stores the result in matrix "C<$matrix1>" (which must
already exist and have the same size as matrix "C<$matrix2>"!).

This operation can also be carried out "in-place", i.e., input and
output matrix may be identical.

=item *

C<$product_matrix = $matrix1-E<gt>multiply($matrix2);>

Calculates the product of matrix "C<$matrix1>" and matrix "C<$matrix2>"
and returns an object reference to a new matrix "C<$product_matrix>" in
which the result of this operation has been stored.

Note that the dimensions of the two matrices "C<$matrix1>" and "C<$matrix2>"
(i.e., their numbers of rows and columns) must harmonize in the following
way (example):

                          [ 2 2 ]
                          [ 2 2 ]
                          [ 2 2 ]

              [ 1 1 1 ]   [ * * ]
              [ 1 1 1 ]   [ * * ]
              [ 1 1 1 ]   [ * * ]
              [ 1 1 1 ]   [ * * ]

I.e., the number of columns of matrix "C<$matrix1>" has to be the same
as the number of rows of matrix "C<$matrix2>".

The number of rows and columns of the resulting matrix "C<$product_matrix>"
is determined by the number of rows of matrix "C<$matrix1>" and the number
of columns of matrix "C<$matrix2>", respectively.

=item *

C<$matrix1-E<gt>negate($matrix2);>

Calculates the negative of matrix "C<$matrix2>" (i.e., multiplies
all elements with "-1") and stores the result in matrix "C<$matrix1>"
(which must already exist and have the same size as matrix "C<$matrix2>"!).

This operation can also be carried out "in-place", i.e., input and
output matrix may be identical.


=item *

C<$matrix_to_power = $matrix1-E<gt>exponent($integer);>

Raises the matrix to the C<$integer> power. Obviously, C<$integer> must
be an integer. If it is zero, the identity matrix is returned. If a negative
integer is given, the inverse will be computed (if it exists) and then raised
the the absolute value of C<$integer>. The matrix must be quadratic.


=head2 Boolean Matrix Operations

=item * 

C<$matrix-E<gt>is_quadratic();>

Returns a boolean value indicating if the given matrix is 
quadratic (also know as "square" or "n by n"). A matrix is 
quadratic if it has the same number of rows as it does columns.

=item * 

C<$matrix-E<gt>is_square();>

This is an alias for C<is_quadratic()>.


=item *

C<$matrix-E<gt>is_symmetric();>

Returns a boolean value indicating if the given matrix is
symmetric. By definition, a matrix is symmetric if and only
if (B<M>[I<i>,I<j>]=B<M>[I<j>,I<i>]). This is equivalent to
C<($matrix == ~$matrix)> but without memory allocation.
Only quadratic matrices can be symmetric.

Notes: A symmetric matrix always has real eigenvalues/eigenvectors.
A matrix plus its transpose is always symmetric.

=item *

C<$matrix-E<gt>is_skew_symmetric();>

Returns a boolean value indicating if the given matrix is
skew symmetric. By definition, a matrix is symmetric if and only
if (B<M>[I<i>,I<j>]=B<-M>[I<j>,I<i>]). This is equivalent to
C<($matrix == -(~$matrix))> but without memory allocation.
Only quadratic matrices can be skew symmetric.


=item *

C<$matrix-E<gt>is_diagonal();>

Returns a boolean value indicating if the given matrix is
diagonal, i.e. all of the nonzero elements are on the main diagonal.
Only quadratic matrices can be diagonal.

=item * 

C<$matrix-E<gt>is_tridiagonal();>

Returns a boolean value indicating if the given matrix is 
tridiagonal, i.e. all of the nonzero elements are on the main diagonal
or the diagonals above and below the main diagonal.
Only quadratic matrices can be tridiagonal.


=item *

C<$matrix-E<gt>is_upper_triangular();>

Returns a boolean value indicating if the given matrix is upper triangular, 
i.e. all of the nonzero elements not on the main diagonal are above it.
Only quadratic matrices can be upper triangular.
Note: diagonal matrices are both upper and lower triangular.


=item *

C<$matrix-E<gt>is_lower_triangular();>

Returns a boolean value indicating if the given matrix is lower triangular,
i.e. all of the nonzero elements not on the main diagonal are below it.
Only quadratic matrices can be lower triangular.
Note: diagonal matrices are both upper and lower triangular.

=item *

C<$matrix-E<gt>is_orthogonal();>

Returns a boolean value indicating if the given matrix is orthogonal.
An orthogonal matrix is has the property that the transpose equals the
inverse of the matrix. Instead of computing each and comparing them, this
method multiplies the matrix by it's transpose, and returns true if this 
turns out to be the identity matrix, false otherwise.
Only quadratic matrices can orthogonal.

=item *

C<$matrix-E<gt>is_binary();>

Returns a boolean value indicating if the given matrix is binary.
A matrix is binary if it contains only zeroes or ones. 

=item *

C<$matrix-E<gt>is_gramian();>

Returns a boolean value indicating if the give matrix is Gramian.
A matrix C<$A> is Gramian if and only if there exists a
square matrix C<$B> such that C<$A = ~$B*$B>. This is equivalent to
checking if C<$A> is symmetric and has all nonnegative eigenvalues, which
is what Math::MatrixReal uses to check for this property.

=item *

C<$matrix-E<gt>is_LR();>

Returns a boolean value indicating if the matrix is an LR decomposition
matrix.

=item *

C<$matrix-E<gt>is_positive();>

Returns a boolean value indicating if the matrix contains only
positive entries. Note that a zero entry is not positive and
will cause C<is_positive()> to return false.

=item *

C<$matrix-E<gt>is_negative();>

Returns a boolean value indicating if the matrix contains only
negative entries. Note that a zero entry is not negative and
will cause C<is_negative()> to return false.

=item *

C<$matrix-E<gt>is_periodic($k);>

Returns a boolean value indicating if the matrix is periodic
with period $k. This is true if C<$matrix ** ($k+1) == $matrix>.
When C<$k == 1>, this reduces down to the C<is_idempotent()>
function. 

=item *

C<$matrix-E<gt>is_idempotent();>

Returns a boolean value indicating if the matrix is idempotent,
which is defined as the square of the matrix being equal to 
the original matrix, i.e C<$matrix ** 2 == $matrix>.

=item *

C<$matrix-E<gt>is_row_vector();>

Returns a boolean value indicating if the matrix is a row vector.
A row vector is a matrix which is 1xn. Note that the 1x1 matrix is
both a row and column vector.

=item *

C<$matrix-E<gt>is_col_vector();>

Returns a boolean value indicating if the matrix is a col vector.
A col vector is a matrix which is nx1. Note that the 1x1 matrix is
both a row and column vector.

=head2 Eigensystems

=over 2

=item *

C<($l, $V) = $matrix-E<gt>sym_diagonalize();>

This method performs the diagonalization of the quadratic
I<symmetric> matrix B<M> stored in $matrix.
On output, B<l> is a column vector containing all the eigenvalues
of B<M> and B<V> is an orthogonal matrix which columns are the
corresponding normalized eigenvectors.
The primary property of an eigenvalue I<l> and an eigenvector
B<x> is of course that: B<M> * B<x> = I<l> * B<x>.

The method uses a Householder reduction to tridiagonal form
followed by a QL algoritm with implicit shifts on this
tridiagonal. (The tridiagonal matrix is kept internally
in a compact form in this routine to save memory.)
In fact, this routine wraps the householder() and
tri_diagonalize() methods described below when their
intermediate results are not desired.
The overall algorithmic complexity of this technique
is O(N^3). According to several books, the coefficient
hidden by the 'O' is one of the best possible for general
(symmetric) matrixes.

=item *

C<($T, $Q) = $matrix-E<gt>householder();>

This method performs the Householder algorithm which reduces
the I<n> by I<n> real I<symmetric> matrix B<M> contained
in $matrix to tridiagonal form.
On output, B<T> is a symmetric tridiagonal matrix (only
diagonal and off-diagonal elements are non-zero) and B<Q>
is an I<orthogonal> matrix performing the tranformation
between B<M> and B<T> (C<$M == $Q * $T * ~$Q>).

=item *

C<($l, $V) = $T-E<gt>tri_diagonalize([$Q]);>

This method diagonalizes the symmetric tridiagonal
matrix B<T>. On output, $l and $V are similar to the
output values described for sym_diagonalize().

The optional argument $Q corresponds to an orthogonal
transformation matrix B<Q> that should be used additionally
during B<V> (eigenvectors) computation. It should be supplied
if the desired eigenvectors correspond to a more general
symmetric matrix B<M> previously reduced by the
householder() method, not a mere tridiagonal. If B<T> is
really a tridiagonal matrix, B<Q> can be omitted (it
will be internally created in fact as an identity matrix).
The method uses a QL algorithm (with implicit shifts).

=item *

C<$l = $matrix-E<gt>sym_eigenvalues();>

This method computes the eigenvalues of the quadratic
I<symmetric> matrix B<M> stored in $matrix.
On output, B<l> is a column vector containing all the eigenvalues
of B<M>. Eigenvectors are not computed (on the contrary of
C<sym_diagonalize()>) and this method is more efficient
(even though it uses a similar algorithm with two phases).
However, understand that the algorithmic complexity of this
technique is still also O(N^3). But the coefficient hidden
by the 'O' is better by a factor of..., well, see your
benchmark, it's wiser.

This routine wraps the householder_tridiagonal() and
tri_eigenvalues() methods described below when the
intermediate tridiagonal matrix is not needed.

=item *

C<$T = $matrix-E<gt>householder_tridiagonal();>

This method performs the Householder algorithm which reduces
the I<n> by I<n> real I<symmetric> matrix B<M> contained
in $matrix to tridiagonal form.
On output, B<T> is the obtained symmetric tridiagonal matrix
(only diagonal and off-diagonal elements are non-zero). The
operation is similar to the householder() method, but potentially
a little more efficient as the transformation matrix is not
computed.

=item *

C<$l = $T-E<gt>tri_eigenvalues();>

This method computesthe eigenvalues of the symmetric
tridiagonal matrix B<T>. On output, $l is a vector
containing the eigenvalues (similar to C<sym_eigenvalues()>).
This method is much more efficient than tri_diagonalize()
when eigenvectors are not needed.

=back

=head2 Miscellaneous 

=item *

C<$matrix-E<gt>zero();>

Assigns a zero to every element of the matrix "C<$matrix>", i.e.,
erases all values previously stored there, thereby effectively
transforming the matrix into a "zero"-matrix or "null"-matrix,
the neutral element of the addition operation in a Ring.

(For instance the (quadratic) matrices with "n" rows and columns
and matrix addition and multiplication form a Ring. Most prominent
characteristic of a Ring is that multiplication is not commutative,
i.e., in general, "C<matrix1 * matrix2>" is not the same as
"C<matrix2 * matrix1>"!)

=item *

C<$matrix-E<gt>one();>

Assigns one's to the elements on the main diagonal (elements (1,1),
(2,2), (3,3) and so on) of matrix "C<$matrix>" and zero's to all others,
thereby erasing all values previously stored there and transforming the
matrix into a "one"-matrix, the neutral element of the multiplication
operation in a Ring.

(If the matrix is quadratic (which this method doesn't require, though),
then multiplying this matrix with itself yields this same matrix again,
and multiplying it with some other matrix leaves that other matrix
unchanged!)

=item *

C<$latex_string = $matrix-E<gt>as_latex( align=E<gt> "c", format =E<gt> "%s", name =E<gt> "" );>

This function returns the matrix as a LaTeX string. It takes a hash as an
argument which is used to control the style of the output. The hash element C<align>
may be "c","l" or "r", corresponding to center, left and right, respectively. The
C<format> element is a format string that is given to C<sprintf> to control the
style of number format, such a floating point or scientific notation. The C<name>
element can be used so that a LaTeX string of "$name = " is prepended to the string.

Example:

    my $a = Math::MatrixReal->new_from_cols([[ 1.234, 5.678, 9.1011],[1,2,3]] );
    print $a->as_latex( ( format => "%.2f", align => "l",name => "A" ) );

    Output:
    $A = $ $
    \left( \begin{array}{ll}
    1.23&1.00 \\
    5.68&2.00 \\
    9.10&3.00
    \end{array} \right)
    $

=item *

C<$yacas_string = $matrix-E<gt>as_yacas( format =E<gt> "%s", name =E<gt> "", semi =E<gt> 0 );>

This function returns the matrix as a string that can be read by Yacas.
It takes a hash as
an an argument which controls the style of the output. The
C<format> element is a format string that is given to C<sprintf> to control the
style of number format, such a floating point or scientific notation. The C<name>
element can be used so that "$name = " is prepended to the string. The <semi> element can
be set to 1 to that a semicolon is appended (so Matlab does not print out the matrix.) 

Example:

    $a = Math::MatrixReal->new_from_cols([[ 1.234, 5.678, 9.1011],[1,2,3]] );
    print $a->as_yacas( ( format => "%.2f", align => "l",name => "A" ) );

Output:

    A := {{1.23,1.00},{5.68,2.00},{9.10,3.00}}

=item *

C<$matlab_string = $matrix-E<gt>as_matlab( format =E<gt> "%s", name =E<gt> "", semi =E<gt> 0 );>

This function returns the matrix as a string that can be read by Matlab. It takes a hash as
an an argument which controls the style of the output. The
C<format> element is a format string that is given to C<sprintf> to control the
style of number format, such a floating point or scientific notation. The C<name>
element can be used so that "$name = " is prepended to the string. The <semi> element can
be set to 1 to that a semicolon is appended (so Matlab does not print out the matrix.) 

Example:

        my $a = Math::MatrixReal->new_from_rows([[ 1.234, 5.678, 9.1011],[1,2,3]] );
        print $a->as_matlab( ( format => "%.3f", name => "A",semi => 1 ) );

Output:
        A = [ 1.234 5.678 9.101;
         1.000 2.000 3.000];


=item *

C<$scilab_string = $matrix-E<gt>as_scilab( format =E<gt> "%s", name =E<gt> "", semi =E<gt> 0 );>

This function is just an alias for C<as_matlab()>, since both Scilab and Matlab have the
same matrix format.

=item *

C<$minimum = Math::MatrixReal::min($number1,$number2);>

Returns the minimum of the two numbers "C<number1>" and "C<number2>".

=item *

C<$maximum = Math::MatrixReal::max($number1,$number2);>

Returns the maximum of the two numbers "C<number1>" and "C<number2>".

=item *

C<$minimal_cost_matrix = $cost_matrix-E<gt>kleene();>

Copies the matrix "C<$cost_matrix>" (which has to be quadratic!) to
a new matrix of the same size (i.e., "clones" the input matrix) and
applies Kleene's algorithm to it.

See L<Math::Kleene(3)> for more details about this algorithm!

The method returns an object reference to the new matrix.

Matrix "C<$cost_matrix>" is not changed by this method in any way.

=item *

C<($norm_matrix,$norm_vector) = $matrix-E<gt>normalize($vector);>

This method is used to improve the numerical stability when solving
linear equation systems.

Suppose you have a matrix "A" and a vector "b" and you want to find
out a vector "x" so that C<A * x = b>, i.e., the vector "x" which
solves the equation system represented by the matrix "A" and the
vector "b".

Applying this method to the pair (A,b) yields a pair (A',b') where
each row has been divided by (the absolute value of) the greatest
coefficient appearing in that row. So this coefficient becomes equal
to "1" (or "-1") in the new pair (A',b') (all others become smaller
than one and greater than minus one).

Note that this operation does not change the equation system itself
because the same division is carried out on either side of the equation
sign!

The method requires a quadratic (!) matrix "C<$matrix>" and a vector
"C<$vector>" for input (the vector must be a column vector with the same
number of rows as the input matrix) and returns a list of two items
which are object references to a new matrix and a new vector, in this
order.

The output matrix and vector are clones of the input matrix and vector
to which the operation explained above has been applied.

The input matrix and vector are not changed by this in any way.

Example of how this method can affect the result of the methods to solve
equation systems (explained immediately below following this method):

Consider the following little program:

  #!perl -w

  use Math::MatrixReal qw(new_from_string);

  $A = Math::MatrixReal->new_from_string(<<"MATRIX");
  [  1   2   3  ]
  [  5   7  11  ]
  [ 23  19  13  ]
  MATRIX

  $b = Math::MatrixReal->new_from_string(<<"MATRIX");
  [   0   ]
  [   1   ]
  [  29   ]
  MATRIX

  $LR = $A->decompose_LR();
  if (($dim,$x,$B) = $LR->solve_LR($b))
  {
      $test = $A * $x;
      print "x = \n$x";
      print "A * x = \n$test";
  }

  ($A_,$b_) = $A->normalize($b);

  $LR = $A_->decompose_LR();
  if (($dim,$x,$B) = $LR->solve_LR($b_))
  {
      $test = $A * $x;
      print "x = \n$x";
      print "A * x = \n$test";
  }

This will print:

  x =
  [  1.000000000000E+00 ]
  [  1.000000000000E+00 ]
  [ -1.000000000000E+00 ]
  A * x =
  [  4.440892098501E-16 ]
  [  1.000000000000E+00 ]
  [  2.900000000000E+01 ]
  x =
  [  1.000000000000E+00 ]
  [  1.000000000000E+00 ]
  [ -1.000000000000E+00 ]
  A * x =
  [  0.000000000000E+00 ]
  [  1.000000000000E+00 ]
  [  2.900000000000E+01 ]

You can see that in the second example (where "normalize()" has been used),
the result is "better", i.e., more accurate!

=item *

C<$LR_matrix = $matrix-E<gt>decompose_LR();>

This method is needed to solve linear equation systems.

Suppose you have a matrix "A" and a vector "b" and you want to find
out a vector "x" so that C<A * x = b>, i.e., the vector "x" which
solves the equation system represented by the matrix "A" and the
vector "b".

You might also have a matrix "A" and a whole bunch of different
vectors "b1".."bk" for which you need to find vectors "x1".."xk"
so that C<A * xi = bi>, for C<i=1..k>.

Using Gaussian transformations (multiplying a row or column with
a factor, swapping two rows or two columns and adding a multiple
of one row or column to another), it is possible to decompose any
matrix "A" into two triangular matrices, called "L" and "R" (for
"Left" and "Right").

"L" has one's on the main diagonal (the elements (1,1), (2,2), (3,3)
and so so), non-zero values to the left and below of the main diagonal
and all zero's in the upper right half of the matrix.

"R" has non-zero values on the main diagonal as well as to the right
and above of the main diagonal and all zero's in the lower left half
of the matrix, as follows:

          [ 1 0 0 0 0 ]      [ x x x x x ]
          [ x 1 0 0 0 ]      [ 0 x x x x ]
      L = [ x x 1 0 0 ]  R = [ 0 0 x x x ]
          [ x x x 1 0 ]      [ 0 0 0 x x ]
          [ x x x x 1 ]      [ 0 0 0 0 x ]

Note that "C<L * R>" is equivalent to matrix "A" in the sense that
C<L * R * x = b  E<lt>==E<gt>  A * x = b> for all vectors "x", leaving
out of account permutations of the rows and columns (these are taken
care of "magically" by this module!) and numerical errors.

Trick:

Because we know that "L" has one's on its main diagonal, we can
store both matrices together in the same array without information
loss! I.e.,

                 [ R R R R R ]
                 [ L R R R R ]
            LR = [ L L R R R ]
                 [ L L L R R ]
                 [ L L L L R ]

Beware, though, that "LR" and "C<L * R>" are not the same!!!

Note also that for the same reason, you cannot apply the method "normalize()"
to an "LR" decomposition matrix. Trying to do so will yield meaningless
rubbish!

(You need to apply "normalize()" to each pair (Ai,bi) B<BEFORE> decomposing
the matrix "Ai'"!)

Now what does all this help us in solving linear equation systems?

It helps us because a triangular matrix is the next best thing
that can happen to us besides a diagonal matrix (a matrix that
has non-zero values only on its main diagonal - in which case
the solution is trivial, simply divide "C<b[i]>" by "C<A[i,i]>"
to get "C<x[i]>"!).

To find the solution to our problem "C<A * x = b>", we divide this
problem in parts: instead of solving C<A * x = b> directly, we first
decompose "A" into "L" and "R" and then solve "C<L * y = b>" and
finally "C<R * x = y>" (motto: divide and rule!).

From the illustration above it is clear that solving "C<L * y = b>"
and "C<R * x = y>" is straightforward: we immediately know that
C<y[1] = b[1]>. We then deduce swiftly that

  y[2] = b[2] - L[2,1] * y[1]

(and we know "C<y[1]>" by now!), that

  y[3] = b[3] - L[3,1] * y[1] - L[3,2] * y[2]

and so on.

Having effortlessly calculated the vector "y", we now proceed to
calculate the vector "x" in a similar fashion: we see immediately
that C<x[n] = y[n] / R[n,n]>. It follows that

  x[n-1] = ( y[n-1] - R[n-1,n] * x[n] ) / R[n-1,n-1]

and

  x[n-2] = ( y[n-2] - R[n-2,n-1] * x[n-1] - R[n-2,n] * x[n] )
           / R[n-2,n-2]

and so on.

You can see that - especially when you have many vectors "b1".."bk"
for which you are searching solutions to C<A * xi = bi> - this scheme
is much more efficient than a straightforward, "brute force" approach.

This method requires a quadratic matrix as its input matrix.

If you don't have that many equations, fill up with zero's (i.e., do
nothing to fill the superfluous rows if it's a "fresh" matrix, i.e.,
a matrix that has been created with "new()" or "shadow()").

The method returns an object reference to a new matrix containing the
matrices "L" and "R".

The input matrix is not changed by this method in any way.

Note that you can "copy()" or "clone()" the result of this method without
losing its "magical" properties (for instance concerning the hidden
permutations of its rows and columns).

However, as soon as you are applying any method that alters the contents
of the matrix, its "magical" properties are stripped off, and the matrix
immediately reverts to an "ordinary" matrix (with the values it just happens
to contain at that moment, be they meaningful as an ordinary matrix or not!).

=item *

C<($dimension,$x_vector,$base_matrix) = $LR_matrix>C<-E<gt>>C<solve_LR($b_vector);>

Use this method to actually solve an equation system.

Matrix "C<$LR_matrix>" must be a (quadratic) matrix returned by the
method "decompose_LR()", the LR decomposition matrix of the matrix
"A" of your equation system C<A * x = b>.

The input vector "C<$b_vector>" is the vector "b" in your equation system
C<A * x = b>, which must be a column vector and have the same number of
rows as the input matrix "C<$LR_matrix>".

The method returns a list of three items if a solution exists or an
empty list otherwise (!).

Therefore, you should always use this method like this:

  if ( ($dim,$x_vec,$base) = $LR->solve_LR($b_vec) )
  {
      # do something with the solution...
  }
  else
  {
      # do something with the fact that there is no solution...
  }

The three items returned are: the dimension "C<$dimension>" of the solution
space (which is zero if only one solution exists, one if the solution is
a straight line, two if the solution is a plane, and so on), the solution
vector "C<$x_vector>" (which is the vector "x" of your equation system
C<A * x = b>) and a matrix "C<$base_matrix>" representing a base of the
solution space (a set of vectors which put up the solution space like
the spokes of an umbrella).

Only the first "C<$dimension>" columns of this base matrix actually
contain entries, the remaining columns are all zero.

Now what is all this stuff with that "base" good for?

The output vector "x" is B<ALWAYS> a solution of your equation system
C<A * x = b>.

But also any vector "C<$vector>"

  $vector = $x_vector->clone();

  $machine_infinity = 1E+99; # or something like that

  for ( $i = 1; $i <= $dimension; $i++ )
  {
      $vector += rand($machine_infinity) * $base_matrix->column($i);
  }

is a solution to your problem C<A * x = b>, i.e., if "C<$A_matrix>" contains
your matrix "A", then

  print abs( $A_matrix * $vector - $b_vector ), "\n";

should print a number around 1E-16 or so!

By the way, note that you can actually calculate those vectors "C<$vector>"
a little more efficient as follows:

  $rand_vector = $x_vector->shadow();

  $machine_infinity = 1E+99; # or something like that

  for ( $i = 1; $i <= $dimension; $i++ )
  {
      $rand_vector->assign($i,1, rand($machine_infinity) );
  }

  $vector = $x_vector + ( $base_matrix * $rand_vector );

Note that the input matrix and vector are not changed by this method
in any way.

=item *

C<$inverse_matrix = $LR_matrix-E<gt>invert_LR();>

Use this method to calculate the inverse of a given matrix "C<$LR_matrix>",
which must be a (quadratic) matrix returned by the method "decompose_LR()".

The method returns an object reference to a new matrix of the same size as
the input matrix containing the inverse of the matrix that you initially
fed into "decompose_LR()" B<IF THE INVERSE EXISTS>, or an empty list
otherwise.

Therefore, you should always use this method in the following way:

  if ( $inverse_matrix = $LR->invert_LR() )
  {
      # do something with the inverse matrix...
  }
  else
  {
      # do something with the fact that there is no inverse matrix...
  }

Note that by definition (disregarding numerical errors), the product
of the initial matrix and its inverse (or vice-versa) is always a matrix
containing one's on the main diagonal (elements (1,1), (2,2), (3,3) and
so on) and zero's elsewhere.

The input matrix is not changed by this method in any way.

=item *

C<$condition = $matrix-E<gt>condition($inverse_matrix);>

In fact this method is just a shortcut for

  abs($matrix) * abs($inverse_matrix)

Both input matrices must be quadratic and have the same size, and the result
is meaningful only if one of them is the inverse of the other (for instance,
as returned by the method "invert_LR()").

The number returned is a measure of the "condition" of the given matrix
"C<$matrix>", i.e., a measure of the numerical stability of the matrix.

This number is always positive, and the smaller its value, the better the
condition of the matrix (the better the stability of all subsequent
computations carried out using this matrix).

Numerical stability means for example that if

  abs( $vec_correct - $vec_with_error ) < $epsilon

holds, there must be a "C<$delta>" which doesn't depend on the vector
"C<$vec_correct>" (nor "C<$vec_with_error>", by the way) so that

  abs( $matrix * $vec_correct - $matrix * $vec_with_error ) < $delta

also holds.

=item *

C<$determinant = $LR_matrix-E<gt>det_LR();>

Calculates the determinant of a matrix, whose LR decomposition matrix
"C<$LR_matrix>" must be given (which must be a (quadratic) matrix
returned by the method "decompose_LR()").

In fact the determinant is a by-product of the LR decomposition: It is
(in principle, that is, except for the sign) simply the product of the
elements on the main diagonal (elements (1,1), (2,2), (3,3) and so on)
of the LR decomposition matrix.

(The sign is taken care of "magically" by this module)

=item *

C<$order = $LR_matrix-E<gt>order_LR();>

Calculates the order (called "Rang" in German) of a matrix, whose
LR decomposition matrix "C<$LR_matrix>" must be given (which must
be a (quadratic) matrix returned by the method "decompose_LR()").

This number is a measure of the number of linear independent row
and column vectors (= number of linear independent equations in
the case of a matrix representing an equation system) of the
matrix that was initially fed into "decompose_LR()".

If "n" is the number of rows and columns of the (quadratic!) matrix,
then "n - order" is the dimension of the solution space of the
associated equation system.

=item *

C<$rank = $LR_matrix-E<gt>rank_LR();>

This is an alias for the C<order_LR()> function. The "order"
is usually called the "rank" in the United States.

=item *

C<$scalar_product = $vector1-E<gt>scalar_product($vector2);>

Returns the scalar product of vector "C<$vector1>" and vector "C<$vector2>".

Both vectors must be column vectors (i.e., a matrix having
several rows but only one column).

This is a (more efficient!) shortcut for

  $temp           = ~$vector1 * $vector2;
  $scalar_product =  $temp->element(1,1);

or the sum C<i=1..n> of the products C<vector1[i] * vector2[i]>.

Provided none of the two input vectors is the null vector, then
the two vectors are orthogonal, i.e., have an angle of 90 degrees
between them, exactly when their scalar product is zero, and
vice-versa.

=item *

C<$vector_product = $vector1-E<gt>vector_product($vector2);>

Returns the vector product of vector "C<$vector1>" and vector "C<$vector2>".

Both vectors must be column vectors (i.e., a matrix having several rows
but only one column).

Currently, the vector product is only defined for 3 dimensions (i.e.,
vectors with 3 rows); all other vectors trigger an error message.

In 3 dimensions, the vector product of two vectors "x" and "y"
is defined as

              |  x[1]  y[1]  e[1]  |
  determinant |  x[2]  y[2]  e[2]  |
              |  x[3]  y[3]  e[3]  |

where the "C<x[i]>" and "C<y[i]>" are the components of the two vectors
"x" and "y", respectively, and the "C<e[i]>" are unity vectors (i.e.,
vectors with a length equal to one) with a one in row "i" and zero's
elsewhere (this means that you have numbers and vectors as elements
in this matrix!).

This determinant evaluates to the rather simple formula

  z[1] = x[2] * y[3] - x[3] * y[2]
  z[2] = x[3] * y[1] - x[1] * y[3]
  z[3] = x[1] * y[2] - x[2] * y[1]

A characteristic property of the vector product is that the resulting
vector is orthogonal to both of the input vectors (if neither of both
is the null vector, otherwise this is trivial), i.e., the scalar product
of each of the input vectors with the resulting vector is always zero.

=item *

C<$length = $vector-E<gt>length();>

This is actually a shortcut for

  $length = sqrt( $vector->scalar_product($vector) );

and returns the length of a given (column!) vector "C<$vector>".

Note that the "length" calculated by this method is in fact the
"two"-norm of a vector "C<$vector>"!

The general definition for norms of vectors is the following:

  sub vector_norm
  {
      croak "Usage: \$norm = \$vector->vector_norm(\$n);"
        if (@_ != 2);

      my($vector,$n) = @_;
      my($rows,$cols) = ($vector->[1],$vector->[2]);
      my($k,$comp,$sum);

      croak "Math::MatrixReal::vector_norm(): vector is not a column vector"
        unless ($cols == 1);

      croak "Math::MatrixReal::vector_norm(): norm index must be > 0"
        unless ($n > 0);

      croak "Math::MatrixReal::vector_norm(): norm index must be integer"
        unless ($n == int($n));

      $sum = 0;
      for ( $k = 0; $k < $rows; $k++ )
      {
          $comp = abs( $vector->[0][$k][0] );
          $sum += $comp ** $n;
      }
      return( $sum ** (1 / $n) );
  }

Note that the case "n = 1" is the "one"-norm for matrices applied to a
vector, the case "n = 2" is the euclidian norm or length of a vector,
and if "n" goes to infinity, you have the "infinity"- or "maximum"-norm
for matrices applied to a vector!

=item *

C<$xn_vector = $matrix-E<gt>>C<solve_GSM($x0_vector,$b_vector,$epsilon);>

=item *

C<$xn_vector = $matrix-E<gt>>C<solve_SSM($x0_vector,$b_vector,$epsilon);>

=item *

C<$xn_vector = $matrix-E<gt>>C<solve_RM($x0_vector,$b_vector,$weight,$epsilon);>

In some cases it might not be practical or desirable to solve an
equation system "C<A * x = b>" using an analytical algorithm like
the "decompose_LR()" and "solve_LR()" method pair.

In fact in some cases, due to the numerical properties (the "condition")
of the matrix "A", the numerical error of the obtained result can be
greater than by using an approximative (iterative) algorithm like one
of the three implemented here.

All three methods, GSM ("Global Step Method" or "Gesamtschrittverfahren"),
SSM ("Single Step Method" or "Einzelschrittverfahren") and RM ("Relaxation
Method" or "Relaxationsverfahren"), are fix-point iterations, that is, can
be described by an iteration function "C<x(t+1) = Phi( x(t) )>" which has
the property:

  Phi(x)  =  x    <==>    A * x  =  b

We can define "C<Phi(x)>" as follows:

  Phi(x)  :=  ( En - A ) * x  +  b

where "En" is a matrix of the same size as "A" ("n" rows and columns)
with one's on its main diagonal and zero's elsewhere.

This function has the required property.

Proof:

           A * x        =   b

  <==>  -( A * x )      =  -b

  <==>  -( A * x ) + x  =  -b + x

  <==>  -( A * x ) + x + b  =  x

  <==>  x - ( A * x ) + b  =  x

  <==>  ( En - A ) * x + b  =  x

This last step is true because

  x[i] - ( a[i,1] x[1] + ... + a[i,i] x[i] + ... + a[i,n] x[n] ) + b[i]

is the same as

  ( -a[i,1] x[1] + ... + (1 - a[i,i]) x[i] + ... + -a[i,n] x[n] ) + b[i]

qed

Note that actually solving the equation system "C<A * x = b>" means
to calculate

        a[i,1] x[1] + ... + a[i,i] x[i] + ... + a[i,n] x[n]  =  b[i]

  <==>  a[i,i] x[i]  =
        b[i]
        - ( a[i,1] x[1] + ... + a[i,i] x[i] + ... + a[i,n] x[n] )
        + a[i,i] x[i]

  <==>  x[i]  =
        ( b[i]
            - ( a[i,1] x[1] + ... + a[i,i] x[i] + ... + a[i,n] x[n] )
            + a[i,i] x[i]
        ) / a[i,i]

  <==>  x[i]  =
        ( b[i] -
            ( a[i,1] x[1] + ... + a[i,i-1] x[i-1] +
              a[i,i+1] x[i+1] + ... + a[i,n] x[n] )
        ) / a[i,i]

There is one major restriction, though: a fix-point iteration is
guaranteed to converge only if the first derivative of the iteration
function has an absolute value less than one in an area around the
point "C<x(*)>" for which "C<Phi( x(*) ) = x(*)>" is to be true, and
if the start vector "C<x(0)>" lies within that area!

This is best verified grafically, which unfortunately is impossible
to do in this textual documentation!

See literature on Numerical Analysis for details!

In our case, this restriction translates to the following three conditions:

There must exist a norm so that the norm of the matrix of the iteration
function, C<( En - A )>, has a value less than one, the matrix "A" may
not have any zero value on its main diagonal and the initial vector
"C<x(0)>" must be "good enough", i.e., "close enough" to the solution
"C<x(*)>".

(Remember school math: the first derivative of a straight line given by
"C<y = a * x + b>" is "a"!)

The three methods expect a (quadratic!) matrix "C<$matrix>" as their
first argument, a start vector "C<$x0_vector>", a vector "C<$b_vector>"
(which is the vector "b" in your equation system "C<A * x = b>"), in the
case of the "Relaxation Method" ("RM"), a real number "C<$weight>" best
between zero and two, and finally an error limit (real number) "C<$epsilon>".

(Note that the weight "C<$weight>" used by the "Relaxation Method" ("RM")
is B<NOT> checked to lie within any reasonable range!)

The three methods first test the first two conditions of the three
conditions listed above and return an empty list if these conditions
are not fulfilled.

Therefore, you should always test their return value using some
code like:

  if ( $xn_vector = $A_matrix->solve_GSM($x0_vector,$b_vector,1E-12) )
  {
      # do something with the solution...
  }
  else
  {
      # do something with the fact that there is no solution...
  }

Otherwise, they iterate until C<abs( Phi(x) - x ) E<lt> epsilon>.

(Beware that theoretically, infinite loops might result if the starting
vector is too far "off" the solution! In practice, this shouldn't be
a problem. Anyway, you can always press <ctrl-C> if you think that the
iteration takes too long!)

The difference between the three methods is the following:

In the "Global Step Method" ("GSM"), the new vector "C<x(t+1)>"
(called "y" here) is calculated from the vector "C<x(t)>"
(called "x" here) according to the formula:

  y[i] =
  ( b[i]
      - ( a[i,1] x[1] + ... + a[i,i-1] x[i-1] +
          a[i,i+1] x[i+1] + ... + a[i,n] x[n] )
  ) / a[i,i]

In the "Single Step Method" ("SSM"), the components of the vector
"C<x(t+1)>" which have already been calculated are used to calculate
the remaining components, i.e.

  y[i] =
  ( b[i]
      - ( a[i,1] y[1] + ... + a[i,i-1] y[i-1] +  # note the "y[]"!
          a[i,i+1] x[i+1] + ... + a[i,n] x[n] )  # note the "x[]"!
  ) / a[i,i]

In the "Relaxation method" ("RM"), the components of the vector
"C<x(t+1)>" are calculated by "mixing" old and new value (like
cold and hot water), and the weight "C<$weight>" determines the
"aperture" of both the "hot water tap" as well as of the "cold
water tap", according to the formula:

  y[i] =
  ( b[i]
      - ( a[i,1] y[1] + ... + a[i,i-1] y[i-1] +  # note the "y[]"!
          a[i,i+1] x[i+1] + ... + a[i,n] x[n] )  # note the "x[]"!
  ) / a[i,i]
  y[i] = weight * y[i] + (1 - weight) * x[i]

Note that the weight "C<$weight>" should be greater than zero and
less than two (!).

The three methods are supposed to be of different efficiency.
Experiment!

Remember that in most cases, it is probably advantageous to first
"normalize()" your equation system prior to solving it!

=head1 OVERLOADED OPERATORS

=head2 SYNOPSIS

=over 2

=item *

Unary operators:

"C<->", "C<~>", "C<abs>", C<test>, "C<!>", 'C<"">'

=item *

Binary (arithmetic) operators:

"C<+>", "C<->", "C<*>", "C<**>",
"C<+=>", "C<-=>", "C<*=>", "C<**=>"

=item *

Binary (relational) operators:

"C<==>", "C<!=>", "C<E<lt>>", "C<E<lt>=>", "C<E<gt>>", "C<E<gt>=>"

"C<eq>", "C<ne>", "C<lt>", "C<le>", "C<gt>", "C<ge>"

Note that the latter ("C<eq>", "C<ne>", ... ) are just synonyms
of the former ("C<==>", "C<!=>", ... ), defined for convenience
only.

=back

=head2 DESCRIPTION

=over 5

=item '-'

Unary minus

Returns the negative of the given matrix, i.e., the matrix with
all elements multiplied with the factor "-1".

Example:

    $matrix = -$matrix;

=item '~'

Transposition

Returns the transposed of the given matrix.

Examples:

    $temp = ~$vector * $vector;
    $length = sqrt( $temp->element(1,1) );

    if (~$matrix == $matrix) { # matrix is symmetric ... }

=item abs

Norm

Returns the "one"-Norm of the given matrix.

Example:

    $error = abs( $A * $x - $b );

=item test

Boolean test

Tests wether there is at least one non-zero element in the matrix.

Example:

    if ($xn_vector) { # result of iteration is not zero ... }

=item '!'

Negated boolean test

Tests wether the matrix contains only zero's.

Examples:

    if (! $b_vector) { # heterogenous equation system ... }
    else             { # homogenous equation system ... }

    unless ($x_vector) { # not the null-vector! }

=item '""""'

"Stringify" operator

Converts the given matrix into a string.

Uses scientific representation to keep precision loss to a minimum in case
you want to read this string back in again later with "new_from_string()".

Uses a 13-digit mantissa and a 20-character field for each element so that
lines will wrap nicely on an 80-column screen.

Examples:

    $matrix = Math::MatrixReal->new_from_string(<<"MATRIX");
    [ 1  0 ]
    [ 0 -1 ]
    MATRIX
    print "$matrix";

    [  1.000000000000E+00  0.000000000000E+00 ]
    [  0.000000000000E+00 -1.000000000000E+00 ]

    $string = "$matrix";
    $test = Math::MatrixReal->new_from_string($string);
    if ($test == $matrix) { print ":-)\n"; } else { print ":-(\n"; }

=item '+'

Addition

Returns the sum of the two given matrices.

Examples:

    $matrix_S = $matrix_A + $matrix_B;

    $matrix_A += $matrix_B;

=item '-'

Subtraction

Returns the difference of the two given matrices.

Examples:

    $matrix_D = $matrix_A - $matrix_B;

    $matrix_A -= $matrix_B;

Note that this is the same as:

    $matrix_S = $matrix_A + -$matrix_B;

    $matrix_A += -$matrix_B;

(The latter are less efficient, though)

=item '*'

Multiplication

Returns the matrix product of the two given matrices or
the product of the given matrix and scalar factor.

Examples:

    $matrix_P = $matrix_A * $matrix_B;

    $matrix_A *= $matrix_B;

    $vector_b = $matrix_A * $vector_x;

    $matrix_B = -1 * $matrix_A;

    $matrix_B = $matrix_A * -1;

    $matrix_A *= -1;

=item '**'

Exponentiation

Returns the matrix raised to an integer power. If 0 is passed,
the identity matrix is returned. If a negative integer is passed,
it computes the inverse (if it exists) and then raised the inverse
to the absolute value of the integer. The matrix must be quadratic.

Examples:

    $matrix2 = $matrix ** 2;

    $matrix **= 2;

    $inv2 = $matrix ** -2;

    $ident = $matrix ** 0;



=item '=='

Equality

Tests two matrices for equality.

Example:

    if ( $A * $x == $b ) { print "EUREKA!\n"; }

Note that in most cases, due to numerical errors (due to the finite
precision of computer arithmetics), it is a bad idea to compare two
matrices or vectors this way.

Better use the norm of the difference of the two matrices you want
to compare and compare that norm with a small number, like this:

    if ( abs( $A * $x - $b ) < 1E-12 ) { print "BINGO!\n"; }

=item '!='

Inequality

Tests two matrices for inequality.

Example:

    while ($x0_vector != $xn_vector) { # proceed with iteration ... }

(Stops when the iteration becomes stationary)

Note that (just like with the '==' operator), it is usually a bad idea
to compare matrices or vectors this way. Compare the norm of the difference
of the two matrices with a small number instead.

=item 'E<lt>'

Less than

Examples:

    if ( $matrix1 < $matrix2 ) { # ... }

    if ( $vector < $epsilon ) { # ... }

    if ( 1E-12 < $vector ) { # ... }

    if ( $A * $x - $b < 1E-12 ) { # ... }

These are just shortcuts for saying:

    if ( abs($matrix1) < abs($matrix2) ) { # ... }

    if ( abs($vector) < abs($epsilon) ) { # ... }

    if ( abs(1E-12) < abs($vector) ) { # ... }

    if ( abs( $A * $x - $b ) < abs(1E-12) ) { # ... }

Uses the "one"-norm for matrices and Perl's built-in "abs()" for scalars.

=item 'E<lt>='

Less than or equal

As with the '<' operator, this is just a shortcut for the same expression
with "abs()" around all arguments.

Example:

    if ( $A * $x - $b <= 1E-12 ) { # ... }

which in fact is the same as:

    if ( abs( $A * $x - $b ) <= abs(1E-12) ) { # ... }

Uses the "one"-norm for matrices and Perl's built-in "abs()" for scalars.

=item 'E<gt>'

Greater than

As with the '<' and '<=' operator, this

    if ( $xn - $x0 > 1E-12 ) { # ... }

is just a shortcut for:

    if ( abs( $xn - $x0 ) > abs(1E-12) ) { # ... }

Uses the "one"-norm for matrices and Perl's built-in "abs()" for scalars.

=item 'E<gt>='

Greater than or equal

As with the '<', '<=' and '>' operator, the following

    if ( $LR >= $A ) { # ... }

is simply a shortcut for:

    if ( abs($LR) >= abs($A) ) { # ... }

Uses the "one"-norm for matrices and Perl's built-in "abs()" for scalars.

=back

=head1 SEE ALSO

Math::VectorReal, Math::PARI, Math::MatrixBool,
DFA::Kleene, Math::Kleene,
Set::IntegerRange, Set::IntegerFast .

=head1 VERSION

This man page documents Math::MatrixReal version 2.01.

The latest version can always be found at
http://leto.net/code/Math-MatrixReal/

=head1 AUTHORS

Steffen Beyer <sb@engelschall.com>, Rodolphe Ortalo <ortalo@laas.fr>,
Jonathan Leto <jonathan@leto.net>.

Currently maintained by Jonathan Leto, send all bugs/patches to me.

=head1 CREDITS

Many thanks to Prof. Pahlings for stoking the fire of my enthusiasm for
Algebra and Linear Algebra at the university (RWTH Aachen, Germany), and
to Prof. Esser and his assistant, Mr. Jarausch, for their fascinating
lectures in Numerical Analysis!

=head1 COPYRIGHT

Copyright (c) 1996-2002 by Steffen Beyer, Rodolphe Ortalo, Jonathan Leto.
All rights reserved.

=head1 LICENSE AGREEMENT

This package is free software; you can redistribute it and/or
modify it under the same terms as Perl itself.

