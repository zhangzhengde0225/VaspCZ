#!/usr/bin/env perl
#;-*- Perl -*-

# todo:
#  incorporate FPIchainfix
#  extend .con format: allow comment fields to give radius + color
#  add numbers to scale
#  allow multiple views on same page?
#  fix labeling so it doesn't change when new atoms added?
#    or autolabel them? or don't add as #1?
#   comm. doesn't work after 'save sel?'
# make labels smaller

#require 'ctime.pl';

$M_PI = 3.14159265358979323846264338327950288;
$sphflag=0;

$progname = 'tcon';
$version = '0.40';
$rcfile = '.tconrc';
$datfile = '.tcondat';

$nd = 9;       # block size for atom  (see %aref)
$defrad=1.0;   # default atomic radius
$defcol=0.333; # default atomic color
$defsphere=0;  # default sphere toggle

%aref = ('x',0,'y',1,'z',2,'tag',3,'type',4,'radius',5,'color',6,'label',7,'sphere',8);

# Added APS June 2, 1994:
$nbd = 9;       # block size for bonds (see %bref)
$maxlength = 2.0; # default maximum bond length
%bref = ('atom1',0,'atom2',1,'length',2,'label',3,'width',4,'shade',5,'red',6,'blue',7,'green',8);
# End of addition

# %con,@box,$nt,$natm,@typen,@atm,@velocs,@rgbs,@comments  <~ current conf vars

# %con: {name,seed,time,shear,pmass,vol,nperm,nharm,masses}

sub bound_vals {
  local($min,$max,@ar) = @_;
  foreach $i (@ar) {$i = $min if $i<$min; $i = $max if $i>$max;}
  wantarray?@ar:$ar[0];
}
sub min {local($m)= 1e32; foreach $i (@_) {$m = $i if $m > $i} $m}
sub max {local($m)=-1e32; foreach $i (@_) {$m = $i if $m < $i} $m}
sub abs {
  local(@foo)=@_;
  foreach $i (@foo) {$i=-$i if $i<0}
  wantarray?@foo:$foo[0];
}
sub getnums {  #todo: make more efficient
  local($n,$i,@ar) = @_;
  for($i=0;$i<$n;$i++) {
    $ar[$i]=$1,s/[^\d.\+-]*$1//
#   ,print STDOUT "${i}th match was _$ar[$i]_\n"
     if /^[^\d.\+-]*([\dEe.\+-]+)/;
  }
  wantarray?@ar:$ar[0];
}
sub arr_to_tbl {     # converts list of atoms to boolean table of belonging
  local($n,@ar,@tbl) = @_;
  $#tbl = $n-1;
  foreach $i (@ar) {$tbl[$i]=1;}
  return @tbl;
}
sub group_reduce    {local(@tbl); grep(!$tbl[$_]++,@_)}
sub group_intersect {local(@tbl); grep( $tbl[$_]++,@_)}
sub group_subtract {
  local($n,@a,@b,@tbl) = @_;  @b = splice(@a,$n);
  @tbl = &arr_to_tbl($natm,@b);
  grep(!$tbl[$_],@a);
}
sub group_xor {local($_,@tbl); grep($tbl[$_]++,@_); grep($tbl[$_]==1,@_);}

sub get_input {
  return shift(@FN_ARGS) if $#FN_ARGS >= 0;
  local($line);
  print $_[0];
  chop($line = <STDIN>);
  print "\n";
  return $line;
}
sub _store_prep {
  local($a,$b,$c,@r,$s,$t) = @_;
  foreach $i (@r) {
    $i =~ s/'/\\'/g if $a;
    $i =~ s/"/\\"/g if $b;
    $i =~ s/,/\\,/g if $c;
    $s .= ',' if $t++;
    $s .= $a;
  }
  $s;
}
sub store_config {
  local($_,$name,$s,$a,$b,@foo);
  $_ = &get_input(" Enter name to store configuration under: ");
  (print "Store aborted.\n"),return if /^\s*$/;
  $name = $_;
  $s = '%con = (';
  foreach $i (keys %con) {
    $a = $con{$i};
    $a =~ s/'/\\'/g;
    $s .= ',' if $b++;
    $s .= "'$i','$a'";
  }
  $s .= ");\n\@box = (" . join(',',@box) . ");\n";
  $s .= "\$nt=$nt; \$natm=$natm; \@typen=(" . join(',',@typen) . ");\n";
  @foo = @atm;
  foreach $i (0..$natm-1) { $foo[$i*$nd+$aref{'label'}] =~ s/'/\\'/g }
  $s .= '@atm=(' . join(',',@foo) . ");\n";
  $s .= '@velocs=(' . join(',',@velocs) . ");\n";
  $s .= '@rgbs=(' . join(',',@rgbs) . ");\n";
  $a = $b = 0;
  $s .= '@comments=(';
  foreach $i (@comments) {
    $a = $i;
    $a =~ s/'/\\'/g;
    $s .= ',' if $b++;
    $s .= "'$a'";
  }
  $s .= ");\n";
  $savedcons{$name} = $s;
  print "Stored as '$name'.\n";
#  print "YOU: Saved is \n--------------------------\n$s\n-------------\n";
}
sub restore_config {
  local($j);
  while(1) {
    $_ = &get_input(" Enter name configuration was stored as (or <CR> to list): ");
    last unless /^\s*$/;
    foreach $i (keys %savedcons) {
      printf "   %15s",$i;
      print "\n" unless ($j=1-$j);
    }
    print "\n" if $j;
  }
  unless($j = $savedcons{$_}) {
    print "No configuration stored as '$_'.\n";
    return;
  }
  eval($j);
  print "\n  $_ restored.\n\n";
}
  
sub sub_loadatm {
  local($nm) = &expand_file($_[0]);
  (warn "Error opening `$nm': $!\n\n"),return -1 unless open(CIN,$nm);
  $con{'name'} = $_[0];
  $plotstyle{'file'} = "$_[0].eps";
  $plotstyle{'title'} = $_[0] if $plotstyle{'title'};
  $_ = <CIN>;
  (warn "Error in '$_[0]': first line does not contain 'random'.\n"),
    return -1 unless /random/i;
  $con{'seed'} = &getnums(1);              $_ = <CIN>;
  $con{'time'} = &getnums(1);              $_ = <CIN>;
  @box = &abs(&getnums(3));                $_ = <CIN>;
  $con{'shear'} = join(',',&getnums(3));   $_ = <CIN>;
  @con{'pmass','vol'} = &getnums(2);       $_ = <CIN>;
  @con{'nperm','nharm'} = &getnums(2);     $_ = <CIN>;
  ($nt) = &getnums(1);                     $_ = <CIN>;
  @typen = &getnums(10);                   $_ = <CIN>;
  @foo = &getnums(10);

  die "Fatal error: < 1 component.. ($nt).\n" if $nt < 1;
  die "Fatal error: too many components ($nt).\n" if $nt > 10;
  $#typen = $#foo = $nt-1;

  $con{'masses'} = join(',',@foo);

  $natm = 0;  $#atm = -1;
  foreach $i (@typen) {$natm += $i;}
  $#atm = $natm*$nd-1;

  die "Fatal error: $natm atoms!\n" if $natm < 0 || $natm > 100000;

  local($i,$j,$n,$t);
  for($t=$n=0; $t < $nt; $t++) {
    chop($comments[$t] = <CIN>);
    $_ = <CIN>;
#    print "Reading $typen[$t] atoms for component $t ~ header was $_";
    for($i=0; $i < $typen[$t]; $i++,$n++) {
      $j = $n*$nd;
      $_ = <CIN>;
      ($atm[$j],$atm[$j+1],$atm[$j+2],$atm[$j+3])
          = /^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
      ( $atm[$j+$aref{'type'}],$atm[$j+$aref{'radius'}],
	$atm[$j+$aref{'color'}],$atm[$j+$aref{'sphere'}] ) = ($t,$defrad,$defcol,$defsphere);
#      print "Atom $n ($i) is at @atm[$j .. $j+4].\n";
    }
  }
  $_ = <CIN>;
  if(defined($_)) {
    $#velocs = $natm*3-1;
    for($t=$n=0; $t < $nt; $t++) {
      $_ = <CIN> if $t;
      for($i=0; $i < $typen[$t]; $i++,$n++) {
        $_ = <CIN>;
        @velocs[$n*3,$n*3+1,$n*3+2] = /^\s*(\S+)\s+(\S+)\s+(\S+)/;
      }
    }
    ($#velocs=0),warn "Warning: Velocities terminated early.\n" if !defined $_;
  } else {
    $#velocs = 0;
  }
  close(CIN);
  @group = (0 .. $natm-1);
  foreach $i (keys %saved) {
    $saved{$i} = join(',',grep($_ < $natm,split(',',$saved{$i})));
  }
  return 0;
}
sub sub_loadnew {
    $nbonds="";			# added 10/4/96 BPU to erase old bonds
				# after a new load
    return unless ($_ = &get_input(" Enter file to load: "));
    return if &sub_loadatm($_) < 0;
    &sub_coninfo(0);
}
sub numsort {$a <=> $b;}
sub sub_saveatm {
  local(@group,$i,$j,$k,$n,$knm,@enums,@foo) = @_;
  $_ = &get_input("Enter filename to save to: [out.con] ");
  $_ = "out.con" if $_ =~ /^\s*$/;
  $_ = &expand_file($_);
  (warn "Error: Could not save to file `$_': $!\n"),return
    unless open(COUT,">$_");
  $knm = ($settings{'save'} eq 'renumbered')?1:0;
  printf(COUT " %-54d Random Number Seed\n %-54g Time\n",@con{'seed','time'});
  printf(COUT " %-17.10g %-17.10g %-17.10g  Box: ax,ay,az\n",@box);
  printf(COUT " %-17.10g %-17.10g %-17.10g  Shear: alpha,beta,gamma\n",
         split(',',$con{'shear'}));
  printf(COUT " %-17.10g %-35.10g  Piston Mass,Volume\n",@con{'pmass','vol'});
  printf(COUT " %-17.10g %-35.10g  #permanent,#harmonic\n",
         @con{'nperm','nharm'});
  printf(COUT " %-54d #Components\n",$nt);
  @group = sort numsort @group;
  $#enums = $nt-1;
  foreach $i (@group) {$enums[$atm[$i*$nd+$aref{'type'}]]++;}
  print COUT " @enums",(" " x (54-length("@enums")))," #Atoms/component\n";
  @foo = split(/,/,$con{'masses'});
  print COUT " @foo",(" " x (54-length("@foo")))," Masses\n";
  ($j,$n) = (-1,0);
  foreach $i (@group) {
    $k = $i*$nd;
    print COUT "$comments[++$j]\n  COORDINATES OF COMPONENT # ",($j+1),"\n"
      while($j != $atm[$k+$aref{'type'}]);
    printf COUT " % -21.17g  % -21.17g  % -21.17g  %3d  %3d\n",
      @atm[$k,$k+1,$k+2,$k+3],($knm?$i:$n)+1;
    $n++;
  }
  if($#velocs > 1) {
    ($j,$n) = (-1,0);
    foreach $i (@group) {
      print COUT "  VELOCITIES FOR COMPONENT # ",(++$j+1),"\n"
        while($j != $atm[$i*$nd+$aref{'type'}]);
      printf COUT " % -23.17g  % -23.17g  % -23.17g  %3d\n",
        @velocs[$i*3,$i*3+1,$i*3+2],($knm?$i:$n)+1;
      $n++;
    }
  }
  close(COUT);
}
sub sub_coninfo {
  local($flag,$i,$j,@bounds) = @_;
  print " System `$con{'name'}' has $nt component" . ($nt==1?"":"s") .
        " and $natm atom" . ($natm==1?"":"s");
  print " total (@typen)" if $nt>1;
  printf(".\n Box: %g, %g, %g\n\n",$box[0], $box[1], $box[2]);
  return unless $flag;
  print " Time = $con{'time'}\t Mass",($nt>1?"es:":" =")," $con{'masses'}\n\n";
}

sub shift_con {
  local(@group,@trans,$i,$j,$k,@xyz) = @_;
  if($#FN_ARGS >= 0 && $FN_ARGS[0] =~ /^drop (\d+),?\s*,?(\d+),?\s*,?(\d+)/) {
    if($1 < 1 || $2 < 1 || $3 < 1 || $1>= $natm || $2>= $natm || $3>= $natm) {
      print " Invalid atom indices.\n";
      return;
    }
    @trans = &_drop_site($1-1,$2-1,$3-1);
    if($trans[0]<-1e6) {
      print " Error dropping to $1,$2,$3.\n";
      return;
    }
    for $i (0..2) { $trans[$i] -= $atm[$_[0]*$nd+$i] }
  } else {
    $_ = &get_input(" Enter dx,dy,dz for shift (or '=atom' for new center): ");
    if(/=\s*atom/) {
      print " You're supposed to enter an atom index, silly!\n\n";
      return;
    }
    s/=//,$j=1 if /=/;
    @trans = &getnums(3);
    if($j && ( ($j=$trans[0]*$nd)<1 || $j>$natm*$nd )) {
      print "Error: Invalid atom index entered ($j).\n";
      return;
    }
    ($j-=$nd),@trans = (-$atm[$j],-$atm[$j+1],-$atm[$j+2]) if $j;
  }
  if($#trans != 2) {
    $#trans = 2;
    print "Using {",join(',',@trans),"} for translation.\n";
  }
  $_ = &get_input(" Use periodic boundary conditions? [y] ");
  $p = (/^\s*n/i)?0:1;
  for $i (@group) { for $j (0..2) {
    $atm[$k = $i*$nd+$j] += $trans[$j];
    if($p) {
      $atm[$k] -= $box[$j] if $atm[$k]>$box[$j]/2;
      $atm[$k] += $box[$j] if $atm[$k]<-$box[$j]/2;
    }
  } }
}
sub scale_con {
  local(@cent,@scal,@group,$f) = @_;  @group = splice(@cent,3);
  $_ = &get_input(" Enter x,y,z for scaling: ");
  @scal = &getnums(3);
  $#scal = 2;
  for $i (0..2) {$f = $scal[$i] = 1 unless $scal[$i]}
  print " Using {",join(',',@scal),"} for scaling..\n" if $f;
  foreach $i (@group) { for $j (0..2) {
    $atm[$i*$nd+$j] = $cent[$j] + ($atm[$i*$nd+$j] - $cent[$j])*$scal[$j];
  } }
  for $i (0..2) {
    $box[$i] = 2*&abs($cent[$i] +
                     (($cent[$i]<0?1:-1)*$box[$i]/2 - $cent[$i])*$scal[$i]);
  }
}
sub printmat {
  local($nr,$nc,@mat,$n) = @_;
  for $i (0 .. $nr-1) {
    for $j (0 .. $nc-1) { printf(" %8.4g",$mat[$n++]); }
    print "\n";
  }
}
#Added APS June 2, 1994
sub sub_bonds {
  local(@group,$len, $i, $j, $k, $n, $x, $s, @d) = @_;
  $_ = &get_input(" Enter maximum bond length: ");
  $len = &getnums(1);
  (warn "Error: bond length???\n"), return unless $len > 0.000001;
  $n = $nbonds;
  $bond_dash=0;
  $_ = &get_input(" Make bonds dashed lines? [n] ");
  $bond_dash = 1 if /^\s*[yY]/;
  $_ = &get_input(" Enter length at which bonds become dashed: ");
  $dash_len=&getnums(1);
  $_=&get_input(" Enter width of this bond: (0 is normal) ");
  $bond_width=&getnums(1);
  $bond_shade=0;
  $_=&get_input(" Shade this bond? [n] ");
  $bond_shade=1 if /^\s*[yY]/;
  if ($bond_shade) {
      $_=&get_input(" Red component of bond: ");
      $bond_red=&getnums(1);
      $_=&get_input(" Green component of bond: ");
      $bond_green=&getnums(1);
      $_=&get_input(" Blue component of bond: ");
      $bond_blue=&getnums(1);
  } 
  foreach $i (@group) {
    foreach $j (@group) {
      next unless ($i < $j);
      $x = $n*$nbd; # Position in the bonds array...
      @d = ($atm[$i*$nd] - $atm[$j*$nd], 
            $atm[$i*$nd+1] - $atm[$j*$nd+1],
            $atm[$i*$nd + 2] - $atm[$j*$nd + 2]);
#      printf("Atoms %d, %d, disp %g %g %g\n", $i, $j, $d[0], $d[1], $d[2]);
      $s = $d[0]*$d[0] + $d[1]*$d[1] + $d[2]*$d[2];
      if ($s < $len*$len) {
        $bonds[$x] = $i;
        $bonds[$x + 1] = $j;
        $bonds[$x + $bref{'length'}] = sqrt($s);
	$bonds[$x + $bref{'width'}] = $bond_width;
	$bonds[$x + $bref{'shade'}] = $bond_shade;
        $bonds[$x + $bref{'red'}] = $bond_red;
        $bonds[$x + $bref{'green'}] = $bond_green;
        $bonds[$x + $bref{'blue'}] = $bond_blue;
        printf("Bond %d, atom %d to atom %d, length %g\n", $n,
		 $bonds[$x]+1, $bonds[$x + 1]+1, $bonds[$x + $bref{'length'}]);
        $n++;
      }
    }
  }
  $nbonds = $n;
}
# End of addition
sub mirror_con {
    local(@group,$plane,$pos,@foo,$x1,$y1,$z1,$i,$x)=@_;
    $x1=1;$y1=1;$z1=1;
    $_ = &get_input(" Enter axis normal to mirror plane [x, y or z] ");
    if (/[xyz]/i) {
	$x1=-1 if /x/i;
	$y1=-1 if /y/i;
	$z1=-1 if /z/i;
    }
    
    $_ = &get_input(" Enter position of mirror plane: ");
    $pos=&getnums(1);
    @foo=((1-$x1)*$pos,(1-$y1)*$pos,(1-$z1)*$pos);
    print $foo[0]." ".$foo[1]." ".$foo[2]." $x1 $y1 $z1\n";
    foreach $i (@group) {
	$x=$i*$nd;
#	print "B: ".$atm[$x]." ".$atm[$x+1]." ".$atm[$x+2];
	@atm[$x,$x+1,$x+2]=($foo[0]+$x1*$atm[$x],$foo[1]+$y1*$atm[$x+1],$foo[2]+$z1*$atm[$x+2]);
#	print "  A: ".$atm[$x]." ".$atm[$x+1]." ".$atm[$x+2]."\n";
    }
}
sub rot_con {
  local(@cent,@group,@foo,@bar,@mat,$i,$x,$y,$z,$ca,$sa) = @_;
  @group = splice(@cent,3);
  $_ = &get_input(" Enter rotation normal (or x/y/z): ");
  if(/[xyz]/i) {
    $x = 1 if /x/i;
    $y = 1 if /y/i;
    $z = 1 if /z/i;
  } else {
    ($x,$y,$z) = &getnums(3);
  }
  $i = sqrt($x*$x + $y*$y + $z*$z);
  (warn "Error: Vector has no magnitude!\n"),return unless &abs($i)>0.000001;
  ($x,$y,$z) = ($x/$i,$y/$i,$z/$i);

  $_ = &get_input(" Enter angle to rotate (counterclockwise): ");
  $i = &getnums(1) * $M_PI/180;
  ($ca,$sa) = (cos($i),sin($i));

  @mat =
    ($x*$x + $ca*(1-$x*$x),   $x*$y*(1-$ca) - $z*$sa,  $x*$z*(1-$ca) + $y*$sa,
     $x*$y*(1-$ca) + $z*$sa,  $y*$y + $ca*(1-$y*$y),   $y*$z*(1-$ca) - $x*$sa,
     $x*$z*(1-$ca) - $y*$sa,  $y*$z*(1-$ca) + $x*$sa,  $z*$z + $ca*(1-$z*$z) );
#  print "\nMATRIX:\n";&printmat(3,3,@mat);print"\n";

  foreach $i (@group) {
    $x = $i*$nd;
    @foo = ($atm[$x]-$cent[0], $atm[$x+1]-$cent[1], $atm[$x+2]-$cent[2]);
    for $j (0..2) {
      $bar[$j] = $mat[$j*3]*$foo[0]+$mat[$j*3+1]*$foo[1]+$mat[$j*3+2]*$foo[2];
    }
    @atm[$x,$x+1,$x+2] = ($bar[0]+$cent[0],$bar[1]+$cent[1],$bar[2]+$cent[2]);
  }
}
sub set_rgbcolors {
  local($n,$foo) = $aref{'color'};
  $foo = 1,shift @FN_ARGS if $FN_ARGS[0] =~ /depth/i;
  foreach $i (0..2) {
    unshift(@FN_ARGS,'by depth') if $foo;
    grep($atm[$_*$nd+$n]=$rgbs[$_*3+$i],@_);
    print (('Red','Green','Blue')[$i]),":\n";
    &set_colors(@_);
    grep($rgbs[$_*3+$i] = $atm[$_*$nd+$n], @_);
    grep($atm[$_*$nd+$n]=-1,@_);
  }
}
sub set_colors {
  local($h,$c,$d,$n,$p) = $aref{'color'};
  if($#FN_ARGS >= 0 && $FN_ARGS[0] =~ /depth/i) {
    shift @FN_ARGS;
    print " By depth:\n";
    $_ = &get_input(" Enter coordinate axis (x,y,z): [z] ");
    $_ = (/([xyz])/)?$1:'z';    $_ =~ tr/xyz/012/;
    $n = $_;
    $_ = &get_input(" Enter min, max intensities: ");
    $p=1 if /[\+-]/;
    ($c,$d) = /^\s*\+?([-\d.]+)[\s,]+\+?([-\d.]+)/;
    local($f,$g,$s) = (1e32,-1e32);
    foreach $i (@_) {$f = &min($f,$atm[$i*$nd+$n]);
                       $g = &max($g,$atm[$i*$nd+$n]);}
    $s = 'foreach $i (@_) { $atm[$i*$nd+$h]' . ($p?'+=':'=') .
                    '(($atm[$i*$nd+$n]-$f)/($g-$f))*($d-$c) + $c;}';
    eval($s);
  } else {
    $_ = &get_input(" Enter intensity: ");
    $p = 1 if /[\+-]/;
    (print "Error: '$1' is not a valid color.\n"),return if
                                                !/^\s*\+?([-\d.]+)/;
    grep($atm[$_*$nd+$h] = $1, @_) unless $p;
    grep($atm[$_*$nd+$h] += $1, @_) if $p;
  }
  if($p) {foreach $i (@_) {@atm[$i*$nd+$h]=&bound_vals(0,1,$atm[$i*$nd+$h])}}
}

sub set_toggle {
    print "Sphere toggled for this selection\n";
    foreach $i (@_) {@atm[$i*$nd+$aref{'sphere'}]=1-@atm[$i*$nd+$aref{'sphere'}];}
}

sub set_tags {
  local($n,$c) = $aref{'tag'};
  $c = &get_input(" Enter number to tag atoms: ");
  grep($atm[$_*$nd+$n] = $c, @_);
}
sub set_types {
  local($n,$c) = $aref{'type'};
  $_ = &get_input(" Enter component for atoms: ");
  $c = &getnums(1) || ((print "Error: $_ is not a valid component.\n"),return);
  grep($atm[$_*$nd+$n] = $c-1, @_);
}
sub set_radii {
  local($c,$d,$e,$f,$r,$b);
  print " Enter radius for atoms OR touch for max\n";
  $_ = &get_input("                        OR depth for by depth: ");
  if (/^\s*depth/i) {		# added 10/4/96 BPU to let user enter
				# depth as an option to get radius
				# by depth
    print " By depth:\n";
    $_ = &get_input(" Enter coordinate axis (x,y,z): [z] ");
    $_ = (/([xyz])/)?$1:'z';    $_ =~ tr/xyz/012/;
    $n = $_;
    $_ = &get_input(" Enter min, max radii: ");
    ($c,$d) = /^\s*\+?([-\d.]+)[\s,]+\+?([-\d.]+)/;
    local($mn,$mx,$s) = (1e32,-1e32);
    foreach $i (@_) {$mn = &min($mn,$atm[$i*$nd+$n]);
                       $mx = &max($mx,$atm[$i*$nd+$n]);}
    $b=$aref{'radius'};
    $s = 'foreach $i (@_) { $atm[$i*$nd+$b] = 
                    (($atm[$i*$nd+$n]-$mn)/($mx-$mn))*($d-$c) + $c;}';
    eval($s);
  } else { 
      if (/^\s*touch/i) {
	  print " Calculating max radius.." if $#_>80;
	  $r = 1e6;
	  for $i (0 .. $#_) { for $j (($i+1) .. $#_) {
	      ($c,$d) = ($nd*$_[$i],$nd*$_[$j]);
	      $e = $atm[$c]-$atm[$d]; # 
	      $f = $atm[$c+1]-$atm[$d+1];
	      $e = $e*$e + $f*$f;
	      $f = $atm[$c+2]-$atm[$d+2];
	      $e += $f*$f;
	      $r = $e if $r>$e;
	  } }
	  $r = sqrt($r)/2;
	  print ".done.  " if $#_>80;
	  printf " Radius was %.4g.\n",$r;
	  $c = $r;
      } else {
	  if(($c = &getnums(1)) < 0.000001) {
	      print "Error: '$_' is not a valid radius.\n";
	      return;
	  }
      }
      grep($atm[$_*$nd+$aref{'radius'}] = $c, @_);
  }
}	

sub set_labels {
  local($k,$i,$j,$_) = ($aref{'label'});
  $_ = &get_input(" Enter label ('' to rm, '=tag' for tag): [atom number] ");
  if(/^\s*''/) {
    grep($atm[$_*$nd+$aref{'label'}] = '', @_);
    return;
  }
  ($j = $1),s/^\^[\d.]+\^// if /^(\^[.\d]+\^)/;
  if(/^\s*=\s*tag/) {
    grep($atm[$_*$nd+$k] = $atm[$_*$nd+$aref{'tag'}], @_);
  } else {
    $i = &psstring($_) unless /^\s*$/;
    grep($atm[$_*$nd+$k] = $i?$i:$_+1,@_);
  }
  grep($atm[$_*$nd+$k] = "$j$atm[$_*$nd+$k]",@_) if $j;
}

sub sub_selbytype {
  local(@sel,$i);
  (print " (only one component present)\n"),return if $nt<2;
  $i = '(s)' if $nt>2;
  $_ = &get_input(" Enter component$i to select: ");
  @sel = &abs(split(/[^\d]/,$_));
  (print " (no components given)\n"),return if $#sel<0;
  @sel = &arr_to_tbl($nt,@sel);
  grep($sel[$atm[$_*$nd+$aref{'type'}]+1], 0 .. $natm-1);
}
sub sub_selbytag {
  local(@sel);
  $_ = &get_input(" Enter tag(s) to select: ");
  @sel = split(/[^\d]/,$_);
  (print " (no tags given)\n"),return if $#sel<0;
  @sel = &arr_to_tbl($nt,@sel);
  grep($sel[$atm[$_*$nd+$aref{'tag'}]], 0 .. $natm-1);
}

sub set_cellsize {
    $_ = &get_input(" Enter lengths x,y,z for sides of cell: ");
    @box=&abs(&getnums(3));
}

sub subbox_con {
  local($centered,@cent,@box,@bnds,@grp,$i,$j) = @_;
  if($centered) {
    $_ = &get_input(" Enter lengths x,y,z for sides of box: ");
    @box=&abs(&getnums(3));
    $i=$box[0],@box=($i,$i,$i) if !$#box;
    $#box=2,printf("Using {%g, %g, %g} for box.\n",$box[0],$box[1],$box[2])
      unless $#box==2;
    for $i (0..2) {
      $bnds[$i+3] = ($bnds[$i] = $cent[$i]-$box[$i]/2) + $box[$i];
    }
  } else {
    $_ =&get_input(" Enter sides of box: (xmin,xmax, ymin,ymax, zmin,zmax): ");
    @box = &getnums(6);
    if($#box!=5) {$#box = 5; print "Using {@box} for box.\n";}
    @bnds = ($box[0],$box[2],$box[4],$box[1],$box[3],$box[5]);
  }
  grep($atm[$_*$nd] >= $bnds[0]   && $atm[$_*$nd] <= $bnds[3] &&
       $atm[$_*$nd+1] >= $bnds[1] && $atm[$_*$nd+1] <= $bnds[4] &&
       $atm[$_*$nd+2] >= $bnds[2] && $atm[$_*$nd+2] <= $bnds[5], 0 .. $natm-1);
}
sub subsph_con {
  local(@cent,$r) = @_;
  while(1) {
    $_ = &get_input(" Enter radius for sphere: ");
    ($r) = &getnums(1);
    last if $r > 0;
    print "\n Radius $r is invalid!\n";
  }
  $r *= $r;
  grep(($atm[$_*$nd]-$cent[0])**2 + ($atm[$_*$nd+1]-$cent[1])**2 +
       ($atm[$_*$nd+2]-$cent[2])**2 < $r, 0 .. $natm - 1);
}
sub sub_prtatm {
  local($npl,@group,$i,$j) = @_;
  $npl = $FN_ARGS[0] unless $#FN_ARGS;
  (print "(No atoms in set)\n\n"),return if $#group<0;
  for($i=$j=0;$i<=$#group;$i++) {
    printf(" %3d",$group[$i]+1);
    (print "\n"),$j=0 if ++$j==$npl;
  }
  print "\n" if $j;
  print "\n";
}
sub sub_atminfolist {
  local(@grp) = @_;
  local($ty,$ta,$ra,$co,$i,$j,$n,$v) = @aref{'type','tag','radius','color'};
  (print "(No atoms given)\n\n"),return if $#_<0;
  shift(@FN_ARGS),$v = 1
              if ($#velocs > 1 && $#FN_ARGS >= 0 && $FN_ARGS[0] =~ /^-[vV]/);
  shift(@FN_ARGS),@grp=splice(@grp,0,$1) if
                             $#FN_ARGS >= 0 && $FN_ARGS[0] =~ /^\s*(\d+)/;
  $i = '  Type' if $nt>1;
  $j = '           Velocity        ' if $v;
  print"   #$i  Tag   Radius     Color           Position       ${j}  Label\n";
  $i = '-' x length($i) if $nt > 1;
  $j = '-' x (length($j)-1) if $v;
  print"----$i--------------------------------------------------${j}-------\n";
  foreach $i (@grp) {
    $n = $i*$nd;
    printf('%4d',$i+1);
    printf(' %5d',$atm[$n+$ty]+1) if $nt>1;
    printf(' %4d  %7.5g ',$atm[$n+$ta],$atm[$n+$ra]);
    $j = $atm[$n+$co];
    printf('%9.5g  ',$j) if $j>=0;
    printf('%3.1g/%3.1g/%3.1g',@rgbs[$i*3,$i*3+1,$i*3+2]) if $j<0;
    printf(' % 7.2g,% 7.2g,% 7.2g',@atm[$n,$n+1,$n+2]);
    printf('   % 7.2g,% 7.2g,% 7.2g',@velocs[$i*3,$i*3+1,$i*3+2]) if $v;
    print "  $atm[$n+$aref{'label'}]\n";
  }
  print "\n";
}
sub _getcenterpt {
  local($_,$i,@cent);
  $_ = &get_input(" Enter atom or x,y,z position to use as center: [0,0,0] ");
  while(1) {
    return (0,0,0) if /^\s*0*\s*$/;
    if(/^\s*(\d+)\s*$/) {
      if(($i = $1) < 1 || $i>$natm) {
        print "Invalid atom ($cent[0]).  $natm atoms are present.\n\n";
        next;
      }
      $i--;
      @cent = @atm[$i*$nd,$i*$nd+1,$i*$nd+2];
      printf " ..using atom ",$i-1," @ (%.4g, %.4g, %.4g).\n",@cent;
      return @cent;
    }
    @cent = &getnums(3);
    $#cent = 2;
    return @cent;
  } continue {
    print " Enter atom or x,y,z position to use as center: [0,0,0] ";
    $_ = <STDIN>;
  }
}

sub add_atoms {
  local(@a,@p,@tbl,$_,$j,$n);
  $_ = &get_input(" Enter x,y,z [,tag,component] for new atom: ");
  @p = &getnums(5);
  $p[4] = &max($p[4],0);
  $nt = &max($nt,$p[4]+1);
  $#a = $nd-1;
  @a[0,1,2,@aref{'tag','type','radius','color','sphere'}]
    = (@p[0,1,2,3,4],$#_>=0?$_[0]:$defrad,$defcol,$defsphere);
  if(&sub_getopt('settings') eq 'on') {
    @group = grep(++$_,@group);
    foreach $i (keys %saved) {
      $saved{$i} = join(',',grep(++$_,split(/,/,$saved{$i})));
    }
    unshift(@atm,@a);
    unshift(@velocs,(0,0,0)) if $#velocs>=0;
    unshift(@rgbs,('1,1,1')) if $#rgbs>0;
    $n = 1;
  } else {
    push(@atm,@a);
    push(@velocs,(0,0,0)) if $#velocs>=0;
    push(@rgbs,('1,1,1')) if $#rgbs>0;
    $n = $natm+1;
  }
  $natm++;
  for $i (0..2) {($j=1),last if $p[$i] < -$box[0]/2 || $p[$i] > $box[0]/2;}
  if($j) {
    $_ = &get_input(" Atom is outside box.  Enlarge box? [y] ");
    if(!/^\s*n/i) {
      for $i (0..2) {$box[$i] = &max($box[$i],&abs(2*$p[$i]));}
    }
  }
  print " New atom is # $n.\n";
}
sub _drop_site {
  local(@i,$r,$a,$b,$c,$d,@j,@k,@l) = @_;
  $r = 1e6;
  for $i (0..$#i) { for $j ($i+1..$#i) {
    ($c,$d) = ($i[$i]*$nd,$i[$j]*$nd);
    $a = $atm[$c]-$atm[$d];
    $b = $atm[$c+1]-$atm[$d+1];
    $a = $a*$a + $b*$b;
    $b = $atm[$c+2]-$atm[$d+2];
    $a += $d*$d;
    $r = $d if $r>$d;
  } }
  if($r < 0.0000001 || $r > 1e3) {
    printf "\n Oops, min distance between %d,%d,%d is %.4g!",@i[0..2],$r;
    return (-1e6,-1e6,-1e6);
  }
  for $i (0..2) {
    $l[$i] = ($atm[$i[0]*$nd+$i]+ $atm[$i[1]*$nd+$i]+ $atm[$i[2]*$nd+$i])/3;
  }
  @j = @atm[$i[0]*$nd, $i[0]*$nd+1, $i[0]*$nd+2];
  @k = @atm[$i[1]*$nd, $i[1]*$nd+1, $i[1]*$nd+2];
  for $i (0..2) {
    $j[$i] -= $l[$i];
    $k[$i] -= $l[$i];
  }
  @j=($j[1]*$k[2]-$j[2]*$k[1],$j[2]*$k[0]-$j[0]*$k[2],$j[0]*$k[1]-$j[1]*$k[0]);
  @j = (-$j[0],-$j[1],-$j[2]) if $j[2]<0;
  $a = $k[0]*$k[0] + $k[1]*$k[1] + $k[2]*$k[2];
  $a = sqrt($r/($r-$a));
  @l=($l[0]+$a*$j[0], $l[1]+$a*$j[1], $l[2]+$a*$j[2]);
  return @l;
}
sub drop_atom {
  local($_,$n,$r,@i,@j,@k,@tbl);
  $_ = &get_input(" Enter at least three atoms for drop site (or 'help'): ");
  if(/help/i) {
    print <<"EOH";
 If you entered n atoms, each unique subgroup of 3 atoms is used to calculate
 a position equidistant from the atoms (the one with the largest z value is
 preferred).  An atom is then placed at the average of all these positions.

EOH
    return;
  }
  @i = grep($_>0 && $_<=$natm && !$tbl[--$_]++,&getnums(100));
  (print "\n Drop aborted.\n"),return if $#i<2;
  for $i (0..$#i) { for $j ($i+1..$#i) { for $k ($j+1..$#i) {
    @j = &_drop_site(@i[$i,$j,$k]);
    next if $j[0]<-1e4 && $j[1]<-1e4 && $j[2]<-1e4;
    $k[0] += $j[0];  $k[1] += $j[1];  $k[2] += $j[2];
    $n++;
  } } }
  (print " All triads failed.  Placement aborted.\n"),return if $n<1;
  for $i (@i) { $r += $atm[$i*$nd + $aref{'radius'}] }
  @k = ($k[0]/$n, $k[1]/$n, $k[2]/$n);
  @FN_ARGS = ("$k[0],$k[1],$k[2]");
  &add_atoms($r/($#i+1));
}
sub atom_distances {
  local(@a,@tbl,$a,$b,$c,$d,$_);
  $_ = &get_input(" Enter atoms to find distance between: ");
  @a = sort numsort grep($_>0 && $_<=$natm && !$tbl[--$_]++, &getnums(100));
  (print " Aborted.\n"),return if $#a<1;
  print "\n";
  for $i (0..$#a) { for $j ($i+1..$#a) {
    ($c,$d) = ($a[$i]*$nd,$a[$j]*$nd);
    $a = $atm[$c]-$atm[$d];
    $b = $atm[$c+1]-$atm[$d+1];
    $a = $a*$a+$b*$b;
    $b = $atm[$c+2]-$atm[$d+2];
    printf "Distance %d <--> %d == %.4g\n",$a[$i]+1,$a[$j]+1,sqrt($a+$b*$b);
  } }
  print "\n";
}

sub get_saved_group {
  local($s) = @_;
  return (0, 0 .. $natm-1) if $s eq "all";
  return (0,@group) if $s eq "sel";
  (warn "Error: No such saved set `$1'.\n"),return (-1) if !$saved{$s};
  return (0,split(',',$saved{$s}));
}
sub load_saved_group {
  local(@foo);
  while(1) {
      $_ = &get_input(" Name of set to retrieve (<cr> to list): ");
      (print "Sets defined: ",join(', ',keys(%saved)),"\n\n"),next
	  unless /^\s*(\S+)/;
      @foo = &get_saved_group($1);
      shift @foo;
      return @foo;
  }
}
sub save_group {
  while(1) {
    $_ = &get_input(" Enter name for set: ");
    (print "Name '$_' is invalid.\n"),next unless /^\s*(\S+)/;
    last;
  }
  $saved{$_} = join(',',@_);
}
sub sub_getopt {
  local($v,@val) = @_;
  return "" unless (@val = split(',',$settings{$v}));
  $val[0];
}
sub _ssconv {
  local(@foo) = @_;
  foreach $i (@foo) {
    next if $i =~ s/\+/any positive number/;
    next if $i =~ s/\*/any number/;
    $i =~ s/0/any nonnegative number/;
  }
  @foo;
}
sub alphabetically { $a cmp $b;}
sub sub_setopts {
  local(@opts,$ok,$p);
  if($#FN_ARGS < 0) {
    foreach $opt (sort alphabetically keys(%settings)) {
      @opts = split(',',$settings{$opt});
      printf " %10s = %-20s ",$opt,shift(@opts);
      print "(domain: ",join(', ',&_ssconv(@opts)),")\n\n";
    }
    return;
  }
  (print "Error: no such option '$FN_ARGS[0]'.\n\n"),return
    if !$settings{$FN_ARGS[0]};
  @opts = split(',',$settings{$FN_ARGS[0]});
  (printf " %10s = %-20s\n\n",$FN_ARGS[0],$opts[0]),return if !$#FN_ARGS;
  shift(@opts);
  $ok = 0;
  $p = $opts[0];
  $ok = 1,$p *=1 if $p =~ /^\*/;
  $ok = 1,$p *=1 if $p =~ /^0/ && $FN_ARGS[1]>=0;
  $ok = 1,$p *=1 if $p =~ /^\+/ && $FN_ARGS[1]>0;
  $ok = 1 if $p !~ /^[\*0\+]/ && grep($_ eq $FN_ARGS[1],@opts);
  (printf "Error: $FN_ARGS[1] not a valid value for $FN_ARGS[0].\n\n"),
    return unless $ok;
  $settings{$FN_ARGS[0]} = join(',',$FN_ARGS[1],@opts);

  $defsetmode = &sub_getopt("selmode");
}
sub do_alias {
  local($name,$val) = @_;
  if(!$name) {
    foreach $name (sort alphabetically keys(%aliases)) {
      printf "%10s = %s\n",$name,$aliases{$name};
    }
    return;
  }
  if(!$val) {
    print "Error: no alias for '$name' present.\n" if !defined $aliases{$name};
    print "%10s = %s\n",$name,$aliases{$name} if defined $aliases{$name};
    return;
  }
  $aliases{$name} = $val;
}
sub do_unalias {
  local($name) = @_;
  (print "Error: no alias defined for '$name'\n"),return
                                          unless defined $aliases{$name};
  delete $aliases{$name};
}

################################
#---------------------
# PRIMARY CODE BEGINS

$| = 1;

%plotstyle = ();                       # plotting setup
@pagesize = (8.5 * 72, 11.0 * 72);     # constant pagesize.. ps: 1pt = 1/72 in.
@margins = (1.0 * 72, 1.0 * 72);       # constant margins..

%settings = ('prompt','plain,plain,dir', 'selmode','set,add,int,set,sub,xor',
             'pager','off,on,off,less,more', 'menu','on,on,off',
             'save','renumbered,renumbered,prenumbered',
             'diamonds','1,+', 'Xsize','1,+', 'postcom',',*',
             'landscape','vertical,horizontal,vertical',
             'bbox','minimal,large,margins,minimal',
             'boxdraw','silent,query,silent',
             'plotfontsz','24.0,+', 'plotfont','Helvetica,*',
             'linewidth','0,*',
             'dropasone','on,on,off');

%aliases = ('ls','ls -CF');

%saved = ();

$undoable = 0;      # last command is undoable (so show menu option)
$showmenu = 0;      # when zero, don't show menu unless it's been $nmenu
$didmenu = 0;       # was initial menu-showing done?
$nmenu = 2;         #   ..times since the menu was last showed.

$defsetmode = $setmode = "set";  # one of (add,int,set,sub,xor)
$defnotmode = 0;          # could be 1
$notmode=0;               # take complement of group before doing modeop

$doingop=0;               # used by interrupt_handler - die on INT if false.

$pager = $ENV{'PAGER'} || $ENV{'MANPAGER'} || 'more';

$home = $ENV{'HOME'} || $ENV{'LOGDIR'} || ".";
$user = $ENV{'USER'} || $ENV{'LOGNAME'} || "";

$rc = $rcfile if !$rc && -o $rcfile;
$rc = "$home/$rcfile" if !$rc && -o "$home/$rcfile";

&load_defaults($datfile);

&show_version(),exit 0 if $#ARGV >= 0 && $ARGV[0] =~ /^-?-v(ers(ion)?)?/;
$rcfile="",shift(@ARGV) if $#ARGV >= 0 && $ARGV[0] eq '-F';
shift(@ARGV),$rc=shift(@ARGV) if $#ARGV >= 0 && $ARGV[0] eq '-rc';
die "Usage: $progname [-F | -rc <file>] <file>\n" if $#ARGV;
exit 1 if &sub_loadatm($ARGV[0]) < 0;
&sub_coninfo(0);

if($rc) {
  open(RCFILE,$rc) || die "Couldn't open $home/$rcfile: $!\n";
  @commands = <RCFILE>;
  close(RCFILE);
}

$SIG{'INT'} = 'interrupt_handler';

do {
  select(STDOUT);
  print $@ if $@ && $@ !~ /aborted/i;
  eval("&mainloop");
} while $@ !~ /time to die/;

exit 0;

# PRIMARY CODE ENDS
#---------------------
################################

sub do_quit { &save_defaults($datfile); die 'time to die'; }

sub interrupt_handler {
  die 'time to die' unless $doingop;
  print STDERR "\nAborted.\n"; die "Aborted.";
}

sub show_menu {
  return if &sub_getopt("menu") eq 'off';
  &force_menu();
}
sub force_menu {
  $undostr = $undoable?"(u) undo":(" " x 8);
  $grppo = $#group+1;

  print <<"EOM";
 sel [add|int|set|sub|xor] [not] _command_:                   [$grppo selected]
  (sph)        spherical region         (save)          save current selection
  (box)        box region               (load | <name>) load saved selection
  (boxp)       box region about point   (list) or ()    list current selection
  (tag/comp)   by tag/component         (info)    list attributes of selection
  (1,2,3,..)   given atoms              (all | *)       all atoms
  (x..y | x-y) range of atoms

 _command_ sel|all|indices..|(name of saved selection):
  (rot)        rotate atoms                (save)        save atoms to file
  (scale)      scale atoms' positions      (plot)        generate plot of atoms
  (shift)      translate atoms             (label)       set atoms' labels
  (color)      color atoms                 (set radius)  set atoms' radii
  (tag)        set atoms' tags             (set comp)    set atoms' component
  (bond)       generate bonds              (toggle)      toggle spherical atoms
  (mirror)     mirror atoms thru plane

 _command_:
  (add)        add new atom(s)             (dist)        distance between atoms
  (drop)       drop an atom

  (load)       load new configuration      $undostr
  (info)       info on configuration       (help)        further features
  (!com)       run com as shell command    (quit | exit) quit

EOM
}

sub mainloop { $doingop = -1; while(1) {
  if($doingop >= 0 && ($foo = &sub_getopt('postcom'))) {
    @commands = ($foo);
    $doingop = -1;
  } else {
    $doingop = 0;
  }
  close(PAGER),select(STDOUT),$paging="" if $paging;
  if(@commands) {
    $o = shift @commands;
  } else {
    $didmenu=1,&show_menu() unless $didmenu;
    if(&sub_getopt("prompt") eq "dir") {
      chop($foo=`pwd`);
      print "[$foo] > ";
    } else {
      print " Command: ";
    }
    $o = <STDIN>;
    $o = "quit\n" unless defined($o);
  }
  $doingop = 1 unless $doingop;
  chop($o);
  if($o =~ /^\s*$/) {
    next if ++$nnsc <= $nmenu;
    &show_menu();
    $nnsc = 0;
    next;
  }
  if(!$pagerbad && ($pageopt = &sub_getopt("pager")) ne "off") {
    if(open(PAGER, ($pageopt eq "on")?"| $pager":"| $pageopt")) {
      select(PAGER); $paging=1;
    } else {
      warn "Warning: Pager '$pager' not found\n";
      $pagerbad = 1;
      $pager = $paging = 0;
    }
  }
  system($1),++$nnsc,next if $o =~ /^!(.*)/;
  if($o =~ /^\s*cd\s*(\S.*)?$/) {
    chdir($o = &expand_file($1)) || warn "Error: could not cd to $o: $!\n";
    $nnsc++; next;
  }
  if($o =~ /^\s*unalias(\s+(\S.*))?\s*$/) {
    (print "Usage error: unalias <name>\n"),next unless $2;
    &do_unalias($2);
    next;
  }
  if($o =~ /^\s*alias(\s+(\S.*))?\s*$/) {
    &do_alias('',''),next unless $2;
    $o = $2;
    $o =~ /^(\S+)\s+(\S.*)$/;
    &do_alias($1,$2);
    next;
  }

  @p = ();
  @o = split(/;/,$o);
  foreach $o (@o) {
    $o =~ /^\s*(\S+)/;
    $o =~ s/^\s*$1/$foo/ if ($foo=$aliases{$1});
    push(@p,split(';',$o));
  }

  $showmenu = 0;
  print "\n";
  foreach $mainloop_com (@p) {
    $o = $mainloop_com;
    @FN_ARGS = ();
    $showmenu++;

    #  Resolve aliases..
    $o =~ s/^\s*color(\W.*)?$/set color\1/;
    $o =~ s/^\s*tag(\W.*)?$/set tag\1/;
    $o =~ s/^\s*typify(\W.*)?$/set comp\1/;
    $o =~ s/^\s*label(\W.*)?$/set label\1/;
    $o =~ s/^\s*toggle(\W.*)?$/set toggle\1/;

    if($o =~ /^\s*sel(\s+(\S.*))?\s*$/) {
	($notmode,$setmode) = ($defnotmode,$defsetmode);
	$undoable = 1;   @undogroup = @group;   @undoatm = @atm;
	$o = $2;

	if($o eq "") {&sub_prtatm(10,@group);} else
	{
	    ($setmode,$o) = ($1,$2) if $o =~ /^\s*(add|int|set|sub|xor)\s+(\S.*)$/;
	    ($notmode,$o) = (1-$notmode,$1) if $o =~ /^\s*not\s+(\S.*)$/;

	    if($o =~ /^\d/) {
#  added range for sel with - or .. as indicator 12/19/96 BPU
		if ($o=~/(-+|\.\.)/) {
		    @group=();
		    @blas=split (/\D+/,$o);
		    for ($blas=$blas[0];$blas<=$blas[1];$blas++) {
			@group=(@group,$blas);
		    }
		} else {
		    @group=split(/\D+/,$o);
#		@group = split(/[^\d\+.]/,$o);
		}
		grep($_--,@group);
	    } else {
		($o,$foo,$args) = ($o =~ /^\s*([\*\w]+)(\s(\S.*))?$/);
		if($o =~ /^(sph|boxp?|tag|comp(onent)?|all|\*|load|save|list|info)$/) {
		    @FN_ARGS = split(/\//,$args);
  
		    @cent = &_getcenterpt() if $o =~ /^sph|boxp/;
		    $o eq "sph" && (@group = &subsph_con(@cent));
		    $o eq "boxp" && (@group = &subbox_con(1,@cent));
		    $o eq "box" && (@group = &subbox_con(0));
		    $o eq "tag" && (@group = &sub_selbytag());
		    $o =~ /^comp(onent)?$/ && (@group = &sub_selbytype());
		    ($o eq "all" || $o eq "*") && ((@group = (0 .. $natm-1)),$showmenu--);
		    $o eq "load" && (@group = &load_saved_group());
		    $o eq "save" && &save_group(@group);
		    $o eq "list" && &sub_prtatm(10,@group);
		    $o eq "info" && &sub_atminfolist(@group);
		} else {
		    last unless (@foo = &get_saved_group($o)) && $foo[0]>=0;
		    shift @foo;
		    @group = @foo;
		}
	    }
	}

      @group = &group_subtract($natm,0 .. $natm-1,@group) if $notmode;

      @group = &group_reduce(@group,@undogroup) if $setmode eq "add";
      @group = &group_intersect(@group,@undogroup) if $setmode eq "int";
      @group = &group_subtract($#undogroup+1,@undogroup,@group)
        if $setmode eq "sub";
      @foo=&group_xor(@group,@undogroup) if $setmode eq "xor";
    }
    elsif($o =~ /^\s*(plot|bond|mirror|rot|scale|shift|save|set\s+(\w+))/) {
      $foo = $1;
      if($o !~ /^\s*$foo\s+([\w,]+)(\s+(\S.*))?\s*$/) {
        print " Usage: $foo 'all'|'sel'|indices..|(name of saved set)\n\n";
        last;
      }
      ($bar,$bash) = ($1,$3);
      if($bar =~ /^\d/) {
        @grp = (1,split(/[^\d.]/,$bar));
        @grp = grep(--$_>=0 && $_<$natm,@grp);
      }
      elsif((@grp = &get_saved_group($bar)) && $grp[0]<0) {
        print " No such group '$bar' exists.\n";
        last;
      }
      shift @grp;
      ($o,@FN_ARGS) = ($foo,($foo eq 'load')?($bash):split(/\//,$bash));
      @cent = &_getcenterpt() if $o =~ /^(rot|scale)/;
      $o eq "save" && &sub_saveatm(@grp);
      $o eq "plot" && &sub_plotatm(@grp);
# Added APS June 2, 1994
      $o eq "bond" && &sub_bonds(@grp);
# End of addition (other "bond" additions above also)
      $o eq "mirror" && &mirror_con(@grp);
      $o eq "rot" && &rot_con(@cent,@grp);
      $o eq "scale" && &scale_con(@cent,@grp);
      $o eq "shift" && &shift_con(@grp);
      $o = $1 if $o =~ /^set\s+(\w+)/;
      $o eq "color" && &set_colors(@grp);
      $o eq "rgbcolor" && &set_rgbcolors(@grp);
      $o eq "tag"   && &set_tags(@grp);
      $o eq "comp" && &set_types(@grp);
      $o eq "radius" && &set_radii(@grp);
      $o eq "label" && &set_labels(@grp);
      $o eq "toggle" && &set_toggle(@grp);
    }
 elsif($o=~/^\s*(add|drop|dist|load|store|recall|info|help|quit|exit|menu|setopt|undo|Version|cellsize)/) {
      ($o,$foo,$args) = ($o =~ /^\s*(\w+)(\s+(\S.*))?$/);
      @FN_ARGS = ($o eq 'load')?($args):split(/[\s\/]/,$args);
      $o eq "cellsize" && &set_cellsize(@grp);
      $o eq "add" && &add_atoms();
      $o eq "drop" && &drop_atom();
      $o eq "dist" && &atom_distances();
      $o eq "load" && &sub_loadnew();
      $o eq "store" && &store_config();
      $o eq "recall" && &restore_config();
      $o eq "info" && &sub_coninfo(1);
      $o eq "help" && &sub_showhelp();
      $o eq "setopt" && &sub_setopts();
      $o eq "menu" && &force_menu();
      $o eq "Eval" && (eval($args),$showmenu--);
      $o eq "Version" && &show_version();
      ($o eq "quit" || $o eq "exit") && &do_quit();
      $o eq "undo" && $undoable && (@group=@undogroup,@atm=@undoatm);
      $undoable = 0 unless $o =~ /^imlp$/;
    }
    elsif($o=~/^Eval\s+(.*)$/) {
      eval($1);
      $showmenu--;
    } else {
      system($o);
      print "\n";
      $showmenu--;
    }
  }
  if($showmenu) {&show_menu();} else {$nnsc++;}
}}


#####################
#-------------------#
# Plotting Command  #
#-------------------#
sub viewpt_sort {
  local($d,$e)=($a*$nd,$b*$nd);
  $d = $atm[$d]*$viewpt[0] + $atm[$d+1]*$viewpt[1] + $atm[$d+2]*$viewpt[2];
  $e = $atm[$e]*$viewpt[0] + $atm[$e+1]*$viewpt[1] + $atm[$e+2]*$viewpt[2];
  $d <=> $e;
}
sub viewpt_diff{
  local($d,$e) = @_;
  $d *= $nd;
  $e *= $nd;
  $d = $atm[$d]*$viewpt[0] + $atm[$d+1]*$viewpt[1] + $atm[$d+2]*$viewpt[2];
  $e = $atm[$e]*$viewpt[0] + $atm[$e+1]*$viewpt[1] + $atm[$e+2]*$viewpt[2];
  return $d - $e;
}
sub sub_plotatm {
  local(@grp,@views,@foo,$foo,$i,$M,@shpnm) = @_;
  @shpnm = ('?','X\'s','diamonds','circles');

  if(!$plotstyle{'init'}) {    # note: $plotstyle{'file'} set by sub_loadatm..
    @plotstyle{'autobounds','style','color',
        'outline','radius','box','layout','views','title'}
        = (1, 2, -1 ,0, 1, 1, 'p', 'xy','');
  }

  $_ = &get_input(" Plot which view(s)? (xy, yz, xz, zx, ..) " .
                                         "[$plotstyle{'views'}] ");
  @views = ();
  $plotstyle{'views'} = $_ unless /^\s*$/;
  $_ = $plotstyle{'views'} unless /(xy|yz|yz|zy|xz|zx)/i;
  push(@views,1) if /xy/i;  push(@views,2) if /yx/i;
  push(@views,3) if /yz/i;  push(@views,4) if /zy/i;
  push(@views,5) if /xz/i;  push(@views,6) if /zx/i;

  if($M = $plotstyle{'init'}) {
    if($#FN_ARGS < 0) {
      push(@foo,"draw as $shpnm[$plotstyle{'style'}]");
      ($i,$j,$k) = @plotstyle{'color','outline','radius'};
      push(@foo,"color " . ($i<0?'by atom':$i));
      push(@foo,"outline " .($j<0?'by atom':"=$j")) if $plotstyle{'style'}!=1;
      push(@foo,"radius " . ($k<0?'by atom':$k)) if $plotstyle{'style'} == 3;
      push(@foo,$plotstyle{'box'}?"draw box":"no box");
      push(@foo,$plotstyle{'layout'} eq 'l'?"landscape":"portrait");
      print " (settings: " . join(',',@foo) . ")\n";
    }
    $_ = &get_input(" Use previous plot settings? [y] ");
    $M = 0 if /^\s*[nN]/;
  }
  if(!$M) {
    while(1) {
      $_ = &get_input(" Plot style? (1=X's,2=diamonds,3=circles) [" .
         "$plotstyle{'style'}] ");
      ($i) = (/^\s*\d/) ? &getnums(1) : ($plotstyle{'style'});
      last if $i > 0 && $i < 4;
      print "Error: Invalid input.\n";
    }
    $plotstyle{'style'} = $i;
    if($i == 2 || $i == 3) {
      $_ = &get_input(" Fill color (black=0 .. 1=white" .
        ($plotstyle{'color'}<0?'':',a=atom color') . ")? [" .
        ($plotstyle{'color'}<0?'a=atom color':$plotstyle{'color'}) . "] ");
      $_ =~ s/a/-1/;
      $plotstyle{'color'} = &getnums(1) unless /^\s*$/;
      $_ = &get_input(" Outline color? [" .
        ($plotstyle{'outline'}<0?'atom color':$plotstyle{'outline'}) ."] ");
      $_ =~ s/a/-1/;
      $plotstyle{'outline'} = &getnums(1)  unless /^\s*$/;
    } else {
      $_ = &get_input(" X color (black=0 .. 1=white" .
        ($plotstyle{'color'}<0?'':',a=atom color') . ")? [" .
        ($plotstyle{'color'}<0?'a=atom color':$plotstyle{'color'}) . "] ");
      $_ =~ s/a/-1/;
      $plotstyle{'color'} = &getnums(1) unless /^\s*$/;
    }
    if($i == 3) {
      $_ = &get_input(" Enter radius" .
        ($plotstyle{'radius'}<0?'':' (a=atom radius)') . ": [" .
        ($plotstyle{'radius'}<0?'atom radius':$plotstyle{'radius'}) . "] ");
      $_ =~ s/a/-1/;
      $plotstyle{'radius'} = &getnums(1) unless /^\s*$/;
    }
    $_ = &get_input(" Put box around boundaries? [" .
                             ($plotstyle{'box'}?'y':'n') . "] ");
    $plotstyle{'box'} = (/^\s*[yY]/)?1:0 unless /^\s*$/;
    $_ = &get_input(" Title: [$plotstyle{'title'}] ");
    $plotstyle{'title'} = $_ unless /^\s*$/;
    $plotstyle{'title'} = '' if /^\s*(none|'')\s*$/;
    $_ = &get_input(" Landscape (l) or Portrait (p) format? [" .
      ($plotstyle{'layout'} eq 'l'?'l':'p') . "] ");
    $plotstyle{'layout'} = (/^\s*[lL]/)?"l":"p" unless /^\s*$/;
  }
  $plotstyle{'init'} = 1;

  local($j,$k,$n,$l,$code,$land);
  local($x,$y,$sym,$col,$ocol,$rad,@trans,$scale,@pbounds,@pgsz,@marg,@bbox);
  ($sym,$col,$ocol,$rad) = @plotstyle{'style','color','outline','radius'};
  $col  = sprintf("%.4g",$col ) if $col>=0;
  $ocol = sprintf("%.4g",$ocol) if $ocol>=0;
  $rad  = sprintf("%.4g",$rad ) if $rad>=0;

  if(!$M) {
    $_ = &get_input(" Fit to page? (or set bounds manually) [" .
      ($plotstyle{'autobounds'}?'y':'n') . "] ");
    $plotstyle{'autobounds'} = (/^\s*y/i)?1:0 unless /^\s*$/;
  }
  $j = 0;
  @bounds = (1e32,1e32,1e32,-1e32,-1e32,-1e32);  # see use of box below..
  foreach $i (@grp) {
    @foo = @atm[$i*$nd,$i*$nd+1,$i*$nd+2];
    $j = ($rad < 0?$atm[$i*$nd+$aref{'radius'}]:$rad) if $sym == 3;
    @bounds[0,1,2] = (&min($bounds[0],$foo[0]-$j),
             &min($bounds[1],$foo[1]-$j),&min($bounds[2],$foo[2]-$j));
    @bounds[3,4,5] = (&max($bounds[3],$foo[0]+$j),
             &max($bounds[4],$foo[1]+$j),&max($bounds[5],$foo[2]+$j));
  }
  @pbounds = @bounds;
  $land = ($plotstyle{'layout'} ne 'l')?0:
                   ((&sub_getopt('landscape') eq 'horizontal')?2:1);

  foreach $i (@views) {
    while(1) {
      $j = " Enter output file (or '|command') for ";
      $j .= "file" unless $#views;
      $j .= ('?','xy','yx','yz','zy','xz','zx')[$i] . "-view" if $#views;
      $j = &get_input("$j: [$plotstyle{'file'}] ");
      $plotstyle{'file'} = $j unless $j =~ /^\s*$/;
      $j = $plotstyle{'file'} if $j =~ /^\s*$/;
      if(/^\s*\|(.*)$/) {
        (print "Could not pipe to '$1': $!\n\n"),next unless open(PRT,$j);
      } else {
        $j = &expand_file($j);
        (print "Could not open file '$j': $!\n\n"),next unless open(PRT,">$j");
      }
      select(PRT);
      last;
    }
    @bounds = @pbounds;
    @pgsz = @pagesize;  @marg = @margins;
    (@pgsz = reverse @pgsz),@marg=reverse @marg if $land == 2;

    $i--;
    @viewpt = split(/,/,('0,0,1','0,0,-1','1,0,0','-1,0,0',
                         '0,1,0','0,-1,0')[$i]);
    ($x,$y) = ((int $i/2)==1?1:0, (int $i/2)==0?1:2);
    ($x,$y) = ($y,$x) if (int $i)%2;

    # take maximum image wrt box in x direction.
    @bounds[$x,$x+3] = (&min($bounds[$x],-$box[$x]/2),
                        &max($bounds[$x+3],$box[$x]/2)) if $plotstyle{'box'};

    # take minimum with box in y direction if user says so, or if
    # <= 50% larger and box being drawn
    if($plotstyle{'box'}) {
      ($k,$j) = (0, $box[$y]/($bounds[$y+3] - $bounds[$y]));
      $k = 1 if $j <= 1.5;
      if(!$k && &sub_getopt('boxdraw') eq 'query') {
        select(STDOUT);
        $j = sprintf("%.4g",$j);
#        $_ = &get_input(" Picture would be scaled by ~1/$j if" .
#                        " the entire box is drawn.\n Draw anyway? [n] ");
#        $k = (/^y/i)?1:0;
	$k=1;
        select(PRT);
      }
      @bounds[$y,$y+3] = (&min($bounds[$y],-$box[$y]/2),
                          &max($bounds[$y+3],$box[$y]/2)) if $k;
    }

    if(!$plotstyle{'autobounds'}) {
      select(STDOUT);
      @foo = ();  $foo[$i+1] = 1;  $#foo = 6;
      @bounds = &_plt_getbounds(@foo,@bounds);
      select(PRT);
    }

    ($x,$y) = ($y,$x) if $land == 1;

    # reverse "x" bounds (about origin) if vertical landscape..
    @bounds[$x,$x+3] = (-$bounds[$x+3],-$bounds[$x]) if $land == '1';

    # make bbox a little (10%) larger than the minimum..
    foreach $j (0 .. 2) {$foo[$j] = ($bounds[$j+3]-$bounds[$j]) * 0.1/2;}
    @bounds = ($bounds[0]-$foo[0],$bounds[1]-$foo[1],$bounds[2]-$foo[2],
               $bounds[3]+$foo[0],$bounds[4]+$foo[1],$bounds[5]+$foo[2]);

    if($plotstyle{'title'}) {      # see below..
      $j = &sub_getopt('plotfontsz');
      if($j < 0.00001) {
        warn "Warning: Font size was '$j'!.  Title removed.\n";
        $plotstyle{'title'} = '';
      } else {
        $pgsz[0] -= $j;
      }
    }

    $scale = &min(
      ($pgsz[0] - 2*$marg[0])/($bounds[$x+3] - $bounds[$x]),
      ($pgsz[1] - 2*$marg[1])/($bounds[$y+3] - $bounds[$y]) );
    @trans = (($pgsz[0]/$scale - $bounds[$x+3] - $bounds[$x])/2,
              ($pgsz[1]/$scale - $bounds[$y+3] - $bounds[$y])/2 );


    if(&sub_getopt('bbox') eq 'large') {
      @bbox = (0,0,@pgsz);
    }
    elsif(&sub_getopt('bbox') eq 'margins') {
      @bbox = (@marg,$pgsz[0]-$marg[0],$pgsz[1]-$marg[1]);
    } else {
      @bbox = (($bounds[$x]+$trans[0])*$scale,($bounds[$y]+$trans[1])*$scale,
            ($bounds[$x+3]+$trans[0])*$scale,($bounds[$y+3]+$trans[1])*$scale);

      # if title, box enlarged to maximum between margins,
      # because we have no way of knowing how long it will be.
      # we also add the size of the font back on top..
      @bbox[0,2,3] = ($marg[0],$pgsz[0]-$marg[0],
                    $bbox[3]+&sub_getopt('plotfontsz')) if $plotstyle{'title'};
    }

# Set ends of bonds according to the sphere radii of the atoms...
    local($dist,$rad1,$rad2,@bondends);
    $#BONDENDS = 6 * $nbonds;
    foreach $i (0..$nbonds-1) {
      $k = $bonds[$i*$nbd];
      $l = $bonds[$i*$nbd + 1];
      $dist = $bonds[$i*$nbd + 2];
      $rad1 = $atm[$k*$nd + $aref{'radius'}]/$dist;
      $rad2 = $atm[$l*$nd + $aref{'radius'}]/$dist;
      foreach $j (0..2) {
        $bondends[$i*6 + $j] = $atm[$k*$nd +$j] + $rad1*($atm[$l*$nd+$j]
					- $atm[$k*$nd + $j]);
        $bondends[$i*6 + 3 + $j] = $atm[$l*$nd +$j] + $rad2*($atm[$k*$nd+$j]
					- $atm[$l*$nd + $j]);
      }
    }

    print "%!PS-Adobe-3.0 EPSF-3.0\n";
    printf "%%%%BoundingBox: %.5g %.5g %.5g %.5g\n",@bbox;
#    print "%%Creator: $progname\n%%CreationDate: ",&ctime(time);
    print "%%DocumentData: Clean7Bit\n";
    $sphere=<<ENDOFSPHERE;

/shadedrod
{ gsave
  /y2 exch def
  /x2 exch def
  /y1 exch def
  /x1 exch def
  /hd exch def
  /rodblue exch def
  /rodgreen exch def
  /rodred exch def
  /rodblue 1 rodblue sub def
  /rodgreen 1 rodgreen sub def
  /rodred 1 rodred sub def
  /cosb .5 def
  x1 y1 translate
  x2 x1 neg add
  y2 y1 neg add
  {atan neg rotate} stopped not {
  85 -5 0 {/angle exch def
  gsave
  newpath
   angle cos 1.0 rodred 0.5 mul neg add mul angle cos 1.0 rodgreen 0.5 mul neg add mul angle cos 1.0 rodblue 0.5 mul neg add mul setrgbcolor
%   nowred nowgreen nowblue setrgbcolor
%   /nowgrey angle cos 1.0 cosb 0.5 mul neg add mul 
%   /nowgrey 1 bangle cos sub cosb mul cosb .5 mul add def
%   cos 1.0 cosb 0.5 mul neg add mul setgray
%   nowgrey setgray
   angle sin 1.0 scale
   1 cosb scale
   0 0 hd 0 180 arcn
   x2 x1 neg add dup mul
   y2 y1 neg add dup mul
   add sqrt
   0 cosb eq {/cosb 1.0 def} if 0 exch cosb div translate
   0 0 hd 180 360 arc
  closepath fill
  grestore } for
  } if
  grestore } def

/latitude
{/lat exch def
newpath
0 lat sin r mul lat cos r mul 0 360 arc
%lat sin maxgray mul mingray add setgray
1 lat cos sub maxgray mul mingray add setgray
fill} def
/ball
{gsave
/tilt exch def
/r exch def
/maxgray exch def
/mingray exch def
/yy exch def
/xx exch def
xx yy translate
/mingray maxgray .5 mul def
gsave
300 rotate
newpath
0 0 r 0 360 arc
closepath
clip
mingray setgray
fill
1 tilt cos scale
-10 5 90 {latitude} for
grestore
newpath
%0 0 r 0 360 arc
%0 setgray stroke
%-xx -yy translate
grestore} def

ENDOFSPHERE

        $colorsphere=<<ENDOFCOLORSPHERE;

/colorlatitude
{/lat exch def
newpath
0 lat sin r mul lat cos r mul 0 360 arc
%lat sin maxgray mul mingray add setgray
/nowred 1 lat cos sub maxred mul maxred .5 mul add def
/nowblue 1 lat cos sub maxblue mul maxblue .5 mul add def
/nowgreen 1 lat cos sub maxgreen mul maxgreen .5 mul add def
nowred nowgreen nowblue setrgbcolor
fill} def
/colorball
{gsave
/tilt exch def
/r exch def
/maxblue exch def
/maxgreen exch def
/maxred exch def
/yy exch def
/xx exch def
xx yy translate
gsave
300 rotate
newpath
0 0 r 0 360 arc
closepath
clip
0 0 0 setrgbcolor
fill
1 tilt cos scale
-10 10 90 {colorlatitude} for
grestore
newpath
%0 0 r 0 360 arc
%0 setgray stroke
%-xx -yy translate
grestore} def

ENDOFCOLORSPHERE
    
    print $sphere;

print $colorsphere;

    if($plotstyle{'title'}) {
      $j = &sub_getopt('plotfont');
      print "%%DocumentFonts: $j\n%%DocumentNeededFonts: $j\n";
    }
    if($user) {
#      $foo = (getpwnam($user))[6];       ..getpwnam undefined on NeXT?!
#      $foo =~ s/,.*// if $foo;
#      print "%%For: $user" . ($foo?" ($foo)":'') . "\n";
       print "%%For: $user\n";
    }
    print "%%Orientation: ",($land==2?'Landscape':'Portrait'),"\n";
    print "%%Pages: 0\n%%Title: ($plotstyle{'title'})\n%%EndComments\n";

    $l = ($sym != 1) && $col != $ocol;
    print "/mv {moveto} bind def /sg {setgray} bind def /lc {2 copy} bind def";
    print "  /sg1 {$col sg} bind def" if $col >= 0;
    print "  /sg2 {$ocol sg} bind def" if $ocol >= 0;
    print " /sg3 {0 sg} bind def";
    print "\n";
    print " /X {2 copy} bind def" if $sym==2 && $l;
    print "/srgb {setrgbcolor} bind def\n" if $#rgbs>0;
    print "/X {2 copy moveto" . ($l?' 2 copy':'') . "} bind def\n" if $sym==3;
    if($sym == 1) {
      $j = sprintf("%.5g",&sub_getopt('Xsize')/2.5);
      $k = sprintf("%.5g",$j/2);
    print "/B  {exch $k add exch $k add 2 copy moveto -$j -$j rlineto stroke" .
        "\n     -$j add moveto -$j $j rlineto stroke} bind def\n";
    }
    elsif($sym == 2) {
      $j = sprintf("%.5g",&sub_getopt('diamonds')/2.5);
      if($j <= 0) {
        warn "Error: Diamond-size $j detected.  Set to 1 for this plot.\n";
        $j = 1;
      }
      print "/B {$j add mv -$j -$j rlineto $j -$j rlineto $j $j rlineto " .
            "closepath fill} bind def\n";
      print "/O {$j add mv -$j -$j rlineto $j -$j rlineto $j $j rlineto " .
            "closepath stroke} bind def\n";
    } else {
      if($rad < 0) {
	  print "/B {0 360 arc fill} bind def\n";
	  print "/O {3 copy add moveto 90 450 arc stroke} bind def\n";
	  print "/B2 {0 0 arc} bind def\n";
	  print "/O2 {ball} bind def\n";
	  print "/O2C {colorball} bind def\n";
      } else {
        print "/B {$rad 0 360 arc fill} bind def\n";
        print "/O {2 copy $rad add moveto $rad 90 450 arc stroke} bind def\n";
      }
    }
    $n = 1000;
    foreach $l (@grp) {
      $n= &min($n,$atm[$l*$nd+$aref{'radius'}]) if $atm[$l*$nd+$aref{'label'}];
    }
    $n = (int ($n*$scale*2/1.2))/2 if $n != 1000;
#   x y (s) -> x y (s) dx -> x y (s) x' -> x y (s) x' y'
    printf "/L {dup stringwidth pop 2 div 3 index exch sub 2 index %.4g sub\n"
         . "    moveto show pop pop} bind def\n",$n/$scale/4 if $n != 1000;
    $k = &sub_getopt('plotfont');
    if($j = $plotstyle{'title'}) {
      $l = &sub_getopt('plotfontsz');
      $j = &psstring($j);
      printf "0 sg /$k findfont %.4g scalefont setfont\n",$l;
      printf "%.4g ($j) stringwidth pop 2 div sub %.4g moveto ($j) show\n",
             ($bbox[0]+($bbox[2]-$bbox[0])/2),($bbox[3]-$l);
    }
    printf "/$k findfont %.4g scalefont setfont\n",$n/$scale if $n != 1000;
    printf "%.5g %.5g scale %.5g %.5g translate",$scale,$scale,@trans;
    $j = &sub_getopt('linewidth');
    print " $j setlinewidth";
    print " sg1" if $col>=0;
    print "\n";
    $code = "";
    $code .= 'foreach $j (sort viewpt_sort @grp) {$k=$j*' . "$nd;\n";
$code .= 'if ($atm[$k+'.$aref{'sphere'}.'] < 1) {';
    $code .= 'printf "%.4g %.4g' . ($l?' X':'') . ' ",' .
             ($land==1?'-':'+') . '$atm[$k+' . $x . '],$atm[$k+' . "$y];\n";
    $code .= "print 'lc ' if \$atm[\$k+$aref{'label'}];\n";
#    $code .= '$l=$n,printf("%.4g sg ",$n) if $l != ($n=$atm[$k+'
#           . $aref{'color'} . "]) && \$n>=0;\n" if $col<0;
    $code .= '$l=$n,printf("%.4g sg ",$n) if \($n=$atm[$k+'
            . $aref{'color'} ."])>=0;\n" if $col<0;
    $code .= 'printf("%.4g %.4g %.4g srgb ",@rgbs[$j*3,$j*3+1,$j*3+2])' .
             ' if $n<0;' if $#rgbs > 0;
    $code .= 'printf("%.4g ",$atm[$k+' . $aref{'radius'} . ']);'
              if $sym == 3 && $rad<0;
    $code .= "print 'B';";
$code.='}'."\n";
    if($l) {
$code .= 'if ($atm[$k+'.$aref{'sphere'}.'] < 1) {';
	    $code .= '$l=$n,printf(" %.4g sg ",$n) if $l != ($n=$atm[$k+'
		. $aref{'color'} . "]);\n" if $ocol<0;
	    $code .= 'print " sg2 ";' if $ocol>=0;
	    $code .= 'print "O";' if $sym!=3 || $rad>=0;
            $code .= 'printf("%.4g O",$atm[$k+'. $aref{'radius'} .']);'
		if $sym==3 && $rad<0;
	    $code .= "\nif(\$k = \$atm[\$k+$aref{'label'}]) {\n";
	    $code .= 'if($k=~/^\^([.\d]+)\^(.*)$/) {printf " %.4g sg ($2) L",$1} ' .
#               'else {$1=0;print " 0 sg ($k) L"} }' . "\n";
		'else {print " 0 sg ($k) L"} }' . "\n";
$code .= '} else {';
            $code .= 'printf("%.4g %.4g %.4g %.4g %.4g %.4g 45 O2C",'.($land==1?'-':'+') . '$atm[$k+' . $x . '],$atm[$k+' . $y.'],@rgbs[$j*3,$j*3+1,$j*3+2],'.'$atm[$k+'. $aref{'radius'} .']) if ($atm[$k+'.$aref{'color'}.'])<0;' if $#rgbs > 0 && $sym==3 && $rad<0;
	    $code .= 'printf("%.4g %.4g 0 %.4g %.4g 45 O2",'.($land==1?'-':'+') . '$atm[$k+' . $x . '],$atm[$k+' . $y.'],$atm[$k+'.$aref{'color'}.'],'.'$atm[$k+'. $aref{'radius'} .']) if ($atm[$k+'.$aref{'color'}.'])>=0;'
		if $col <0 && $sym==3 && $rad<0;
	    $code .= "\n\$kk=\$k;if(\$k = \$atm[\$k+$aref{'label'}]) {\n";
	    $code .= 'printf " %.4g %.4g ",'.($land==1?'-':'+') . '$atm[$kk+' . $x . '],$atm[$kk+' . $y.']; ';
	    $code .= 'if($k=~/^\^([.\d]+)\^(.*)$/) {printf " %.4g sg ($2) L",$1} ' .
		'else {printf "  0 sg ($k) L"} }' . "\n";
$code .='}';
#      $code .= "\nif(\$k = \$atm[\$k+$aref{'label'}]) {\n";
#      $code .= 'if($k=~/^\^([.\d]+)\^(.*)$/) {printf " %.4g sg ($2) L",$1} ' .
#               'else {$1=0;print " 0 sg ($k) L"} }' . "\n";
#               'else {print " 0 sg ($k) L"} }' . "\n";
      $code .= 'print " sg1 "; $l ' . "= $col;\n" if $col>=0;
      $code .= "\$l = \$k?\$1:$ocol;\n" if $col<0;
} else {
      $code .= '$k = $atm[$k+$aref{\'label\'}]; ';
      $code .= 'if($k) {if($k=~/^\^([.\d]+)\^(.*)$/) {printf ' .
               '"%.4g sg ($2) L",$1} else {$1=0;print "0 sg ($k) L"}}' . "\n";
#               '"%.4g sg ($2) L",$1} else {print "0 sg ($k) L"}}' . "\n";
      $code .= 'print "sg1" if $k;' . "\n" if $col>=0;
      $code .= '$l = $1 if $k;' . "\n" if $col<0;
}
    $code .= 'print "\n";' . "\n";
# Print out bonds of this atom to higher atoms first:
# Added APS June 2, 1994
    $code .= 'foreach $i (0..$nbonds-1) { ' . "\n";
    $code .= 'if ((($bonds[$i*$nbd] == $j)&&' .
	'(&viewpt_diff($bonds[$i*$nbd+1],$j) >= 0)) ||' .
	    '(($bonds[$i*$nbd+1] == $j)&&' .
		'(&viewpt_diff($bonds[$i*$nbd],$j) > 0))) {';
    $code .= "\n";

# line bonds
    $code .= 'if ($bonds[$i*$nbd+'.$bref{'shade'}.']!=1) {';
    $code .= 'printf(" %.4g setlinewidth ",$bonds[$i*$nbd+'.$bref{'width'}.']);';
    $code .= 'if (($bond_dash) '. 
	'&& ($bonds[$i*$nbd + $bref{\'length\'}] > $dash_len)) {'.
	    'printf("[0.05 .2] 0 setdash\n");}';
    $code .= 'printf("sg3 %.4g %.4g moveto ", ' . ($land==1? '-':'+') .
			' $bondends[$x + 6*$i],' .
			' $bondends[$y + 6*$i]);' . "\n";
    $code .= 'printf("%.4g %.4g lineto\n", ' . ($land==1?'-':'+') .
			' $bondends[$x + 6*$i + 3],' .
			' $bondends[$y + 6*$i + 3]);' . "\n";
    $code .= 'print  "closepath stroke\n";';
    if ($bond_dash) {$code .= 'printf("[] 0 setdash\n");';}
    $j = &sub_getopt('linewidth');
    $code.= 'print "'. $j.' setlinewidth "';

# shaded bonds
   
    $code.=' ;} else { ';
    $code.= 'printf("%.4g %.4g %.4g %.4g %.4g %.4g %.4g %.4g shadedrod ",$bonds[$i*$nbd+'.$bref{'red'}.'],$bonds[$i*$nbd+'.$bref{'green'}.'],$bonds[$i*$nbd+'.$bref{'blue'}.'],$bonds[$i*$nbd+'.$bref{'width'}.'],' . ($land==1? '-':'+') .
			' $bondends[$x + 6*$i],' .
			' $bondends[$y + 6*$i],'. ($land==1?'-':'+') .
			' $bondends[$x + 6*$i + 3],' .
			' $bondends[$y + 6*$i + 3]);' . "\n";
    $code.= ";}}}}\n";
#    print STDOUT "Code is:\n-------------------\n$code-------------------\n";
    eval($code);
    if($@) {
      warn "Error evaluating printing instructions:\n  $@\n";
      warn "\nCode was:\n$code\n";
      return;
    }
    if($plotstyle{'box'}) {
      print  "0 sg ";
      printf("%.4g %.4g moveto ", -$box[$x]/2,-$box[$y]/2);
      printf("%.4g %.4g lineto\n",-$box[$x]/2, $box[$y]/2);
      printf("%.4g %.4g lineto ",  $box[$x]/2, $box[$y]/2);
      printf("%.4g %.4g lineto ",  $box[$x]/2,-$box[$y]/2);
      print  "closepath stroke\n";
      $j = &min($box[$x]/10,10/$scale);
      printf("%.4g 0 moveto %.4g 0 rlineto stroke\n",-$box[$x]/2, $j);
      printf("%.4g 0 moveto %.4g 0 rlineto stroke\n", $box[$x]/2,-$j);
      $j = &min($box[$y]/10,10/$scale);
      printf("0 %.4g moveto 0 %.4g rlineto stroke\n",-$box[$y]/2, $j);
      printf("0 %.4g moveto 0 %.4g rlineto stroke\n", $box[$y]/2,-$j);
    }
    print "showpage\n";  # Importers are supposed to annul _showpage_.
    select(STDOUT);  close(PRT);
  }
}
sub _plt_getbounds {
  local(@views,@bounds,@foo,@xyz,$j) = @_;
  @bounds = splice(@views,7,6);
  $xyz[0] = ($views[0] || $views[1] || $views[4] || $views[5])?1:0;
  $xyz[1] = ($views[0] || $views[1] || $views[2] || $views[3])?1:0;
  $xyz[2] = ($views[2] || $views[3] || $views[4] || $views[5])?1:0;
  (print "No bounds to edit!\n"),return unless grep($_,@xyz);
  while(1) {
    grep($_ = sprintf("%.5g",$_), @bounds);
    @foo = ();
    foreach $i (0..2) {
      push(@foo,sprintf("% 5g <= %s <= %g",$bounds[$i],('x','y','z')[$i],
                $bounds[$i+3])) if $xyz[$i];
    }
    print 'Bounds are {',join(",\n            ",@foo),"}\n";
    $_ = &get_input("Change these bounds? [y] ");
    last if /^\s*n/i;
    foreach $i (0..2) {
      next unless $xyz[$i];
      $j = ('x','y','z')[$i];
      $_ =&get_input("\nEnter ${j}min, ${j}max [$bounds[$i],$bounds[$i+3]]: ");
      @bounds[$i,$i+3] = &getnums(2) unless /^\s*$/;
      @bounds[$i,$i+3] = (&min(@bounds[$i,$i+3]),&max(@bounds[$i,$i+3]))
    }
  }
  @bounds;
}
sub sub_showhelp {
  local($p) = 0;
  if(&sub_getopt('pager') eq 'off') {
    if(open(PGR, "| $pager") || open(PGR, "| more")) {
      select(PGR); $p=1;
    }
  }
  print <<"EOH";
 General Idea:

  This program allows one to select - by position, tag, or component -
  groups of atoms,  change  these attributes plus radius,color,etc.  ,
  and save or plot the results.

  'labels' are what  go in the center  of atoms when they are plotted.
  They are by default blank.  'sel info' lists attributes of each atom
  in the current  selection, including its label  (if any).   They are
  drawn in black by default, but if you precede the label with ^grey^,
  'grey' (e.g. 1 for white) will be used instead of black (0).

 Commands not listed on menu:

  alias [name [value]]:
    Displays all aliases if <name> not given.
    Displays value of alias <name> if value not given.
    Otherwise,   sets alias for  <name> to  <value> (preserved between
       sessions)
    References   to other  aliases inside of     an alias will  not be
       resolved.

  unalias name - removes the alias for 'name'.

  menu - forces display of menu (even when 'setopt menu off', below).

  setopt [option [value]]:
    options are:
     bbox - plots have minimal bounding  box ('minimal'), one  bounded
                  by the  default margins (~1  inch) ('margins', or by
                  the full-page  bounding box ('large') (all have same
                  size image)
     boxdraw - 'query' if program  should query  user about  including
                  all of box if drawn   (see [*] below).  'silent'  by
                  default.
     diamonds -  scaling for size  of diamonds (plotting mark). '1' by
		  default.
     landscape - true landscape ('horizontal') or landscape rotated to
                  'vertical'.  'vertical' by default.
     menu - display menu? one of (on,off).
     pager - pager type: one of (on,off,less,more).  'off' by default.
                  if 'on', uses the PAGER environment variable.
     plotfont - name   of postscript font plot   uses  - Helvetica  by
                  default.
     plotfontsz   - size of   postscript font  plot uses   - 24pt.  by
                  default.
     linewidth - width of lines drawn for bonds and edges.
     postcom - a command to be evaluated after every command input.
     prompt - prompt type: one of (plain,dir).  'plain' by default.
     save - if 'renumbered', atoms  are  renumbered  in the  outputted
                  .con  file   (so    that   they are    consecutive).
                  Otherwise, they maintain  the same numbering (in the
                  last column) as in the original file.
     selmode - mode selections use (add,int,set,sub,xor) if one is not
                  given.  'set' by default.
     Xsize - scaling for size of X's (plotting mark). '1' by default.

   If <option>  is  not given,  values  of all options  are displayed.
   Option-settings (and  other   information) are  preserved   between
   sessions in the file ~/$datfile.

  Version - displays brief information on this version of $progname.


 Features not  noted  in  menu:  

  The  Interrupt signal (INT,  usually   generated by ^C)  aborts  the
  current function.  If you aren't in a function  (i.e.  you're at the
  main   prompt), the program quits  (without   saving defaults -  see
  below).

  The  'undo'  command, when  listed,  undoes  the effect of  the last
  change of the selection, and/or change in atomic positions.

  If a command is not recognized, it will be run via /bin/sh (e.g. ls,
  date,more,etc.).

  Commands can be semicolon-separated (e.g. sel list; save sel).

  If  a command   prompts for data,  you  can  give  the data  on  the
  command-line immediately after  the  command.  If multiple  data are
  prompted for, they can  be  seperated with  '/' on the  command-line
  (excess data are ignored).  E.g. 'sel save sph8; sel int sph 0,0,0 /
  8; save sel out.con /n'

  Commands (sph,box,boxp,rot,scale) which prompt  for an atom or point
  for the center will accept 0 for the origin (as opposed to regarding
  the single number as an atom: there  is no atom  0).  Thus, 'sel int
  sph 0,0,0 / 8' above could be 'sel int sph 0/8'.

  'color <sel>'  accepts  a plus  sign (+) or   a minus sign (-)  as a
  character before the color to indicate a relative change.

  'color  <sel>' accepts 'by depth' on  the command-line, allowing one
  to  color  atoms by  a gradient  in  x,y,or z.  This   is useful for
  perceiving depth in pictures.  (+ and - are allowed as above)

  'set rgbcolor <sel>' will  let you color  atoms  in color.   It will
  prompt for Red,Green,and Blue components, and accepts the 'by depth'
  flag   on the  command-line  (as  above).   To   set atoms  back  to
  greyscale, just use 'color' (or 'set color').

  'set radius' accepts 'touch' (or 'touching')  instead of a radius as
  instruction  to calculate the largest radius  such that no two atoms
  extend into eachother (this may be time-consuming to calculate).

  'set  radius'  also accepts  depth,  allowing one to  set the radius
  according to  the  depth along  the chosen  axis.  As with  color by
  depth, you  can reverse the  ordering of the  extrema  to get larger
  atoms further  away.   Unlike color by  depth,  you must  wait until
  prompted for a radius before entering depth.

  'sel info' takes two optional, command-line-only options:
    -v         : enables velocity display (if applicable)
    an integer : limits display to at most this many entries.

  'sel list'  takes  an optional, command-line-only option  indicating
  the  number of  indices to print  per  line (default: 10). e.g. 'sel
  list 20'

  'sel save' preserves saved selections between sessions.

  If   the file  '.$rcfile'  is present  in your   home directory, the
  program will act as if the  user is typing  the contents of .$rcfile
  in.  (The author has 'setopt menu off; menu' in his)

  In the 'plot' command:

   Given one chooses the view 'rs' with the  other coordinate being t,
   one will  get a plot of  the selection towards  the origin in the t
   direction such that the r coordinate increases as one goes right on
   the screen,  and the t coordinate  increases as one  goes up on the
   screen.  (so that yx looks from below the configuration, and zy (in
   landscape mode) is the same as yz after a 180-degree rotation about
   the z axis)

   If one chooses automatic  bounding,  only the transverse  (wrt  the
   plot)  coordinate  of the  box  will be taken  into account, unless
   including the entire vertical  face would shrink  the plot by  less
   than 50% [*].   e.g., when a  surface is present,  and one does not
   want to see the entire z length of the box.  If  one did, one could
   plot zy instead of yz (if one  didn't want the picture rotated, one
   could do zy as a landscape plot)

   [*] If the box  isn't being drawn, it isn't  taken  into account at
      all.  If it is, and you setopt 'boxdraw' to 'query', you will be
      asked whether you want to  include the vertical coordinate if it
      would shrink the plot by >%50.

  This program is written in Perl, a scripting language by Larry Wall.
  Hence you can edit the source..

EOH
  select(STDOUT),close(PGR) if $p;
}

sub show_version {
  print " $progname version $version,\n written by Daniel Faken",
        "  ~  (email) absinthe\@u.washington.edu\n\n";
}

sub expand_file {  # tilde-expansion
  local($f,$g) = @_;
  if($f =~ /^~(\w+)/) {
    (warn "Couldn't get pwd info on user $f: $!\n"),return
                                   unless $g = (getpwnam($1))[7];
    $f =~ s/^~$1/$g/;
  } else {
    $f =~ s/^~/$home/;
  }
  $f;
}

sub psstring {
  local(@s) = @_;
  foreach $i (@s) {
    $i =~ s/\\/\\\\/g;
    $i =~ s/\(/\\\(/g;
    $i =~ s/\)/\\\)/g;
  }
  wantarray?@s:$s[0];
}

sub load_defaults {
  local($f,$str) = @_;
  return unless -e "$home/$f";
  open(DAT,"$home/$f") ||
         ((warn "Error opening datafile '$home/$f' for reading: $!\n"),return);
  $str = join('',<DAT>);
  eval($str);
  close(DAT);
# Added APS June 7, 1994: In case data file from previous tcon version.
  &_merge_defaults();
}
sub _merge_defaults {
  local($i,@opts);
  foreach $i (keys(%dsettings)) {
    @opts = split(',',$dsettings{$i});
# Fake a user call to setopt
    $#FN_ARGS = 2;
    $FN_ARGS[0] = $i;
    $FN_ARGS[1] = $opts[0];
    &sub_setopts();
  }
}
sub _save_assocar {
  local($name,*assoc,$j,$k) = @_;
  print "%$name = (";
  foreach $i (sort alphabetically keys(%assoc)) {
    $k = $assoc{$i};
    $k =~ s/'/\\'/g;
    print ',' if $j++;  print "'$i','$k'";
  }
  print ");\n";
}
sub save_defaults {
  local($f) = @_;
  open(DAT,">$home/$f") ||
                ((warn "Error opening '$home/$f' for writing: $!\n"),return);
  select(DAT);
  print "# DO NOT CHANGE THIS FILE!  $progname maintains it.\n";
# Changed APS June 7, 1994:
  &_save_assocar('dsettings',*settings);
  print "\n";
  &_save_assocar('plotstyle',*plotstyle);
  print "\n";
  &_save_assocar('aliases',*aliases);
  print "\n";
  &_save_assocar('saved',*saved);
  print "\n# chksum: 2398139\n";
  select(STDOUT);
  close(DAT);
}
