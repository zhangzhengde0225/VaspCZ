#!/usr/bin/env perl
#;-*- Perl -*-

@args=@ARGV;
@args==2 || die "usage: chgparavg.pl <reference CHGCAR> <CHGCAR2>\n";


open (IN1,$args[0]) || die ("Can't open file $!");
open (IN2,$args[1]) || die ("Can't open file $!");
open (OUT,">PARCHG_avg");

for ($i=0;$i<2;$i++) {
    $line1=<IN1>;
    $line2=<IN2>;
    $header1.=$line1;
}

for ($i=0;$i<3;$i++) {
    $line1=<IN1>;
    $line1 =~ s/^\s+//;
    @line1=split(/\s+/,$line1);
    $vector[$i]=$line1[$i+1];
    $line2=<IN2>;
    $line2 =~ s/^\s+//;
    $header1.=$line1;
}

print "Note: periodic averaging of ionic coordinates assumes orthogonal cell.\n";
print "Vectors: ".$vector[0].", ".$vector[1].", ".$vector[2]."\n";

$atoms1=<IN1>;
$header1.=$atoms1;
$atoms2=<IN2>;

$atoms1 =~ s/^\s+//;
$atoms2 =~ s/^\s+//;
@atoms1=split(/\s+/,$atoms1);
@atoms2=split(/\s+/,$atoms2);

$sum1 += $_ for @atoms1;
$sum2 += $_ for @atoms2;

print "Atoms in file1: ".$sum1.", Atoms in file2: ".$sum2."\n";

$header1.=<IN1>;
$dummy=<IN2>;

for ($i=0;$i<$sum1;$i++) {
   $coordinates1.=<IN1>;
}
for ($i=0;$i<$sum2;$i++) {
   $coordinates2.=<IN2>;
}

@coordinates1=split(/\n/,$coordinates1);
@coordinates2=split(/\n/,$coordinates2);

for ($i=0;$i<@coordinates1;$i++) {
   @line1=split(/\s+/,$coordinates1[$i]);
   @line2=split(/\s+/,$coordinates2[$i]);
   for ($j=1;$j<4;$j++) {
      $delta=$line2[$j]-$line1[$j];
      if ($delta>0.5) {
         $delta=$delta-1;
      } elsif ($delta<-0.5) {
         $delta=$delta+1;
      }
      print "$i $delta\n";
      $line1[$j]=$line1[$j]+$delta/2;
      if ($line1[$j]>1) {
         $line1[$j]=$line1[$j]-1;
      } elsif ($line1[$j]<0) {
         $line1[$j]=$line1[$j]+1;
      }
   }
   $coordinates1[$i]=join(" ",@line1);
}

$coordinates1=join("\n",@coordinates1);

$dummy=<IN1>;
$dummy=<IN2>;

$header1.=$coordinates1."\n";
$points1=<IN1>;
$header1.=$points1;
$points2=<IN2>;

$points1 =~ s/^\s+//;
$points2 =~ s/^\s+//;
@points1 = split(/\s+/,$points1);
@points2 = split(/\s+/,$points2);

$psum1 = $points1[0]*$points1[1]*$points1[2];
$psum2 = $points2[0]*$points2[1]*$points2[2];

print "Points in file1: ".$psum1.", Points in file2: ".$psum2."\n";

if ($psum1 != $psum2) {die ("Number of points not same in two files!");}

print OUT $header1;

for ($i=0;$i<$psum1/10;$i++) {
    $line1=<IN1>;
    $line1 =~ s/^\s+//;
    $line2=<IN2>;
    $line2 =~ s/^\s+//;
    @line1=split(/\s+/,$line1);
    @line2=split(/\s+/,$line2);
    for ($j=0;$j<@line1;$j++) {
         $line1[$j]=($line2[$j]+$line1[$j])/2;
    }
#    printf OUT "%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n",$line1[0],$line1[1],$line1[2],$line1[3],$line1[4],$line1[5],$line1[6],$line1[7],$line1[8],$line1[9];  
    printf OUT " %20.10E" x @line1 . "\n", @line1;

}

close(IN1);
close(IN2);
close(OUT);
