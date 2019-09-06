#!/usr/bin/env perl
#;-*- Perl -*-

use FindBin qw($Bin);
use lib "$Bin";
use Vasp;

#################################################################################
# This script resets the "done" and "bad" saddles in the st.dat file to the status
# of "quenching" and "promising" so that the saddle evaluation can be performed
# with a new energy and distance tolerance. There may be other uses.
# (written by Lijun Xu and Graeme Henkelman in the University of Texas at Austin.
# last modifed by Lijun on July 20, 2007 )
#################################################################################
@ARGV==1 || die "Usage: akmc_reset_bad.pl <st.dat> \n";
$min="min";
$min1="min1";
$min2="min2";
$dynmat="dynmat";
$quench="mins";
$ststatus="status";
$numprst="numprst";
$filename=$ARGV[0];
$output="st_reset_bad.dat";
open(OUTPUT, ">$output") || die "cannnot open $output.\n";
print "************************************\n";
print "Read the current st data ......";
($header,$st)=ReadSt($filename);
print "done\n";
for($i=1;$i<=$st->{$numprst}[0];$i++) {
    $prdir=DirName("pr",$i);
    if(lc($st->{$prdir}[0]) eq "done"){
      if(lc($st->{$prdir}[1]) eq "bad"){
        print "*************************************************\n";
        print "Reset $prdir which was done and bad ... ";
        $initialst="$prdir/$quench/$min1/final/CONTCAR";
        $finalst="$prdir/$quench/$min2/final/CONTCAR";
        if(-e $initialst && -e $finalst){
          $st->{$prdir}[0]="quench";
          $st->{$prdir}[1]="promising";
          $st->{$prdir}[3]="na";
          $st->{$prdir}[4]="na";
        }
        print "done\n";
      }
    }#end of the outmost if
}
print "Save the new st data in $output ......";
WriteSt($header, $st, $output);
print "done\n";
# ------------------------------------
# Directory name
# ------------------------------------
sub DirName {
  $pre=shift;
  $num=shift;
  $snum=sprintf "%04d",$num;
  $dirname=$pre.$snum;
  return $dirname;
}

# -----------------------------------------------------------------------
# Read state information from st.dat [%st=ReadSt($StFile)]
# -----------------------------------------------------------------------
sub ReadSt{
  my $StFile=shift;
  my $stfilename=$StFile;
  my %st=();
  my $line="";
  my @line=();
  my ($i, $j, $numprst, $header);
  %st=();
  if(-e $stfilename){
    open(ST,"<$stfilename");
  }else{
    $st{"status"}=[("NOST")];
    $st{"numprst"}=[(0)];
    $st{"NumSearchesLeft"}=[("unknown")];
    return(\%st);
  }
  $i=0;
  while($line=<ST>){
    $i++;
    $line=~s/^\s+//;
    if($line eq "\n"){ $i--; next; } # skip empty lines
    @line=split(/\s+/,$line);
    if($i==1){# the first line must be st and st energy etc.
      $header=$line[0];
      $st{$header}=[()]; 
      for ($j=1;$j<@line;$j++) {push @{$st{$header}}, $line[$j];}
    }elsif($i==2){
      $st{"status"}=[()];
      for ($j=1;$j<@line;$j++) {push @{$st{"status"}}, lc($line[$j]);}
    }elsif($i==3){
      $st{"dynmat"}=[()];
      for ($j=1;$j<@line;$j++) {push @{$st{"dynmat"}}, lc($line[$j]);} 
    }elsif($i==4){
      $st{"NumSearchesLeft"}=[()];
      for ($j=1;$j<@line;$j++) {push @{$st{"NumSearchesLeft"}}, lc($line[$j]);} 
    }elsif($i==5){
      $st{"numprst"}=[()];
      for ($j=1;$j<@line;$j++) {push @{$st{"numprst"}}, lc($line[$j]);} 
    }else{
      if($i==6 && lc($line[0]) eq "process"){ #an optional description line
        $st{"prinfolist"}=[()];
        for ($j=0;$j<@line;$j++) {push @{$st{"prinfolist"}}, lc($line[$j]);}
        $i--;
        next;
      }
      $numprst=$i-5;
      if (!(lc($line[0]) eq DirName("pr",$numprst))){die "in $StFile, st00xx has a wrong format.\n";}
      $st{lc($line[0])}=[()];
      for ($j=1;$j<@line;$j++) {push @{$st{lc($line[0])}}, lc($line[$j]);}
      if($st{lc($line[0])}[2]=~/^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$/) {
        $st{lc($line[0])}[2] = $st{lc($line[0])}[2] + $st{$header}[0]; # convert barriers into energy.
      }
    }
  }
  if($st{"numprst"}[0] != $numprst){die "The numbers of st00xx states conflict in $stfilename\n";}
  close ST;
  if(exists($st{""})) {delete $st{""}};
  return($header,\%st);
}

# --------------------------------------------------------------------
# Write state information to st.dat [WriteSt($header,$st,$stfilename)]
# --------------------------------------------------------------------
sub WriteSt{
  my ($curdir,$st,$stfilename)=@_;
  my ($energy, $force);
  my $i=0;
  my $j=0; 
  open (OUT,">$stfilename") || die "In WriteSt: cannot open $stfilename\n";
  if(exists($st->{$curdir})){
    $energy=$st->{$curdir}[0];
    if($curdir ne substr($curdir, -6)) {
      print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
      print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
      print "!!!! Warning In WriteSt: bad $curdir information.\n";
      print "!!!! It is supposed to have the st000x (6 letters/numbers) format.\n";
      print "!!!! Make sure it is what you wanted!\n";
      print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
      print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";      
    }
    print OUT $curdir." ".$energy."\n";
  }else{die "No $curdir information in the st hash. Check it.\n";}
  if(exists($st->{"status"})){
    print OUT "status"." ".$st->{"status"}[0]."\n";
  }else{die "No status information in the st hash. Check it.\n";}
  if(exists($st->{"dynmat"})){
    print OUT "dynmat"." ".$st->{"dynmat"}[0]."\n";
  }else{die "No dynmat information in the st hash. Check it.\n";}
  if(exists($st->{"NumSearchesLeft"})){
    print OUT "NumSearchesLeft"." ".$st->{"NumSearchesLeft"}[0]."\n";
  }else{die "No NumSearchesLeft information in the st hash. Check it.\n";}
  if(exists($st->{"numprst"})){
    print OUT "numprst"." ".$st->{"numprst"}[0]."\n";
  }else{die "No numprst information in the st hash. Check it.\n";}
  if(exists($st->{"prinfolist"})){ # discription line is not required
    print OUT "@{$st->{'prinfolist'}}\n";
  }#else{die "No prinfolist information at all in the st hash. Check it.\n";}
  for $i (sort keys %$st){
    if($i eq $curdir || lc($i) eq "status" || lc($i) eq "numprst" || lc($i) eq "dynmat" || lc($i) eq "numsearchesleft" || lc($i) eq "prinfolist"){
      next;
    }
    if($st->{$i}[2]=~/^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$/) {
      $dummy=$st->{$i}[2];
      $st->{$i}[2] = $st->{$i}[2] - $st->{$curdir}[0];
      print OUT $i." "."@{$st->{$i}}"."\n";
      $st->{$i}[2] = $dummy;
    }else{
      print OUT $i." "."@{$st->{$i}}"."\n";
    }
    $j++;
  }
  if($j != $st->{"numprst"}[0]) {die "In WriteST: Recorded number of state processes is not right.\n";}
  close OUT;
}

