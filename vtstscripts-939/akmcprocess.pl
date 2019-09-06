#!/usr/bin/env perl
#;-*- Perl -*-

# 08-31-2008: GH added $Bin

use FindBin qw($Bin);
use lib "$Bin";
use Vasp;

chomp($prdir=`pwd`);
$prdir=~s/^\s+//;
@prdir=split(/\//,$prdir);
$prdir=$prdir[-1];
print "Try to a make movie for the process in $prdir ... ...\n";

$mins="mins";
$min1="min1";
$min2="min2";
$movie=$prdir.".xyz";
$StFile="st.dat";
$quadruple=0;
if(@ARGV >0) {
  $quadruple=$ARGV[0];
}
if(@ARGV >1){
  $StFile=$ARGV[1];
}
if(!(-e "../$StFile")) {die "You must run this script inside the process directory: e.g. cd pr0012; akmc_process.pl <1 or 0 (quadruple)> <st.dat>. Note: 0 is the default, so is st.dat\n";}
$st=ReadSt("..",$StFile);
#@www=(keys %$st);
#print "@www\n";
$POSCAR_sp="POSCAR_sp";
if(exists($st->{$prdir}) && (lc($st->{$prdir}[0]) eq "done" || lc($st->{$prdir}[0]) eq "killed") && lc($st->{$prdir}[1]) eq "good") {
  $final=$st->{$prdir}[3];
  if($final eq "mins/min1") {
    $initial="mins/min2";
  }elsif($final eq "mins/min2") {
    $initial="mins/min1";
  }else{
    die "fatal error: $final is unrecognizable.\n";
  }
}elsif(exists($st->{$prdir}) && lc($st->{$prdir}[0]) eq "done" && lc($st->{$prdir}[1]) eq "bad"){
  $initial="mins/min1";
  $final="mins/min2";
}else{
  die "we have no idea about the status of this process. do nothing.\n";
}
if(-e "$initial/medgamma") {
    if($quadruple) {
      $dummy=`cd $initial; $Bin/quad.pl POSCAR out; $Bin/pos2con.pl out; $Bin/con2xyz.pl out.con; mv out.xyz poscar_init.xyz; rm out out.con;`;
    }else{
      $dummy=`cd $initial; $Bin/pos2con.pl POSCAR; $Bin/con2xyz.pl POSCAR.con; rm POSCAR.con; mv POSCAR.xyz poscar_init.xyz`;
    }
}elsif(-e "$initial/CONTCAR") {
    if($quadruple) {
      $dummy=`cd $initial; $Bin/quad.pl CONTCAR out; $Bin/pos2con.pl out; $Bin/con2xyz.pl out.con; mv out.xyz poscar_init.xyz; rm out out.con;`;
    }else{
      $dummy=`cd $initial; $Bin/pos2con.pl CONTCAR; $Bin/con2xyz.pl CONTCAR.con; rm CONTCAR.con; mv CONTCAR.xyz poscar_init.xyz`;
    }
}else{
    die "fatal error: $initial is not correctly set up\n";
}
if(-e "$final/medgamma") {
    if($quadruple) {
      $dummy=`cd $final; $Bin/quad.pl POSCAR out; $Bin/pos2con.pl out; $Bin/con2xyz.pl out.con; mv out.xyz poscar_final.xyz; rm out out.con;`;
    }else{
      $dummy=`cd $final; $Bin/pos2con.pl POSCAR; $Bin/con2xyz.pl POSCAR.con; rm POSCAR.con; mv POSCAR.xyz poscar_final.xyz`;
    }
}elsif(-e "$final/CONTCAR") {
    if($quadruple) {
      $dummy=`cd $final; $Bin/quad.pl CONTCAR out; $Bin/pos2con.pl out; $Bin/con2xyz.pl out.con; mv out.xyz poscar_final.xyz; rm out out.con;`;
    }else{
      $dummy=`cd $final; $Bin/pos2con.pl CONTCAR; $Bin/con2xyz.pl CONTCAR.con; rm CONTCAR.con; mv CONTCAR.xyz poscar_final.xyz`;
    }
}else{
    die "fatal error: $final is not correctly set up\n";
}
if(-e "POSCAR_sp") {
    if($quadruple) {
      $dummy=`$Bin/quad.pl POSCAR_sp out; $Bin/pos2con.pl out; $Bin/con2xyz.pl out.con; mv out.xyz POSCAR_sp.xyz; rm out out.con;`;
    }else{
      $dummy=`$Bin/pos2con.pl POSCAR_sp; $Bin/con2xyz.pl POSCAR_sp.con; rm POSCAR_sp.con`;
    }
}else{
    die "No saddle $POSCAR_sp file.\n";
}
system("cat $initial/poscar_init.xyz POSCAR_sp.xyz $final/poscar_final.xyz > $movie; rm $initial/poscar_init.xyz POSCAR_sp.xyz $final/poscar_final.xyz");
print "... ... done\n";

sub DirName{
  my $pre=shift;
  my $num=shift;
  my $snum="";
  my $dirname="";
  $snum=sprintf "%04d",$num;
  $dirname=$pre.$snum;
  return $dirname;
}

sub ReadSt{
  my $curdir=shift;
  my $StFile=shift;
  my $stfilename=$curdir."/".$StFile;
  my %st=();
  my $line="";
  my @line=();
  my ($i, $j, $repeat,$numprst);
  %st=();
  if(-e $stfilename){
    open(ST,"<$stfilename");
  }else{
    $st{"status"}=[("NOST")];
    $st{"numprst"}=[(0)];
    $st{"NumSearchesLeft"}=[(0)];
    return(\%st);
  }
  $i=0;
  while($line=<ST>){
    $i++;
    $line=~s/^\s+//;
    @line=split(/\s+/,$line);
    if($i==1){
      $st{$curdir}=[()]; 
      for ($j=1;$j<@line;$j++) {push @{$st{$curdir}}, $line[$j];}
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
    }else{ # Read process information
      $numprst=$i-5;
      if (!(lc($line[0]) eq DirName("pr",$numprst))){die "in $StFile, st00xx has a wrong format.\n";}
      $st{lc($line[0])}=[()];
      for ($j=1;$j<@line;$j++) {push @{$st{lc($line[0])}}, lc($line[$j]);}
    }
  }
  if($st{"numprst"}[0] != $numprst){die "The numbers of st00xx states conflict in $stfilename\n";}
  close ST;
  if(exists($st{""})) {delete $st{""}};
  return(\%st);
}

