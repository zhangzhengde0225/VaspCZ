#!/usr/bin/env perl
#;-*- Perl -*-
use FindBin qw($Bin);
use lib "$Bin";
use Vasp;
#####################################################################################
print "This update will read in a config file and add 4 parameter flags for diminit.pl.\n";
$configfile="config";
if(-e $configfile){
  open (CONFIG,"<$configfile") or die "Can't open $config file.\n";
  $content="";
  while($line=<CONFIG>){
    $line=~s/^\s+//;
    if($line eq "\n"){next;}
    @line=split(/\s+/,$line);
    if($line[1] eq '=') {
      $content.=$line[0]." = ".$line[2]."\n";
      $parms->{lc($line[0])}=$line[2];
    }
  }
  close CONFIG;
  $DisplaceAlgo="DisplaceAlgo";
  $DisplaceRange="DisplaceRange";
  $NN_rcut = "NN_rcut";
  $MaxCoordNum = "MaxCoordNum";
  print "flags to be added (if not already there): $DisplaceAlgo $DisplaceRange $NN_rcut $MaxCoordNum\n";
  $i=0;
  if(!exists($parms->{lc($DisplaceAlgo)})){$content.=$DisplaceAlgo." = 1\n";$i++;}
  if(!exists($parms->{lc($DisplaceRange)})){$content.=$DisplaceRange." = 3\n";$i++;}
  if(!exists($parms->{lc($NN_rcut)})){$content.=$NN_rcut." = 2.6\n";$i++;}
  if(!exists($parms->{lc($MaxCoordNum)})){$content.=$MaxCoordNum." = 8\n";$i++;}
  print "***********************************\n";
  if($i){
    print $content;
    print "***********************************\n";
    print "$i flags added to the end of the $configfile file\n";
    open (CONFIG,">$configfile") or die "Can't open $config file.\n";
    print CONFIG $content;
    close CONFIG;
  }else{
    print "the $configfile file was already updated. no changes made.\n";
  }
}else{
  die "please run it in the akmc directory where the config file resides.\n"
}


