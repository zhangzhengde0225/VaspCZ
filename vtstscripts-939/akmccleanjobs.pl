#!/usr/bin/env perl
#;-*- Perl -*-
#########################################################
# This patch cleans up killed jobs in akmc. It is not a crtical patch.
# jobs that were killed due to high barriers were left intact in the old akmc.pl
# It cleans up the killed run into "killed" folder to save some disk space.
# In the updated akmc.pl, the killed runs will be handled as unconverged and 
# cleaned into inter* folders.
# last updated by Lijun Xu at Oak Ridge 06/20/2008
#########################################################
use FindBin qw($Bin);
use lib "$Bin";
$maindir=".";
if(@ARGV>0){$maindir=$ARGV[0];}
opendir MAINDIR, $maindir or die "can't open this $maindir dir!";
@InMainDir=readdir MAINDIR;
@stdirs = grep /^st\d{4}$/, @InMainDir;
@stdirs=(sort @stdirs);
closedir MAINDIR;
for $stdir (@stdirs){
  if(!(-e "$stdir/st.dat")){ 
     print "=====================================\n";
     print "=========WARNING !!!!!!==============\n";
     print "no st.dat found in the current folder. will do nothing\n";
     print "=====================================\n";
  }else{
    @killedprs=`cd $stdir; grep killed st.dat | sort`;
    for $prdir (@killedprs){
      $prdir=~s/^\s+//;
      @dummy=split(/\s+/, $prdir);
      $target=$stdir."/".$dummy[0];
      if(lc($dummy[1]) eq "killed" && (-e $target)){
       if( -e "$target/killed"){
         print "folder $target already existed. we will skip cleanup here\n";
       }else{
         if(-e "$target/OUTCAR"){
           print "---------clean folder---$target------\n";
           system("cd $target; $Bin/vfin.pl killed");
         }else{
           print "No OUTCAR in $target. We asssume it has been cleaned up\n";
         }
       }
      }
    }
  }
}
