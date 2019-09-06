#!/usr/bin/env perl
#;-*- Perl -*-

#####################################################################################
# Adaptive kinetic Monte Carlo Simulation 2.0 (integrated with VASP and vtsttools)  #
# by Lijun Xu and Graeme Henkelman (The University of Texas at Austin)              #
# (last modified by Lijun on Sept. 30, 2011)                                        #
# version 1.0: two-image dimer method                                               #
# version 2.0: single-image dimer method replaces the  two-image method             #
# version 3.0 (proposed): add LANCZOS (single image) to the code                    #
#                                                                                   #
# Algorithms are described in the following publications:			    #
# L. Xu and G. Henkelman, J. Chem. Phys. 129, 114104 (2008)                         #
# G. Henkelman and H. Jonsson, J. Chem. Phys. 115, 9657-9666 (2001)  		    #
#####################################################################################
use FindBin qw($Bin);
use lib "$Bin";
use Vasp;

if($/ ne "\n"){
  print "Warning: input record separator(".$/.") is not the default (new line). Reset it.\n";
  $/="\n"; #now, chomp will remove the last "\n" in a string.
}
# Note: $Bin will be a global variable and thus do not use it for other variables   #
#$dummy=substr($Bin, -1, 1);
#print "$dummy\n";
#if($dummy eq "\/") {chop($Bin);}
#print "The \$Bin is ".$Bin."\n";
print "Script directory is ".$Bin."\n";
# Note: %UPDATEPROCESS is a global HASH holding information for one process         #
# Set the random seed
srand();

###############################################################################
# Read config file
###############################################################################
$configfile = "config";
print "Read $configfile file:\n";
print "===========================================================\n";
$MaxJobs            = ReadParm($configfile, "MaxJobs",          0);
$RandSeed           = ReadParm($configfile, "RandSeed",         0);
$NumSearches        = ReadParm($configfile, "NumSearches",      20);
$AkmcSteps          = ReadParm($configfile, "AkmcSteps",        -1);
$SimR               = ReadParm($configfile, "SimR",             0.1); #RT More descriptive variable name.
$Ediffmax           = ReadParm($configfile, "Ediffmax",         0.05);
$Rdiffmax           = ReadParm($configfile, "Rdiffmax",         0.04);
$NumKT              = ReadParm($configfile, "NumKT",            20);
$Temperature        = ReadParm($configfile, "Temperature",      77);
$SearchesAlgo       = ReadParm($configfile, "SearchesAlgo",     "Dimer");
$BarrierMax         = ReadParm($configfile, "BarrierMax",       10.0);
$JobFile            = ReadParm($configfile, "JobFile",          "jobs.dat");
$StFile             = ReadParm($configfile, "StFile",           "st.dat");
$AkmcFile           = ReadParm($configfile, "AkmcFile",         "akmc.dat");
$DynMatFile         = ReadParm($configfile, "DynMatFile",       "freq.dat");
$Prefactor          = ReadParm($configfile, "Prefactor",        0);
$RateTableFile      = ReadParm($configfile, "RateTableFile",    "RateTableFile");
$OldRunDirFile      = ReadParm($configfile, "OldRunDirFile",    "OldRunDirFile");
$StEnergyFile       = ReadParm($configfile, "StEnergyFile",     "StEnergyFile");
$Population         = ReadParm($configfile, "Population",       0);
$GrabQuenched       = ReadParm($configfile, "GrabQuenched",     0);
$PrRecycleAll       = ReadParm($configfile, "PrRecycleAll",     1);
$PrRecycle          = ReadParm($configfile, "PrRecycle",        0);
$PrRecycleShift     = ReadParm($configfile, "PrRecycleShift",   0);
$SurfaceRecShift    = ReadParm($configfile, "SurfaceRecShift",  0);
$ConvergenceAlgo    = ReadParm($configfile, "ConvergenceAlgo",  1);
$Equivalency        = ReadParm($configfile, "Equivalency",      0);
$Screenoutput       = ReadParm($configfile, "ScreenOutput",     "out.dat");
$DisplaceAlgo       = ReadParm($configfile, "DisplaceAlgo",     1);
$DisplaceRange      = ReadParm($configfile, "DisplaceRange",    3);
$NN_rcut            = ReadParm($configfile, "NN_rcut",          2.6);
$MaxCoordNum        = ReadParm($configfile, "MaxCoordNum",      8);
$BoltzmanEqu        = ReadParm($configfile, "BoltzmanEqu",      0.999999); 
$UseKDB             = 0; #ReadParm($configfile, "UseKDB",           0);
$KDBcutoff          = ReadParm($configfile, "KDBcutoff",        0.1);
$WaitForQuenchDyn   = ReadParm($configfile, "WaitForQuenchDyn", 0);
###############################################################################
# Check for required parameters
###############################################################################
($MaxJobs >= 0) || die "Error: MaxJobs should be a positive or zero number\n";
($NumSearches > 0 && $NumSearches < 10000) || die " NumSearches must be less than 10000.\n";
($SimR>0) || die "Error: negative SimR, the distance that the dimer images are pushed away to initialize minimizations\n";
if($Screenoutput) {
  $Screenoutput_tmp=$Screenoutput."_tmp";
  $Screenoutput_old=$Screenoutput."_old";
  if(-e $Screenoutput) {system "mv $Screenoutput $Screenoutput_old";}
}
print "MaxJobs = $MaxJobs, SearchesAlgo = $SearchesAlgo.\n";

$turnoff_dynmat=0;
if($Prefactor) {$turnoff_dynmat=1;}
#if($turnoff_dynmat) {print "Skip dynamical matrix calculations by assuming a prefactor of $Prefactor 1/s.\n";}
if($turnoff_dynmat) {print "Assuming a prefactor of $Prefactor/s (skip calculation).\n";}
$min="min";
$dynmat="dynmat";
$quench="mins";
$ststatus="status";
$RateTableFileBolEqu=$RateTableFile."_BolEqu";
$enforced="";
$time_saved=0.0;
if($UseKDB){print "\nUsing kdbquery results, if they exist.\n";}
###############################################################################
# Read job information.
###############################################################################
print "Reading old job information: ";
$jobs=ReadJobs($JobFile);
#for $j (sort keys %$jobs){ print "$j: @{$jobs->{$j}}\n"; }
###############################################################################
# Check and update running jobs
###############################################################################
($numrunjobs,$jobs)=CheckJobs($jobs);
print "Found $numrunjobs running jobs.\n";
#for $j (sort keys %$jobs){ print "$j: @{$jobs->{$j}}\n"; }
###############################################################################
# Get current state and state pool information
###############################################################################
print "Read state pool information: ";
$stpool=ReadStPool($OldRunDirFile);
if(@$stpool>0){
  print "done.\n";
  print "include old states in folders: ";
  for ($j=0;$j< @$stpool;$j++){ print "$stpool->[$j] ";}
  print "\n";
}else{
  print "No old states available.\n";
}
print "===========================================================\n";
print "Number of states: ";
#Note: we keep . (i.e., akmc directory) as the default directory
#RT Going to be changed...
($good,$stdirs,$numst,$curst)=CheckSystemConfig(".");
print $numst.". ";
if($curst){ 
  print "Current: $curst. Configuration is ";
  if($good) { $CurrentStateIsLast=1;
  }else{
    print "not ";
    $CurrentStateIsLast=0;
  }
#  print "the same as the one at home.\n";
  print "good.\n";  #GH: I don't understand what this means
}else{
  print "\n";
}
###############################################################################
# Read akmc step information.
###############################################################################
print "Reading $AkmcFile: ";
$akmc_step = ReadAkmcStep($AkmcFile);
$step = 0;
@lasttwosteps = (sort {$a <=> $b} keys %$akmc_step);#RT No sorting required here. Need only read the last line.
if(@lasttwosteps == 0){print "NONE. It hasn't started yet\n";
}else{print "Step number ";}
for $j (@lasttwosteps){ # there should be only one member i.e., the last step
  $step = $j;
  #print "step $j: @{$akmc_step->{$j}}\n";
  #print "$j. Info: state ".${$akmc_step->{$j}}[0].", ".${$akmc_step->{$j}}[1].", ".${$akmc_step->{$j}}[2].", energy ";
  print "$j, state ".${$akmc_step->{$j}}[0].", ".${$akmc_step->{$j}}[1].", energy ";
  printf ("%3.3f",${$akmc_step->{$j}}[3]);
  print " eV.\n";
}

if($step == 0) { #no $AkmcFile, i.e., minimize the initial configuration
  print "Initialize the first configuration ...\n";
  ($InitDone, $numrunjobs, $jobs, $akmc_step) = Initializer($min, $numrunjobs, $MaxJobs, $jobs, $akmc_step);
  if($InitDone) {
    #copy minimized initial configuration to akmc root directory
    system "cp $min/POSCAR POSCAR"; #Note: min must have been cleaned up.
    #akmc_step hash should have a key "1" and a list value starting from "2bcreated"
    $curstep = "1";
    $Done = 0;
  }else{
    $curstep = "0";
    $Done = 1;
  }
}else{
  # RT: Minimization is done.
  # Get the current step based on the previous step information read from the AkmcFile
  if($curst eq $akmc_step->{$step}[2]) { # RT: Check the state number instead.
    # either first or last since the previous state is same as the current state
    $cur_eq_pre=1;
    if($CurrentStateIsLast) {
      $prest = "";
      $curstep = $step;
      print "We are at the first step.\n";
    }else{
      $prest = $akmc_step->{$step}[2];
      $curstep = $step + 1;
      $akmc_step->{$curstep} = [("2bcreated","na","na",$akmc_step->{$step}[4],"na","na","na",$akmc_step->{$step}[7])];
      print "We are trying to continue from a finished or stopped akmc run.\n";
    }
  }else{ # Create the akmc_step for the current step
    $cur_eq_pre = 0;
    if($CurrentStateIsLast) { # We are working on an unfinished state
      $prest = $akmc_step->{$step}[2];
      $curstep = $step+1;
      $akmc_step->{$curstep} = [($curst, "new", $curst, $akmc_step->{$step}[4], "na", "na", "na", $akmc_step->{$step}[7])];
      print "We are working on an unfinished state.\n";
    }else{
      $prest = $akmc_step->{$step}[2];
      $curstep = $step + 1;
      $akmc_step->{$curstep}=[("2bcreated","na","na",$akmc_step->{$step}[4],"na","na","na",$akmc_step->{$step}[7])];
      print "It stopped short of making a new st. We will do it in this run.\n";
    }
  }
  $Done = 0; # RT: Maybe move this to the top.
}
###############################################################################
# Check if enough steps have been done.
###############################################################################
if($curstep > $AkmcSteps) {
  print "===========================================================\n";
  print "However, we already had maximum(=$AkmcSteps) steps. Do nothing until AkmcSteps in $configfile is increased.\n";
  $Done = 1;
}
print "===========================================================\n";
###############################################################################
# Main loop
###############################################################################
while(!$Done){ # RT: Rename $Done to something signifying that the Main loop is done. (e.g., $MainLoopDone)
  #---------------------------------------------------------------------------
  # Get the current state
  #---------------------------------------------------------------------------
  print "Current step is: $curstep. ";
  $curdir = $akmc_step->{$curstep}[0];
  $repeat = 0;
  if($curdir eq "2bcreated"){ # working on the home poscar
    # check the current state pool to see if this is a new state.
    ($repeat,$akmc_step) = CheckStpool($curstep,$stpool,$StEnergyFile,$Equivalency,$Ediffmax,$Rdiffmax,$akmc_step);
    if(!$repeat) {
      print "Current (new) state is: ";
      $numst++;
      $curdir = DirName("st", $numst);
      ($curdir eq $akmc_step->{$curstep}[0]) || die "error: akmc_step{currentstep}[0] is not same as curdir\n";
      MakeNewSt($curdir);
      $akmc_step->{$curstep}[2] = $curdir;
      print "$curdir\n";
    }else{
      # RT: Getting state from State Pool.
      $curdir = $akmc_step->{$curstep}[2];
      print " current st(repeated): same as $akmc_step->{$curstep}[0] in $akmc_step->{$curstep}[2]".".\n";
    }
  }elsif($akmc_step->{$curstep}[1] eq "new" && $curdir eq DirName("st", $numst)){ # working on the latest state
    if(-e $StEnergyFile){ # RT Redo the format of this file.  Rename "States.dat" or somesuch.
      chomp($lxtwj = `grep $curdir $StEnergyFile | tail -1`);
      if($lxtwj ne "") {$akmc_step->{$curstep}[1] = "repeat";}
    }
    if($akmc_step->{$curstep}[1] eq "repeat") {
      $repeat=1;
      print "Current state (repeated) is: same as $akmc_step->{$curstep}[0] in $akmc_step->{$curstep}[2]".".\n";
    }else{
      print "Current state (latest) is: $curdir.\n";
    }
  }else{
    die "Fatal error in KMC $curstep: Unrecognizable st folder name $curdir.\n";
  }
  #---------------------------------------------------------------------------
  # Read current state data
  #---------------------------------------------------------------------------
  $st = ReadSt($curdir, $StFile);
  $ststatus = $st->{"status"}[0];
  $numprst = $st->{"numprst"}[0];
  if($repeat) {
    if($ststatus eq "done") {
      if(!$turnoff_dynmat) {
        if(!(-e "$curdir/$dynmat")) { die "Warning: we need dynmat calculations which are missing in the repeated state $curdir. Please run $Bin/akmcdynmat.pl for it.\n;"}
      }
    }else{
      die "We are dealing with a repeated st ($curdir), so it must have been marked as done.\n";
    }
  }
  #---------------------------------------------------------------------------
  # Update current state information
  #---------------------------------------------------------------------------
  $curstenergy = $akmc_step->{$curstep}[3];
  if(exists($st->{$curdir})) {
    if(abs($curstenergy - $st->{$curdir}[0]) > 0.01 ) { # 0.01 ev is the threshold, normally they should be almost exactly same.
      die "The states are messed up. Check the next state energy in the previous state akmc_step\nand the config at home which should be the current state configuration.\n";
    }
  }
  if($ststatus eq "NOST"){
    print "State $curdir is new with no processes.  Create a new state and quit.\n";
    # prepare sp and dynmat jobs and submit them.
    ($jobs, $st) = StartSt($curdir,$MaxJobs,$NumSearches,$SearchesAlgo,$quench,$dynmat,$turnoff_dynmat,$curstep,$stpool,$PrRecycleAll,$PrRecycle,$UseKDB,$KDBcutoff,$PrRecycleShift,$SurfaceRecShift,$StEnergyFile,$StFile,$Equivalency,$Ediffmax,$Rdiffmax,$DisplaceAlgo,$DisplaceRange,$NN_rcut,$MaxCoordNum,$akmc_step,$jobs,$st); 
    $numprrun = $st->{"numprst"}[0];
    $numprfin = 0;
    $numrunjobs = $MaxJobs;
    if($curstep == 1) {WriteAkmcStep($AkmcFile, $akmc_step, 1);}
    $Done = 1;
  }elsif($ststatus eq "running"){
    print "State $curdir is $ststatus with $numprst processes and ".$st->{"NumSearchesLeft"}[0]." searches left.\n";
    # save a copy of the old state information
    $st_old = StHashCopy($st);
    ($unique_sp,$total_sp,$good_saddles_inprocess,$belowTenKT,$belowTenKT_total,$MinSaddle,$degeneracy,$sp_energy) = CountSP($curdir,$SearchesAlgo,$NumKT,$Temperature,$st_old);
    if($MinSaddle){
      $MaxSaddle=$st->{$MinSaddle}[2] + $NumKT * $Temperature / 11604.5;
    }else{ # no good saddle has been found
      $MaxSaddle="";
    }
    # check the state
    $numprrun=$numprfin=$numprkilled=0;
    $status_change=1;
    if(!$turnoff_dynmat) {
      # update dynmat calculation
      print "===========================================================\n";
      print "Update dynamical matrix for $curdir: ";
      $curdyndir=$curdir."/".$dynmat;
      if(lc($st->{$dynmat}[0]) eq "done") {
        $numprfin++;
        $status_change=0;
      }else{
        if(!(exists($jobs->{$curdyndir}) && exists($st->{$dynmat}))) {
          die "There is one dynmat directory in $curdir without any job or state information.\n";
        }
        print "HandleDynmat ... ";
        ($numrunjobs,$jobs,$comment)=HandleDynmat($curdyndir,$DynMatFile,$numrunjobs,$MaxJobs,$jobs);
        # determine st information from job information
        if(lc($jobs->{$curdyndir}[0]) eq "completed") {
          $st->{$dynmat}[0]="done";
          $numprfin++;
        }elsif(lc($jobs->{$curdyndir}[0]) eq "2bsubmitted") {
          $st->{$dynmat}[0]="2bsubmitted";
          $numprrun++;
        }elsif(lc($jobs->{$curdyndir}[0]) eq "killed") {
          $st->{$dynmat}[0]="killed";
          $numprkilled++;
        }else{
          $st->{$dynmat}[0]="running"; # after handling the state is either "running" or "done" or "killed".
          $numprrun++;
        }
      }
      print "now it is $st->{$dynmat}[0]\n";
    }else{# do not calculate dynamical matrix
      $status_change=0;
      #$numprfin++;      
    }
    if($status_change) {
      #print "Save the updated jobs($JobFile) and st($curdir/$StFile) information ... ...";
      WriteJobs($JobFile,$jobs);
      WriteSt($curdir,$StFile,$st);
      #print "done\n";
    }
    # update processes
    print "Update process table for $curdir:\n";
    print "===========================================================\n";
    printUPDATEPROCESS_header($Screenoutput,$Screenoutput_tmp);
    $SPEnergyMax=$st->{$curdir}[0] + $BarrierMax;  # the maximum saddle point energy allowed.
    $nowait=1;
    for($i=1; $i<=$numprst; $i++){
      if($st->{"NumSearchesLeft"}[0]==0) {# state is running with no searches
       if($WaitForQuenchDyn){# wait for quenching and dynmat calculations to be finished
         $nowait=0;
         #print "It looks like enough saddle points have been found, but there may be a few unknown saddles.";
         #print "We may have to halt searches until the current quench and dynmat processes are finished.\n";
       }else{
         print "Stop: the convergence criteria for this state has been met.\n";
         print "Waiting for the dynamical matrix (if any) for the state configuration to finish.\n"; 
         last;
       }
      }
      $prdir=DirName("pr",$i);
      $curspdir=$curdir."/".$prdir;
      $curspdyndir=$curspdir."/".$dynmat;
      $curspquendir=$curspdir."/".$quench;
      #update process status
      $status_change=1;
      #print "***********************\n";
      #print "In process $prdir: ";
      PrepareUPDATEPROCESS(); # UPDATEPROCESS is a global hash for holding updating information for one process
      if(!(exists($st->{$prdir}))){
        die "There is a $prdir directory in $curdir without any state information.\n";
      }
      $prstatus=lc($st->{$prdir}[0]);
      $UPDATEPROCESS{"was"}=$prstatus;
      if($prstatus eq "done"){
        $status_change=0;
        if(!(exists($jobs->{$curspdir}))){ # one time issue. make the job list complete
          $jobs = PatchupRecycledJobInfo($jobs,$curspdir,$curspdyndir,$curspquendir,$quench,$dynmat,$turnoff_dynmat,$SearchesAlgo);
          # Assuming that it is a recycled(retrievd) saddle, add dummy for missing job info.
        }
        #print " $prdir is already done.\n";
        #print "***********************\n";
        $numprfin++;
        $UPDATEPROCESS{"now"}="";
        $UPDATEPROCESS{"energy"}=$st->{$prdir}[2];
        $UPDATEPROCESS{"comment"}=$st->{$prdir}[1];
        if(lc($st->{$prdir}[1])=~/^repeat\*/){
          @repeat_prdir=split(/\*/, $st->{$prdir}[1]);
          $UPDATEPROCESS{"comment"}=$UPDATEPROCESS{"comment"}."(".$st->{$repeat_prdir[1]}[1].")";
        }elsif(lc($st->{$prdir}[1]) eq "bad"){
          if(-e $curspquendir){$UPDATEPROCESS{"comment"}="Bad mins";}
          else{$UPDATEPROCESS{"comment"}="Negative Barrier";}
        }
      }elsif($prstatus eq "dynmat"){
        #print " $prdir was doing dynmat calculation.\n";
        #print "***********************\n";
        #print "updating dynmat...\n";
        ($numrunjobs,$jobs,$comment)=HandleDynmat($curspdyndir,$DynMatFile,$numrunjobs,$MaxJobs,$jobs);
        #print "... now the dynmat job is ";
        # determine state information from job information
        if(lc($jobs->{$curspdyndir}[0]) eq "completed"){
          $st->{$prdir}[0]="done";
          $numprfin++;
          #print "done\n";
        }elsif(lc($jobs->{$curspdyndir}[0]) eq "2bsubmitted"){
          $numprrun++;
          #print "2bsubmitted\n";
          $UPDATEPROCESS{"jobstatus"}="2bsubmitted";
        }elsif(lc($jobs->{$curspdyndir}[0]) eq "killed"){
          $st->{$prdir}[0]="killed";
          $st->{$prdir}[1]=$comment."_dynmat";
          $numprkilled++;
          #print "killed\n";
        }else{
          #print "still $jobs->{$curspdyndir}[0] in $curdyndir\n";
          $numprrun++;
          $UPDATEPROCESS{"jobstatus"}=$jobs->{$curspdyndir}[0];
        }
        $UPDATEPROCESS{"now"}=$st->{$prdir}[0];
        $UPDATEPROCESS{"comment"}=$comment;        
      }elsif($prstatus eq "quench"){
        #print " $prdir was doing quench.\n";
        #print "***********************\n";
        #print "updating quench...\n";
        $spenergy=$st->{$prdir}[2];
        ($numrunjobs,$quenchstatus,$jobs,$comment)=HandleQuench($curspquendir,$SearchesAlgo,$numrunjobs,$MaxJobs,$SimR,$Ediffmax,$Rdiffmax,$spenergy,$jobs);
        #print "... now it is ";
        if($quenchstatus eq "completed"){
          #print "completed.\n";
          %minsinfo=("$quench/min1"=>$jobs->{$curspquendir."/min1"}[3], "$quench/min2"=>$jobs->{$curspquendir."/min2"}[3]);
          $st=EvalSP($curdir,$prdir,$quench,$SearchesAlgo,\%minsinfo,$Ediffmax,$Rdiffmax,$st);
          if($st->{$prdir}[1] eq "bad"){
            #print "but the saddle is not good. this $prdir is halted\n";
            $st->{$prdir}[0]="done";
            $numprfin++;
            if($UseKDB){ KdbAddPr($curspdir); }
            $UPDATEPROCESS{"comment"}="BadMins"; 
          }elsif($st->{$prdir}[1] eq "good"){
            #print "and the saddle is good.";
            if(!$turnoff_dynmat){
              #print " So, start a dynmat job ... ";
              ($numrunjobs,$jobs)=StartDynmat($curspdyndir,$numrunjobs,$MaxJobs,$jobs);
              $st->{$prdir}[0]="dynmat";
              $UPDATEPROCESS{"jobstatus"}=$jobs->{$curspdyndir}[0];
              if(lc($jobs->{$curspdyndir}[0]) eq "killed"){
                $st->{$prdir}[0]="killed";
                $st->{$prdir}[1]=$jobs->{$curspdyndir}[1]."_dynmat";
                $UPDATEPROCESS{"comment"}=$jobs->{$curspdyndir}[1];
              }
              $numprrun++;
              #print "done.\n";
            }else{
              #print " And no dynmat job is needed. done\n";
              $st->{$prdir}[0]="done";
              $numprfin++;
            }
            if($UseKDB){ KdbAddPr($curspdir); }
          }else{
            die "In $curspquendir, unreasonable $st->{$prdir}[0] after quench.\n"
          }
        }elsif($quenchstatus eq "pushmore"){
          #print "needs be re-started by pushing two images further more\n";
          delete($jobs->{$curspquendir."/min1"});
          delete($jobs->{$curspquendir."/min2"});
          ($numrunjobs,$jobs,$quenchstatus2,$comment)=StartQuench($curspdir,$SearchesAlgo,$numrunjobs,$MaxJobs,2*$SimR,$jobs);
          $UPDATEPROCESS{"jobstatus"}="run";
          if($quenchstatus2 eq "killed"){
            $st->{$prdir}[0]="killed";
            $st->{$prdir}[1]=$comment."_quench";
            $UPDATEPROCESS{"jobstatus"}="---"; # killed quench has no status
          }else{
            $UPDATEPROCESS{"comment"}="min1=min2=>restart";
          }
          $numprrun++;          
        }elsif($quenchstatus eq "running"){
          #print "still running(or to be submitted).\n";
          $numprrun++;
          $UPDATEPROCESS{"jobstatus"}="run";
        }elsif($quenchstatus eq "killed"){
          $st->{$prdir}[0]="killed";
          $st->{$prdir}[1]=$comment."_quench";
          $numprkilled++;
          #print "killed.\n";
          #print "Warning from HandleQuench: quench stopped fatally due to (at least) one job that was killed by code or user.\n";
        }else{
          die "In $curspquendir, unreasonable $st->{$prdir}[0] after quench.\n"
        }
        $UPDATEPROCESS{"now"}=$st->{$prdir}[0]; 
      }elsif($prstatus eq "search"){
        #print " $prdir was doing search.\n";
        #print "***********************\n";
        #print "updating ...\n";
        ($numrunjobs,$jobs,$comment)=HandleSP($curdir,$prdir,$SearchesAlgo,$SPEnergyMax,$numrunjobs,$MaxJobs,$nowait,$jobs);
        #print " ... now it is\n";
        if(lc($jobs->{$curspdir}[0]) eq "completed"){
          #print "completed.\n";
          $st->{$prdir}[2] = $jobs->{$curspdir}[3]; # get the saddle point energy;
	  ($recycled, $st) = FilterSP($curdir, $prdir, $numprst, $BarrierMax, $SearchesAlgo,
					$quench, $dynmat, $StFile, $turnoff_dynmat, $Equivalency,
					$stpool, $GrabQuenched, $Ediffmax, $Rdiffmax,
					$DynMatFile, $stdirs, $st);
          @repeat_prdir=split(/\*/, $st->{$prdir}[1]);
          if(lc($repeat_prdir[0]) eq "repeat"){
            #print "but it turns out to be repeating process ".$repeat_prdir[1]." The process is halted.\n";
            $st->{$prdir}[0]="done";
            $numprfin++;
          }elsif(lc($st->{$prdir}[1]) eq "promising"){
            #print "it is a promising saddle. Start a quench ...";
            ($numrunjobs,$jobs,$quenchstatus2,$comment)=StartQuench($curspdir,$SearchesAlgo,$numrunjobs,$MaxJobs,$SimR,$jobs);
            $st->{$prdir}[0]="quench";
            $UPDATEPROCESS{"jobstatus"}="run";
            if($quenchstatus2 eq "killed"){
              $st->{$prdir}[0]="killed";
              $st->{$prdir}[1]=$comment."_quench";
              $UPDATEPROCESS{"jobstatus"}="---"; # killed quench has no status
            }elsif($quenchstatus2 eq "2bchecked"){
              $comment="2bchecked";
            }else{
              $comment="a promising saddle";
            }
            $numprrun++;
            #print "done\n";
          }elsif(lc($st->{$prdir}[1]) eq "bad"){
            #print "but it turns out to be a bad saddle due to bad mins or negative barrier. The process is halted.\n";
            if($recycled){#fill the job information for possible re-evaluation
              $jobs->{$curspquendir."/min1"}=[("stopped","jobid_marker","minimization","na","na")];
              $jobs->{$curspquendir."/min2"}=[("stopped","jobid_marker","minimization","na","na")];
            }
            $st->{$prdir}[0]="done";
            $numprfin++;
            $comment="BadSaddle(previously quenched)";
          }elsif(lc($st->{$prdir}[1]) eq "highenergy"){
            #print "but it turns out to be a high energy saddle. The process is halted.\n";
            $st->{$prdir}[0]="done";
            $numprfin++;
            $comment="HighEnergySaddle";
          }elsif(lc($st->{$prdir}[1]) eq "good"){
            #print "it is a good saddle in this state, based on the quenched identical sp in previous states.";
            if($recycled){#fill the job information for possible re-evaluation
              $jobs->{$curspquendir."/min1"}=[("stopped","jobid_marker","minimization","na","na")];
              $jobs->{$curspquendir."/min2"}=[("stopped","jobid_marker","minimization","na","na")];
            }
            if(lc($st->{$prdir}[0]) eq "done") {
              if($recycled && !$turnoff_dynmat){
                $jobs->{$curspdyndir}=[("stopped","jobid_marker","dynmat","na","na")];
              }
              $numprfin++;
              #print "Now, $prdir is done.\n";
              $comment="GoodSaddle(previously quenched)";
            }elsif(lc($st->{$prdir}[0]) eq "dynmat") { # this is set only if the dynmat is turned on
              #print "Now, start a dynmat job ... \n";
              ($numrunjobs,$jobs)=StartDynmat($curspdyndir,$numrunjobs,$MaxJobs,$jobs);
              $UPDATEPROCESS{"jobstatus"}=$jobs->{$curspdyndir}[0];
              if(lc($jobs->{$curspdyndir}[0]) eq "killed"){
                $st->{$prdir}[0]="killed";
                $st->{$prdir}[1]=$jobs->{$curspdyndir}[1]."_dynmat";
                $comment=$jobs->{$curspdyndir}[1];
              }
              $numprrun++;
              #print "done.\n";
            }elsif(lc($st->{$prdir}[0]) eq "killed") {
              $numprkilled++;
              #print "Now, $prdir is killed.\n"; 
              $comment="Killed(previously quenched)";
            }else{
              die "In $curspdir, unreasonable st status: $st->{$prdir}[0] for a good saddle.\n";
            }            
          }else{
            die "In $curspdir, unreasonable saddle quality: $st->{$prdir}[1] after filtering the saddle.\n";
          }
        }elsif(lc($jobs->{$curspdir}[0]) eq "killed") { # either by code or manually
          #print "killed by the user or by the code: bad curvature/barrier/broken job, etc.\n";
          $st->{$prdir}[0]="killed";
          $st->{$prdir}[1]=$comment."_search";
          $numprkilled++;
        }elsif(lc($jobs->{$curspdir}[0]) eq "2bsubmitted") {
          #print "2bsubmitted\n";
          $UPDATEPROCESS{"jobstatus"}=$jobs->{$curspdir}[0];
        }else{
          #print "still running.\n";
          $numprrun++;
          $UPDATEPROCESS{"jobstatus"}=$jobs->{$curspdir}[0];
        }
        $UPDATEPROCESS{"now"}=$st->{$prdir}[0];
        $UPDATEPROCESS{"comment"}=$comment;
      }elsif($prstatus eq "killed"){
        #print "$prdir has been killed (such as broken jobs/crazy saddle energy/disallowed/killed by the user).\n";
        $numprkilled++;
        $UPDATEPROCESS{"now"}="";
        $UPDATEPROCESS{"comment"}=$st->{$prdir}[1]; 
      #}elsif($prstatus eq "2breassessed"){# to be removed. No change st status into "search pending" to be re-assessed
      #  #print "$prdir is to be re-evaluated: update the barrier, etc information.\n";
      #  $st=ReEvaluateProcess($curdir,$prdir,$quench,$SearchesAlgo,$turnoff_dynmat,$Ediffmax,$Rdiffmax,$jobs,$st);
      }else{
        die "Unrecognizable status value: $prstatus for process $prdir.\n";
      }
      if($status_change) {
        #print "Save the updated jobs($JobFile) ... ...";
        WriteJobs($JobFile,$jobs);
        if(lc($st->{$prdir}[0]) ne lc($st_old->{$prdir}[0])) { 
          #print "and st($curdir/$StFile) ... ...";
          WriteSt($curdir,$StFile,$st);
        }
        #print "done\n";
      }
      if($ConvergenceAlgo == 1 || $ConvergenceAlgo == 0){ # check the state status at every step
        if($status_change && (lc($st->{$prdir}[0]) ne lc($st_old->{$prdir}[0]))) { # process status is different
         ($unique_sp,$total_sp,$good_saddles_inprocess,$belowTenKT,$belowTenKT_total,$MinSaddle,$MaxSaddle,$degeneracy,$sp_energy,$st)=SaddleStat($curdir,$prdir,
          $SearchesAlgo,$ConvergenceAlgo,$NumKT,$Temperature,$NumSearches,$unique_sp,$total_sp,$good_saddles_inprocess,$belowTenKT,$belowTenKT_total,$MinSaddle,$MaxSaddle,$degeneracy,$sp_energy,$st);
         #print "NumSearchesLeft number at this step is ".$st->{"NumSearchesLeft"}[0]."\n";
         if($st->{"NumSearchesLeft"}[0] == 0){
           $UPDATEPROCESS{"NumSLeft"}=0;
           printUPDATEPROCESS($prdir,$curdir,$curspdyndir,$curspquendir,$st,$jobs,$quench,$dynmat,$turnoff_dynmat,$SearchesAlgo,$Screenoutput,$Screenoutput_tmp);
           if($WaitForQuenchDyn){ $nowait=0;
           }else{
             $nowait=1;
             last;
           }
         }else{
           $nowait=1;
         }
        }
      }#end of "check the state"-if
      $UPDATEPROCESS{"NumSLeft"}=$st->{"NumSearchesLeft"}[0];
      printUPDATEPROCESS($prdir,$curdir,$curspdyndir,$curspquendir,$st,$jobs,$quench,$dynmat,$turnoff_dynmat,$SearchesAlgo,$Screenoutput,$Screenoutput_tmp);
    }#end of the "for" loop.
    if(-e $Screenoutput_tmp) {system "mv $Screenoutput_tmp $Screenoutput";}
    # check if the current st have been finished and evaluate the found good saddle point
    # and submit more saddle point searches if necessary and possible
    print "===========================================================\n";
    if($ConvergenceAlgo == 1) {
      if($MaxSaddle eq "") { print "Found 0 counted unique saddles.\n";}
      else{$allowedmaxbarrier=$MaxSaddle - $st->{$curdir}[0];
       print "Found $belowTenKT unique saddles with barriers < ";
       printf("%3.3f",$allowedmaxbarrier);
       print " eV.\n";
      #print "Found $unique_sp unique saddles; $belowTenKT of them are below MinSaddle+$NumKT KT, i.e., $allowedmaxbarrier eV(barrier)\n";
      }
    }elsif($ConvergenceAlgo == 0) {
      print "Found total $total_sp saddle points; $unique_sp of them are unique saddle points.\n";
    }
    ($numprfin,$numprrun,$numrunjobs,$jobs,$st)=CheckSt($curdir,$quench,$dynmat,$numprfin,$numprrun,$numrunjobs,$MaxJobs,$NumSearches,$UseKDB,$SearchesAlgo,$ConvergenceAlgo,$Population,$NumKT,$Temperature,$DisplaceAlgo,$DisplaceRange,$NN_rcut,$MaxCoordNum,$nowait,$jobs,$st,$st_old);
    # information of a running state has been updated.
  }elsif($ststatus eq "done") {
    #print "$curdir is already done.\n";
  }else{
    die "Unrecognizable state status $ststatus.\n";
  } # end of the "Update state information" block. 
  #---------------------------------------------------------------------------
  # Decide what to do; quit or do one KMC simulation step
  #---------------------------------------------------------------------------
  if(lc($st->{"status"}[0]) eq "done"){
    print "State $curdir is done.  ";
    if($prest){print "Previous state is $prest.\n";}
    else{print "No previous state.\n";}
    print "Calculating rates.\n";
    if($repeat){# a repeated state; get the RateTableFile only
      print "State $curdir is a repeated state.\n";
      $rt=ReadEventTable($curdir,$Temperature,$RateTableFile);
      if(!$rt){
        print "No rate information found; it will be created from the state information.\n";
        $rt=WriteEventsTable($curdir,$DynMatFile,$turnoff_dynmat,$Temperature,$RateTableFile,$st,$Prefactor);
      }
    }else{# a new state
      # build event table and rate table and put barriers and prefactors into state hash.
      print "State $curdir is a new state.\n";
      if(-e "$curdir/$RateTableFile") { # This provides an opportunity for the user to disallow some phony processes found.
        print "Rate file found: $RateTableFile.  Reading information.\n";
        $rt=ReadEventTable($curdir,$Temperature,$RateTableFile);
        if(!$rt) {
          print "No rate information found; it will be created from the state information.\n";
          $rt=WriteEventsTable($curdir,$DynMatFile,$turnoff_dynmat,$Temperature,$RateTableFile,$st,$Prefactor);
        }
      }else{
        print "No rate file found, creating a new one.\n";
        $rt=WriteEventsTable($curdir,$DynMatFile,$turnoff_dynmat,$Temperature,$RateTableFile,$st,$Prefactor);
      }  
      # save a copy of the job file in the state directory
      system "cp $JobFile $curdir/";
    }
    #select an event
    $rdnum = rand();
    $rdnum2 = rand();
    ($selected,$elapse,$enforced,$time_saved) = SelectEvent($rt,$curdir,$prest,$rdnum,$rdnum2,$enforced,$time_saved,$Rdiffmax,$RateTableFile,$Temperature,$RateTableFileBolEqu,$akmc_step,$curstep,$BoltzmanEqu);
#-------------------------------------------------------------------------------------------------------------------------
    print "The chosen process for $curdir is $selected.\n";
    if($enforced ne ""){print "Process $enforced is enforced.\n";}
    # save the state and update the akmc_step after making a move:
    # setup akmc_step for next state, but the new state has to be set up later.
    ($curstep,$st,$akmc_step)=UpdateSt($curstep,$repeat,$selected,$rdnum,$elapse,$StFile,$StEnergyFile,$AkmcFile,$st,$akmc_step);
    # clear the job file (already copied into the previous state in update()
    %$jobs=();
    if($curstep>$AkmcSteps){
      $Done=1;
    }else{
      $Done=0;
    }
    print "===========================================================\n";
  }else{ # cur step is not finished, but this also means that the current state is not a repeated state
    #print "Record state information for $curdir to file $StFile.\n";
    WriteSt($curdir,$StFile,$st);
    #print "Summary of St: $numprrun running processes, $numprfin finished processes, and $numrunjobs running (and submitted) jobs.\n";
#    print "Total of $numprrun running processes and ",$numrunjobs-$numprrun," submitted jobs.\n";
    print "Total of $numprrun running processes and ",$numprfin," completed.\n";
#    print "$curdir is not done yet.\n";
    $Done=1;
  }
}# end of the "while" loop
WriteJobs($JobFile,$jobs); #always do this.
if($curstep > $AkmcSteps) {
  print "This aKMC simulation is finished.\n";
}else{
  print "This aKMC simulation is at step $curstep.\n";
}
print "===========================================================\n";


###############################################################################
# FUNCTIONS                                                                   #
###############################################################################
# ------------------------------------
# Read a paramter from the config file
# ------------------------------------
sub ReadParm{
  my ($configfile,$parm,$default)=@_;
  my $line="";
  my @line=();
  my @dummy=();
  my $value="";
  $value=$default;
  open (CONFIG,"<$configfile") or die "Can't open $configfile \n";
  while($line=<CONFIG>){
    $line=~s/^\s+//;
    if($line eq "" || $line eq "\n" || $line=~/^#/){ next;}
    @line=split(/=/,$line);
    $line[0]=~s/\s+$//;
    if(lc($line[0]) eq lc($parm)){
      @dummy=split(/\#/,$line[1]);
      $dummy[0]=~s/^\s+//;
      chomp($dummy[0]);
      $dummy[0]=~s/\s+$//; # get rid of possible trailing whitespace
      if($dummy[0] ne ""){ $value=$dummy[0];}
    }
  }
  close CONFIG;
  #print $parm."=".$value."\n";
  return $value;
}

# ------------------------------------
# Directory name
# ------------------------------------
sub DirName{
  my ($pre,$num)=@_;
  my $snum="";
  my $dirname="";
  $snum=sprintf "%04d",$num;
  $dirname=$pre.$snum;
  return $dirname;
}

# ------------------------------------
# CheckSystemConfig [($good,$stdirs,$numst,$curst)=CheckSystemConfig($maindir)]
# ------------------------------------
sub CheckSystemConfig{
  my $maindir=shift;
  my @stdirs=();
  my $numst=0;
  my $good=0;
  my $curst="";
  my @dist=();
  my ($homeposcar, $curst_poscar, $dist,@InMainDir);
  opendir MAINDIR, $maindir or die "can't open this $maindir dir!";
  @InMainDir=readdir MAINDIR;
  @stdirs = grep /^st\d{4}$/, @InMainDir;
  @stdirs=(sort @stdirs);
  closedir MAINDIR;
  $numst=@stdirs;
  if($numst>0){
    $curst=DirName("st", $numst);
    if($curst ne $stdirs[$numst-1]) {
      die "In the current folder, st folders are not ranked properly as st000x.\n";
    }else{
      $homeposcar="POSCAR";
      $curst_poscar=$curst."/POSCAR";
      if((-e $homeposcar) && (-e $curst_poscar)) {
        @dist=`diff $homeposcar $curst_poscar | sort`; # Do not mess up those two POSCARs. Keep them as they are.
        $dist=@dist;
        if(!$dist) {
          $good=1;
        }else{ # as long as "diff" ignores extra blank lines
          chomp($dist=`$Bin/dist.pl $homeposcar $curst_poscar | tail -1`);
          if($dist < 0.01) {$good =1;} # weaker since both files are supposed to be identical if they are not different
        }
      }
    }
  }else{
    $good=1;
  }
  return ($good,\@stdirs,$numst,$curst);
}

# ------------------------------------
# Make a new state directory
# ------------------------------------
sub MakeNewSt{
  my $curstdir=shift;
  if(-e "$curstdir"){
    die " Error in MakeNewSt: State directory $curstdir already exists;\n";
  }else{
    system "mkdir $curstdir"; 
  }
  system "cp POSCAR KPOINTS POTCAR $curstdir";
  system "cp INCAR_sp INCAR_min INCAR_dynmat $curstdir";
  system "cp DISPLACECAR_sp DISPLACECAR $curstdir"; 
  system "cp akmc_sp.sub akmc_min.sub akmc_dynmat.sub $curstdir";
}

# --------------------------------------
# Copy a st Hash
# --------------------------------------
sub StHashCopy{
  my $st=shift;
  my %stold=();
  my ($j,$i,$length);
  for $j (keys %$st) {
    $length=@{$st->{$j}};
    $stold{$j}=[()];
    for($i=0;$i<$length;$i++) {
      push @{$stold{$j}}, $st->{$j}[$i];
    } 
  }
  return \%stold;
}

# -----------------------------------------------
# Read akmc_step [%akmc_step=ReadAkmcStep($AkmcFile)]
# -----------------------------------------------
sub ReadAkmcStep{ # read the last two steps information
  my $AkmcFile=shift;
  my $line="";
  my @line=();
  my %akmc_step=();
  my $j=0;
  my $step=0;
  my @steps=();
  my $linesize=0;
  if(-e $AkmcFile){
    @line=`tail -1 $AkmcFile | sort`;
    for($step=0; $step< 1;$step++) {
      $line[$step]=~s/^\s+//;
      @steps=split(/\s+/,$line[$step]);
      if($steps[0]=~/^\d+$/) {
        $linesize++;
        $akmc_step{$steps[0]}=[()];
        for ($j=1;$j<@steps;$j++) {
          push @{$akmc_step{$steps[0]}}, lc($steps[$j]);
        }
      }
    }
    if(!$linesize) {die "Error In ReadAkmcStep, the last line of $AkmcFile is blank. Go fix it.\n";}
  }
  return \%akmc_step;
}

# -----------------------------------------------
# Read jobs [%jobs=ReadJobs($jobfilename)]
# -----------------------------------------------
sub ReadJobs{
  my $jobfilename=shift;
  my %jobs=();
  my $path;
  my $i=0;
  my @line=();
  my $line="";
  if(-e $jobfilename){
    open(JOB,"<$jobfilename") || die "In ReadJobs: cannot open file $jobfilename\n";
  }else{
    return \%jobs;
  }
  while($line=<JOB>){
    $line=~s/^\s+//;
    if($line eq "" || $line eq "\n"){ next; } # skip empty lines
    $i++;
    if($i!=1) {
      @line=split(/\s+/,$line);
      if($line[0] eq "") {next;}
      $jobs{$line[0]}=[()];
      for ($j=1;$j<@line;$j++) {
        push @{$jobs{$line[0]}}, lc($line[$j]);
      }
    }
  }
  close JOB; 
  return \%jobs;
}

# ---------------------------------------------------
# Read the paths storing states [ReadStPool($OldRunDirFile)]
# ---------------------------------------------------
sub ReadStPool{
  my $OldRunDirFile=shift;
  my @stpools=();
  my $line="";
  my @line=();
  if(-e $OldRunDirFile){
    open(JOB,"<$OldRunDirFile") || die "In ReadStPool: cannot open file $OldRunDirFile\n";
  }else{
    return \@stpools;
  }
  while($line=<JOB>){
    $line=~s/^\s+//;
    @line=split(/\s+/,$line);
    if($line[0] eq "") {next;}
    push @stpools, $line[0];
  }
  close JOB; 
  return \@stpools;
}

# ---------------------------------------------------
# Initialize aKMC: Initializer($min,$numrunjobs,$MaxJobs,$jobs,$akmc_step)
# ---------------------------------------------------
sub Initializer{
  my ($min,$numrunjobs,$MaxJobs,$jobs,$akmc_step)=@_;
  my $curdir=$min;
  my $InitDone=0;
  my $energy="na";
  my $force="na";
  my $dummy="";
  if(-e $curdir) {# there is already a min folder
    print "There is one initial min and no states.\nCheck the minimization status.\n";
    if(!exists($jobs->{$curdir})) {die "There is one min directory without any job information in the jobfile.\n";}
    #update jobs status -- check minimization,wrap up or resubmit the job.
    ($numrunjobs,$jobs)=HandleMins($curdir,$numrunjobs,$MaxJobs,$jobs);
    if(lc($jobs->{$curdir}[0]) eq "completed") {
      print "Initial min has been optimized.\nStart a new st for KMC.\n";
      $energy=$jobs->{$curdir}[3];
      $force=$jobs->{$curdir}[4];
      $InitDone=1;
    }elsif(lc($jobs->{$curdir}[0]) eq "killed"){
      print "Initial min has been killed due to fatal errors or by the user.\n We can't move on without initial min.\n";
    }
  }else{# need to create a new min
    MakeNewSt($curdir);
    ($numrunjobs,$jobs)=StartMins($curdir,$numrunjobs,$MaxJobs,$jobs);
  }
  $dummy="1";
  $akmc_step->{$dummy}=[()];
  push @{$akmc_step->{$dummy}}, "2bcreated","unknown","unknown",$energy,"na","0.0";
  return ($InitDone, $numrunjobs,$jobs, $akmc_step);
}

# ------------------------------------------------------------------------------------------------------------
# Read and write st energy file
# ------------------------------------------------------------------------------------------------------------
sub ReadStEnergyFile{
  my $curfile=shift;
  my %stenergy=();
  my @line=();
  my $line="";
  if( -e $curfile){
    open(JOB, "<$curfile") || die "In ReadStEnergyFile: cannot open $curfile\n";
  }else{
    return \%stenergy;
  }
  while($line=<JOB>){
    $line=~s/^\s+//;
    if($line eq "" || $line eq "\n"){next;} # ignore blank lines
    @line=split(/\s+/,$line);
    $stenergy{$line[0]}=$line[1]; 
    # it is the users reponsibilty to make sure these numbers are correct
  }
  close JOB;
  return \%stenergy;
}

sub AppendStEnergyFile{
  my $curdir=shift;
  my $StEnergyFile = shift;
  my $energy=shift;
  open(OUT, ">>$StEnergyFile");
  print OUT $curdir."  ".$energy."\n";
  close OUT;
}

# ------------------------------------------------------------------------------------------------------------
# Check if the current state is repeated in the state pool 
# ($repeat,$akmc_step)=CheckStpool($curstep,$stpool,$StEnergyFile,$Equivalency,$Ediffmax,$Rdiffmax,$akmc_step);
# ------------------------------------------------------------------------------------------------------------
sub CheckStpool{
  my ($curstep,$stpool,$stenergyfile,$Equivalency,$Ediffmax,$Rdiffmax,$akmc_step)=@_;
  my $repeat=0;
  my $i=0;
  my $j="";
  my ($stenergy,$filepath,$stpoolsize,$energy,$same,$POSCAR,$curst,@curst,$numst,@totalst,$linker);
  $Equivalency=1; # The equivalency checking is enforced here. This is worth the effort.
  $energy=$akmc_step->{$curstep}[3];
  $POSCAR="POSCAR";
  $stpoolsize=@$stpool;
  $stenergy=ReadStEnergyFile($stenergyfile);
  @totalst=(sort keys %$stenergy);
  $numst=@totalst;
  if($numst > 0){
    $curst=$totalst[$numst-1];
    #print $curst. "  ".$numst."\n";
    ($curst eq DirName("st",$numst)) || die "Error in CheckStpool: The state directories must be listed in order as st00xx\n";
  }
  for $j (@totalst){
    if (abs($energy-$stenergy->{$j}) < $Ediffmax) {
      $POSCAR_j=$j."/POSCAR";
      $same=CompareTwoPOSCAR($POSCAR,$POSCAR_j,$Rdiffmax);
      if($Equivalency && !$same) {
        ($same,$linker)=indistinguishable($POSCAR,$POSCAR_j,$Rdiffmax);
      }
      if($same) {
        $repeat=1;
        $akmc_step->{$curstep}[0]=$j;
        $akmc_step->{$curstep}[1]="repeat";
        $akmc_step->{$curstep}[2]=$j;
        last;
      }
    }
  }
  if($repeat){ return ($repeat, $akmc_step);}
  for ($i=0;$i<$stpoolsize;$i++){
    $filepath=$stpool->[$i]."/".$stenergyfile;
    $stenergy=ReadStEnergyFile($filepath);
    for $j (keys %$stenergy){
      if (abs($energy-$stenergy->{$j}) < $Ediffmax) {
        $POSCAR_j=$stpool->[$i]."/".$j."/POSCAR";
        $same=CompareTwoPOSCAR($POSCAR,$POSCAR_j,$Rdiffmax);
        if($Equivalency && !$same) {
          ($same,$linker)=indistinguishable($POSCAR,$POSCAR_j,$Rdiffmax);
        } 
        if($same) {
          $repeat=1;
          $akmc_step->{$curstep}[0]=$j;
          $akmc_step->{$curstep}[1]="repeat";
          $akmc_step->{$curstep}[2]=$stpool->[$i]."/".$j;
         last;
        }
      }
    }
  }
  if(!$repeat){
    $akmc_step->{$curstep}[0]=DirName("st", $numst+1);
    $akmc_step->{$curstep}[1]="new";
    $akmc_step->{$curstep}[2]="na";
  }
  return ($repeat, $akmc_step);
}

# ----------------------------------------------------------------------
# Submit jobs [($jobid, $jobstatus)=SubmitJobs($curdir,$jobtype)]
# We assume the existence of regular files for vasp; the submit script will check them.
# ----------------------------------------------------------------------
sub SubmitJobs{
  my ($curdir,$jobtype)=@_;
  my $jobid="";
  my $jobstatus="";
  my $jobinfo="";
  my $difference="";
  my @POSCAR_size=();
  my $POSCAR_size=0;
  my ($presentworkdir,$akmchome);
  my $filemissing="";
  if($jobtype eq "minimization") {
     if(-e "$curdir/INCAR_min" && -e "$curdir/INCAR"){
        $difference=`diff $curdir/INCAR_min $curdir/INCAR | tail -1`;
        if ($difference ne ""){
	  system "cp $curdir/INCAR_min $curdir/INCAR";
        }
     }elsif(-e "$curdir/INCAR_min"){
        system "cp $curdir/INCAR_min $curdir/INCAR";
     }else{
        #die "Error: in $curdir: no INCAR_min\n";
        $filemissing="INCAR_min";
     }
     if(-e "$curdir/akmc_min.sub"){
        system "cp $curdir/akmc_min.sub $curdir/akmc.sub";
     }else{
        #die "Error: in $curdir: no akmc_min.sub\n";
        $filemissing=$filemissing."akmc_min.sub";
     }
     @POSCAR_size=`cat $curdir/POSCAR | sort`;
     $POSCAR_size=@POSCAR_size;
  }elsif($jobtype eq "dimer"){
     if(-e "$curdir/INCAR_sp" && -e "$curdir/INCAR"){
        $difference=`diff $curdir/INCAR_sp $curdir/INCAR | tail -1`;
        if ($difference ne ""){
          system "cp $curdir/INCAR_sp $curdir/INCAR";
        }
     }elsif(-e "$curdir/INCAR_sp"){
        system "cp $curdir/INCAR_sp $curdir/INCAR";
     }else{
        #die "Error: in $curdir: no INCAR_sp\n";
        $filemissing="INCAR_sp";
     }
     if(-e "$curdir/akmc_sp.sub"){
        system "cp $curdir/akmc_sp.sub $curdir/akmc.sub";
     }else{
        #die "Error: in $curdir: no akmc_sp.sub\n";
        $filemissing=$filemissing."akmc_sp.sub";
     }
     @POSCAR_size=`cat $curdir/POSCAR | sort`;
     if(@POSCAR_size <= 6) {
       $POSCAR_size=@POSCAR_size;
     }else{
       @POSCAR_size=`cat $curdir/POSCAR | sort`;
       $POSCAR_size=@POSCAR_size;
     }
     if(!(-e "$curdir/MODECAR")) {$filemissing=$filemissing."MODECAR";} # we pay attention to special files
  }elsif($jobtype eq "lanczos"){
     if(-e "$curdir/INCAR_sp" && -e "$curdir/INCAR"){
        $difference=`diff $curdir/INCAR_sp $curdir/INCAR | tail -1`;
        if ($difference ne ""){
          system "cp $curdir/INCAR_sp $curdir/INCAR";
        }
     }elsif(-e "$curdir/INCAR_sp"){
        system "cp $curdir/INCAR_sp $curdir/INCAR";
     }else{
        #die "Error: in $curdir: no INCAR_sp\n";
        $filemissing="INCAR_sp";
     }
     if(-e "$curdir/akmc_sp.sub"){
        system "cp $curdir/akmc_sp.sub $curdir/akmc.sub";
     }else{
        #print "Error: in $curdir: no akmc_sp.sub\n";
        $filemissing=$filemissing."akmc_sp.sub";
     }
     @POSCAR_size=`cat $curdir/POSCAR | sort`;
     $POSCAR_size=@POSCAR_size;
     if(!(-e "$curdir/MODECAR")) {$filemissing=$filemissing."MODECAR";}
  }elsif($jobtype eq "dynmat"){
     if(-e "$curdir/INCAR_dynmat" && -e "$curdir/INCAR"){
        $difference=`diff $curdir/INCAR_dynmat $curdir/INCAR | tail -1`;
        if ($difference ne ""){
          system "cp $curdir/INCAR_dynmat $curdir/INCAR";
        }
     }elsif(-e "$curdir/INCAR_dynmat"){
        system "cp $curdir/INCAR_dynmat $curdir/INCAR";
     }else{
        #print "Error: in $curdir: no INCAR_dynmat\n";
        $filemissing="INCAR_dynmat";
     }
     if(-e "$curdir/akmc_dynmat.sub"){
        system "cp $curdir/akmc_dynmat.sub $curdir/akmc.sub";
     }else{
        #print "Error: in $curdir: no akmc_sp.sub\n";
        $filemissing=$filemissing."akmc_dynmat.sub";
     }
     @POSCAR_size=`cat $curdir/POSCAR | sort`;
     $POSCAR_size=@POSCAR_size;
     if(!(-e "$curdir/DISPLACECAR")) {$filemissing="DISPLACECAR";}
  }else{
     die "In SubmitJobs: unknown job type.\n";
  }
  if(!$filemissing && $POSCAR_size > 6){
    chomp($akmchome=`pwd`);
    chomp($presentworkdir=`cd $curdir; pwd`);
    chomp($jobinfo=`cd $curdir; $Bin/akmc_submit.pl $presentworkdir $akmchome | tail -1;`);
    $jobinfo=~s/^\s+//;
    @jobinfo=split(/\s+/,$jobinfo);
    $jobid=$jobinfo[0];
    if($jobid=~/^\d+$/) { # for the time being, we require a digital id. This can be changed in the future.
      $jobstatus="submitted";
    }elsif(lc($jobid) eq "2bsubmitted") {# This comes from an unsuccessful job submission due to all kinds of reasons
      $jobstatus="2bsubmitted";
      $jobid="FailedSubmission"; # maybe the MaxJobs is larger than the number of allowed jobs in the system queue?
      print "Warning from SubmitJobs: in $curdir, just tried to submit a job in $curdir, but it failed. Better check it(marked as $jobstatus)\n";
    }else{# in case that the user passes back something else. To avoid endless repeated trials, we mark it as killed.
      $jobstatus="killed";
      $jobid="UnrecognizableJobid";
    }
  }elsif($filemissing){
    #print "Warning from SubmitJobs: in $curdir, $filemissing is missing. We need that file to run jobs. It is set as killed. Pls submit it manually if it is still needed.\n";
    $jobstatus="killed";
    $jobid="Missing$filemissing";
  }else{
    #print "Warning from SubmitJobs: in $curdir, POSCAR is no good. Probably a failed job that was cleaned up. It is set as killed. Pls submit it manually if it is still needed.\n";
    $jobstatus="killed";
    $jobid="BadPOSCAR";
  }
  return ($jobid, $jobstatus);
}

# ------------------------------------------
# Check jobs [%jobs=CheckJobs(%jobs)]
# ------------------------------------------
sub CheckJobs{
  my $jobs=shift;
  my $jobstatus="";
  my $jobid="";
  my $jobinfo="";
  my $dir="";
  my $dummy="";
  my $numrunjobs=0;
  for $dir (keys %$jobs){
    if($dir eq "") {
       delete $jobs->{$dir};
       next;
    }
    $jobstatus=lc($jobs->{$dir}[0]);
    if($jobstatus eq "completed" || $jobstatus eq "stopped" || $jobstatus eq "2bsubmitted"){
      next;
    }elsif($jobstatus eq "running" || $jobstatus eq "queue"  || $jobstatus eq "submitted"){
      $jobid=$jobs->{$dir}[1];
      chomp($jobinfo=`$Bin/akmc_check.pl $jobid | tail -1;`);
      #print "***$jobinfo***\n";
      $jobinfo=~s/^\s+//;
      @jobinfo=split(/\s+/,$jobinfo);
      $dummy=lc($jobinfo[0]);
      #if($dummy eq "running" || $dummy eq "queue" || $dummy eq "stopped"){
      if($dummy eq "running" || $dummy eq "queue") {
        $jobs->{$dir}[0]=$dummy;
        $numrunjobs++;
      }elsif($dummy eq "stopped") {
        $jobs->{$dir}[0]=$dummy;
      }else{
        die "In CheckJobs: Unrecognizable job status from $Bin/akmc_check.pl.\n";
      }
    }elsif($jobstatus eq "killed"){
      #print "Warning from CheckJobs: the job in $dir has been killed. No action.\n";
    }else{
        die "In CheckJobs: Unrecognizable job status from the job hash.\n";
    }
  } 
  return ($numrunjobs,$jobs);
}

# -----------------------------------------------------------------
# Kill a job : $succeed=KillJobs($curspdir, $jobs->{$curspdir}[1]);
# -----------------------------------------------------------------
sub KillJobs{
  my ($curspdir,$jobid)=@_;
  my $killstatus=0;
  my $line="";
  my @line=();
  my $akmchome="";
  chomp($akmchome=`pwd`);
  chomp($line=`cd $curspdir; $Bin/akmc_kill.pl $jobid $akmchome | tail -1;`);
  $line=~s/^\s+//;
  @line=split(/\s+/,lc($line));
  if($line[0] eq "killed"){
    $killstatus=1;
  }
  return $killstatus;
}

# --------------------------------------------------------------------------------
# Submit a job : ($numrunjobs,%jobs)=HandleJobSubmission($numrunjobs, $MaxJobs, $jobdir, $jobtype, %jobs)
# A robust way to wrap the job submission under the context of maxJobs limitation
# --------------------------------------------------------------------------------
sub HandleJobSubmission{
  my ($numrunjobs,$MaxJobs,$jobdir,$jobtype,$jobs)=@_;
  my $jobid="na";
  my $jobstatus="2bsubmitted";
  if(!(exists($jobs->{$jobdir}))) {# This is the first time for the job to be submitted in jobdir
    $jobs->{$jobdir}=[()];
    if($jobtype ne "minimization" || $jobtype ne "dynmat"){ # it is saddle searching job
      push @{$jobs->{$jobdir}}, $jobstatus, $jobid, $jobtype, "na", "na", "na";
    }else{
      push @{$jobs->{$jobdir}}, $jobstatus, $jobid, $jobtype, "na", "na";
    }
  } # we could use a flag called "firsttime" to replace the hash element judgment
  if($numrunjobs < $MaxJobs) {
    ($jobid,$jobstatus)=SubmitJobs($jobdir,$jobtype);
    $jobs->{$jobdir}[0]=$jobstatus;
    $jobs->{$jobdir}[1]=$jobid;
    if($jobstatus eq "submitted") { # a job was successfully submitted
      $numrunjobs++;
    }
  }else{
    $jobs->{$jobdir}[0]=$jobstatus;
    $jobs->{$jobdir}[1]=$jobid;
  }
  return ($numrunjobs,$jobs);
}
# -----------------------------------------------------------
# Get the POSCAR for the next state from curdir
# -----------------------------------------------------------
sub GetStPOSCAR{
  my $curdir=shift;
  my $StPOSCAR="";
  if(-e "$curdir/final") {# a cleaned-up folder
     $StPOSCAR=$curdir."/POSCAR";
  }elsif(-e "curdir/CONTCAR"){
     $StPOSCAR=$curdir."/CONTCAR";
  }
  system "cp $StPOSCAR POSCAR";
  return
}
# -----------------------------------------------------------
#PrepareUPDATEPROCESS() and printUPDATEPROCESS() and writeUPDATEPROCESS()
# ----------------------------------------------------------- 
sub PrepareUPDATEPROCESS{
  %UPDATEPROCESS=(
    "was" => "", "now" => "---  ",
    "energy" => " ---  ","force" => "  ---  ","curvature" => "  ---  ",
    "jobstatus" => "---  ","NumSLeft" => "  ---  ","comment" => ""
  );
}

sub printUPDATEPROCESS_header{
   my ($Screenoutput,$Screenoutput_tmp)=@_;
    if($ARGV[0] eq "amy"){ 
      $Screenoutput_header="process\twas\tenergy\tforce\tcurvat.\tnow\tjob\t\tNSL\tcomments\n";
    }else{
#      $Screenoutput_header="process\tenergy\tforce\tcurvat.\tstatus\tcomments\n";
      $Screenoutput_header ="process    energy     force   curvature  status   comment\n";
      $Screenoutput_header.="--------  --------  --------  --------  --------  --------\n";
    }
   print $Screenoutput_header;
   if($Screenoutput){
     open(SCREENOUTHEADER, ">$Screenoutput_tmp") || die "cannot open $Screenoutput_tmp file.\n";
     print SCREENOUTHEADER $Screenoutput_header;
     close(SCREENOUTHEADER);
   }
}

sub printUPDATEPROCESS_simple{
  my ($prdir,$stenergy,$Screenoutput,$Screenoutput_tmp)=@_;
  my $line="";
  my $barrier="";
  $line=$prdir."\t".$UPDATEPROCESS{"was"}."\t";
  if($UPDATEPROCESS{"energy"}=~/^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$/){
    $barrier=$UPDATEPROCESS{"energy"}-$stenergy;
    $line=$line.sprintf("%3.3f",$barrier)."\t";
  }else{$line=$line.$UPDATEPROCESS{"energy"}."\t";}
  if($UPDATEPROCESS{"force"}=~/^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$/){
    $line=$line.sprintf("%3.5f",$UPDATEPROCESS{"force"})."\t";
  }else{$line=$line.$UPDATEPROCESS{"force"}."\t";}
  if($UPDATEPROCESS{"curvature"}=~/^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$/){
    $line=$line.sprintf("%3.3f",$UPDATEPROCESS{"curvature"})."\t";
  }else{$line=$line.$UPDATEPROCESS{"curvature"}."\t";}
  if($UPDATEPROCESS{"now"}){
    $line=$line.$UPDATEPROCESS{"now"}."\t".$UPDATEPROCESS{"jobstatus"}."\t\t";
  }else{
    $line=$line."---   \t".$UPDATEPROCESS{"jobstatus"}."\t\t";
  }
  $line=$line.$UPDATEPROCESS{"NumSLeft"}."\t";
  $line=$line.$UPDATEPROCESS{"comment"}."\n";
  print $line;
  if($Screenoutput ne ""){
    open(OUTSCREENFILE, ">>$Screenoutput_tmp") || die "in printUPDATEPROCESS_simple, cannot open file $Screenoutput_tmp\n";
    print OUTSCREENFILE $line;
    close(OUTSCREENFILE);
  }
}
sub printUPDATEPROCESS{
  my ($prdir,$curdir,$curspdyndir,$curspquendir,$st,$jobs,$quench,$dynmat,$turnoff_dynmat,$SearchAlgo,$Screenoutput,$Screenoutput_tmp)=@_;
  my ($curspdir,$curminj,$energy,$force,$curvature,$barrier,$j,$minj_short,$now,$was,$line);
  if($ARGV[0] eq "amy"){
    printUPDATEPROCESS_simple($prdir,$st->{$curdir}[0],$Screenoutput,$Screenoutput_tmp); 
    return;
  }
  $curspdir=$curdir."/".$prdir;
  $line=" ".$prdir."  ";
  if(exists($jobs->{$curspdir})){
    $energy=$jobs->{$curspdir}[3];
    $force=$jobs->{$curspdir}[4];
    $curvature=$jobs->{$curspdir}[5];
  }else{
    ($energy,$force,$curvature)=GetSPEF($curspdir,$SearchAlgo);
  }
  if($energy=~/^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$/){
    $barrier=$energy - $st->{$curdir}[0];
    $line.=sprintf("%9.3f",$barrier)." ";
  }else{$line.="    --    ";}
  if($force=~/^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$/){
    $line.=sprintf("%9.3f",$force)." ";
  }else{$line.="    --    ";}
  if($curvature=~/^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$/){
    $line.=sprintf("%9.3f",$curvature)." ";
  }else{$line.="    --    ";}
  
  $line.=sprintf("%8s",$UPDATEPROCESS{"was"})."  ";
  $was=lc($UPDATEPROCESS{"was"});
  $now=lc($st->{$prdir}[0]);
  if($UPDATEPROCESS{"now"} ne ""){
    if($was ne $now){$line.=$UPDATEPROCESS{"now"}." ";}
  }
  if($UPDATEPROCESS{"comment"} ne ""){
    $line.=$UPDATEPROCESS{"comment"}."\n";
  }else{
    if($now eq "search"){ $line.=$jobs->{$curspdir}[0];
    }elsif($now eq "dynmat"){$line.=$jobs->{$curspdyndir}[0];}
    $line.="\n";
  }

  if((-e "$curspquendir/min1")&&(-e "$curspquendir/min2")){# there are quench mins
     for($j=1;$j<3;$j++){
       $minj_short=$quench."/min".$j;
       $curminj=$curspdir."/".$minj_short;
       $line.="  min".$j."   ";
       if(exists($jobs->{$curminj})){
         $energy=$jobs->{$curminj}[3];
         $force=$jobs->{$curminj}[4];
       }else{
         ($energy,$force)=GetEF($curminj);
       }
       if($energy=~/^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$/){
         $barrier=$energy - $st->{$curdir}[0];
         $line.=sprintf("%9.3f",$barrier)." ";
       }else{$line.="    --    ";}
       if($force=~/^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$/){
          $line.=sprintf("%9.3f",$force)." ";
       }else{$line.="    --    ";}
       $line.="    --    ";
       if(lc($st->{$prdir}[0]) eq "quench"){ #quenching
         $line.="     run  ";
         $line.=$jobs->{$curminj}[0];
       }elsif(lc($st->{$prdir}[0]) eq "done" || lc($st->{$prdir}[0]) eq "dynmat"){ # quenched
         $line.="    done  ";
         if(lc($st->{$prdir}[1]) eq "good"){
           if($minj_short eq $st->{$prdir}[3]){ # minj is the final state
             $line.="final state";
           }else{
             $line.="initial state";
           }
         }elsif(lc($st->{$prdir}[1]) eq "bad" || lc($st->{$prdir}[1]) eq "highenergy"){
           $line.="";
         }
       }elsif(lc($st->{$prdir}[0]) eq "killed"){
         if(exists($jobs->{$curminj}) && lc($jobs->{$curminj}[0]) eq "killed"){
           $line.=" killed   ";
           if($jobs->{$curminj}[1]=~/^\d+$/){ # must use digits as job id.
             $line.="by user";
           }else{
             $line.="$jobs->{$curminj}[1]";
           }
         }else{
           $line.=$jobs->{$curminj}[0]." ";
         }
       }else{
         die "In printUPDATEPROCESS: after a search status, it has to be either quench, dynmat, done or killed\n";
       }
       $line.="\n";
     }#end of the for loop
     if(!$turnoff_dynmat && (-e "$curspdyndir")){# doing dynmat , but we need to make sure
         $line.=$dynmat."   ";
         $line.="    --        --        --    ";
         if(lc($st->{$prdir}[0]) eq "done"){
           $line.=" done   ";
         }elsif(lc($st->{$prdir}[0]) eq "dynmat"){
           $line.=" dynmat   ";
           $line=$line.$jobs->{$curspdyndir}[0];
         }elsif(lc($st->{$prdir}[0]) eq "killed"){
           $line.=" killed   ";
           if(exists($jobs->{$curspdyndir}) && lc($jobs->{$curspdyndir}[0]) eq "killed"){
             if($jobs->{$curminj}[1]=~/^\d+$/){
               $line.="by user";
             }else{
               $line=$line."$jobs->{$curminj}[1]";
             }
           }else{
             print $st->{$prdir}[1];
           }
         }elsif(lc($st->{$prdir}[0]) eq "quench" || lc($st->{$prdir}[0]) eq "search"){ # well, it is weird
           $line.=" ahead   ";
         }else{
           die "In printUPDATEPROCESS: after a search status, a dynmat has to be either dynmat, done, or killed\n";
         }
         $line.="\n";
     }
  }#end of "min1 min2 " if
  print $line;
  #print "------------------------------------------------------\n";
  if($Screenoutput){
    open(OUTSCREENFILE, ">>$Screenoutput_tmp") || die "In printUPDATEPROCESS, cannot open file $Screenoutput_tmp\n";
    print OUTSCREENFILE $line;
    #print OUTSCREENFILE "------------------------------------------------------\n";
    close(OUTSCREENFILE);
  } 
}
# -----------------------------------------------------------
# $jobs=PatchupRecycledJobInfo($jobs,$curspdir,$curspdyndir,$curspquendir,$quench,$dynmat,$turnoff_dynmat,$SearchesAlgo);
# -----------------------------------------------------------
sub PatchupRecycledJobInfo{
  my ($jobs,$curspdir,$curspdyndir,$curspquendir,$quench,$dynmat,$turnoff_dynmat,$SearchesAlgo)=@_;
  my ($energy,$force,$curvature,$j,$minj);
  #get saddle search information
  ($energy,$force,$curvature)=GetSPEF($curspdir,$SearchesAlgo);
  if($energy eq "na" || $force eq "na" || $curvature eq "na"){
    print "Serious warning: a recovered process($curspdir) has no energy/force/curvature. Check it\n";
  }
  $jobs->{$curspdir}=[("completed","jobid_marker",$SearchesAlgo,$energy,$force,$curvature)]; 
  #get quench information
  for($j=1;$j<3;$j++){
    $minj=$curspquendir."/min".$j;
    ($energy,$force)=GetEF($minj);
    if($energy eq "na" || $force eq "na"){
      print "Serious warning: a recovered process($curspdir) has no energy/force for quench. Check it\n";
    }
    if(!(exists($jobs->{$minj}))){
      $jobs->{$minj}=[("completed","jobid_marker".$j,"minimization",$energy,$force)];
    }
  }
  #get dynmat information
  if(-e "$curspdyndir/final"){# it is assumed done here
    if(!(exists($jobs->{$curspdyndir}))){
      $jobs->{$curspdyndir}=[("completed","jobid_marker","dynmat","na","na")];
    }
  }
  return $jobs;
}
# -----------------------------------------------------------
# Check if a stopped vasp process is converged: $converged=CheckMins($prdir,$energy,$force)]
# -----------------------------------------------------------
sub CheckMins{
  my ($prdir,$energy,$force)=@_;
  my $check="";
  my $line="";
  my @line=();
  my $dummy="";
  my @dummy=();
  my $fmax="";
  my $outzipped=0;
  my ($converged,$ediff);
  my $zip = $ENV{'VTST_ZIP'};
  if($zip eq ''){ $zip = 'gzip' ; }
  # "final" folder rules!
  # GetEF/GetSPEF checks $prdir/final/OUTCAR* first,then check $prdir/OUTCAR* and never checks $pddir/inter
  $converged=0;
  if((-e "$prdir/final/OUTCAR.gz") || (-e "$prdir/final/OUTCAR.bz2") || (-e "$prdir/final/OUTCAR")){
    #print "NO outcar(.gz,.bz2) in $prdir; jobs.dat said this job was still running, but we just found final/OUTCAR* in the directory. We ASSUME it is converged.\n";
    if($energy eq "na" || $force eq "na"){ die "In CheckMins, in $prdir, no energy/force with the existence of final/OUTCAR*.\n";}
    $converged=1;
  }elsif((-e "$prdir/OUTCAR") || (-e "$prdir/OUTCAR.gz") || (-e "$prdir/OUTCAR.bz2")) {
    if($energy eq "na" || $force eq "na"){ return $converged; }
    if(!(-e "$prdir/OUTCAR")){
      if(-e "$prdir/OUTCAR.gz"){ $outzipped=1; system "cd $prdir; gunzip OUTCAR.gz";}
      if(-e "$prdir/OUTCAR.bz2"){ $outzipped=1; system "cd $prdir; bunzip2 OUTCAR.bz2";}
    }
    chomp($check=`cd $prdir; grep "reached required accuracy" OUTCAR | tail -1`);
    if($check){
      $converged=1;
    }elsif(-e "$prdir/INCAR"){
      $ediff=0.0001; # vasp default value
      open (incar,"<$prdir/INCAR") || die "Can't open INCAR file in $prdir\n";      
      while($line=<incar>){
        $line=~s/^\s+//;
        if($line eq "" || $line eq "\n" || $line=~/^#/){next;}
        @line=split(/=/,$line);
        $line[0]=~s/\s+$//;
        if(!(lc($line[0])=~/^ediff/)){ next; }
        @dummy=split(/\#/,$line[1]);
        $dummy2=$dummy[0];
        $dummy2=~s/^\s+//;
        @dummy2=split(/\s+/,$dummy2); 
        $dummy2[0]=~s/^\s+//;
        chomp($dummy2[0]);
        $dummy2[0]=~s/\s+$//; # get rid of possible trailing whitespace
        if(lc($line[0]) eq "ediffg" && $dummy2[0]=~/^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$/ && $dummy2[0]< 0)  {
          $fmax=abs($dummy2[0]);
        }elsif(lc($line[0]) eq "ediff" && $dummy2[0]=~/^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$/){
          $ediff=$dummy2[0];
        }
      }
      close incar;
      if($fmax eq ""){
        print "Warning In CheckMins: There is no EDIFFG(negative) in the INCAR in the directory $prdir. Use vasp default.\n";
        $fmax=10*$ediff;
      }
      if($force < $fmax) {$converged=1;}
    }
    if($outzipped) { system "cd $prdir; $zip OUTCAR";}
  }elsif((-e "$prdir/inter/OUTCAR") || (-e "$prdir/inter/OUTCAR.gz") || (-e "$prdir/inter/OUTCAR.bz2")){
    #print "No OUTCAR (.gz,.bz2) in $prdir; jobs.dat said this job was still running, but we just found inter (no final) in the directory. We ASSUME it is not converged.\n";
    $converged=0;
  }else{
    $converged=0;
    #print "In CheckMins: There is no OUTCAR(.gz, .bz2) file in the directory $prdir(/final or /inter). Job might never started", "\n";
  }
  return $converged;
}

# -------------------------------------------------------------
# Check if a stopped DynMat job is completed [$converged=CheckDynmat($curpredir)]
# -------------------------------------------------------------
sub CheckDynmat{
  my $curpredir=shift;
  my $converged=0;
  my $line="";
  my @line=();
  my $dummy="";
  my @dummy=();
  my $steps="";
  my @loops=();
  my $loops="";
  my $outzipped=0;
  my $zip = $ENV{'VTST_ZIP'} ;
  if($zip eq ''){ $zip = 'gzip' ; }
  if((-e "$curpredir/final/OUTCAR") || (-e "$curpredir/final/OUTCAR.gz") || (-e "$curpredir/final/OUTCAR.bz2")){
    #print "There is a final folder in the dynmat calculations. We assume this dynmat has converged.\n";
    $converged=1;
  }elsif((-e "$curpredir/OUTCAR") || (-e "$curpredir/OUTCAR.gz") || (-e "$curpredir/OUTCAR.bz2")) {
    if(!(-e "$curpredir/OUTCAR")){
      if(-e "$curpredir/OUTCAR.gz") {$outzipped=1; system "cd $curpredir; gunzip OUTCAR.gz";}
      if(-e "$curpredir/OUTCAR.bz2"){$outzipped=1; system "cd $curpredir; bunzip2 OUTCAR.bz2";}
    }
    $loops=`cd $curpredir; grep "FORCES: max atom, RMS" OUTCAR;`;
    @loops=split(/\n/, $loops);
    open(incar,"<$curpredir/INCAR") or die "Can't open INCAR file in $curpredir.\n";
    while($line=<incar>){
       $line=~s/^\s+//;
       if($line eq "" || $line eq "\n" || $line=~/^#/){next;}
       @line=split(/=/,$line);
       $line[0]=~s/\s+$//;
       @dummy=split(/\#/,$line[1]);
       $dummy[0]=~s/^\s+//;
       chomp($dummy[0]);
       $dummy[0]=~s/\s+$//;
       if(lc($line[0]) eq "nsw" && $dummy[0]=~/^\d+$/ && $dummy[0] > 0)  {$steps=$dummy[0];}
    }
    close incar;
    if($steps eq "") {die "In CheckDynmat: There is no NSW in the INCAR in the directory $curpredir.\n";}
    if ($steps == @loops) {
       $converged=1;
    }else{
       $converged=0;
    }
    if($outzipped) {system "cd $curpredir; $zip OUTCAR";}
  }else{
    $converged=0;
    #print "In CheckDynmat: There is no OUTCAR(.gz,.bz2) file in the directory $curpredir. Job might never start.\n";
  }
  return $converged;
}

# ---------------------------------------------------------------------------------
# Check if a stopped job has finished at least one ionic cycle.  FinishOneIonicCycle()
# ---------------------------------------------------------------------------------
sub FinishOneIonicCycle{
  my ($curdir,$jobtype)=@_;
  my $yes=0;
  my ($energy,$force,$curvature);
  if($jobtype eq "dimer"){
    #($energy,$force,$curvature)=GetSPEF($curdir,"dimer");
    ($energy,$force)=GetEF($curdir);
    if($energy ne "na" && $force ne "na"){
       #if(-e "$curdir/CENTCAR") {$yes=1;}
       $yes=1;
    }
  }elsif($jobtype eq "minimization"){
    ($energy,$force)=GetEF($curdir);
    if($energy ne "na" && $force ne "na") {$yes=1;}
  }elsif($jobtype eq "dynmat"){
    $yes=0;
  }elsif($jobtype eq "lanczos"){
    die "In FinishOneIonicCycle lanczos is not ready yet\n";
  }else{
    die "In FinishOneIonicCycle unrecoganizable $jobtype\n";
  }
  return $yes;
}

# -----------------------------------------------------
# Cleanup jobs [CleanupJobs($curdir,$converged,$jobtype)]
# -----------------------------------------------------
sub CleanupJobs{
  my ($curdir,$good,$jobtype)=@_;
  my ($foldername, $fname);
  my $update_poscar=1;
  if($good){
    $foldername="final";
    if(-e "$curdir/final"){ 
      #print "Warning:in $curdir folder final already exists and we will do nothing.\n";
      return
    }
  }else{
    #print "it is me $curdir\n";
    $foldername="inter";
    if(-e "$curdir/inter"){
      $fname=GetInter($curdir);
      #print "Warning:in $curdir folder inter already exists and will be renamed as $fname.\n";
      system "cd $curdir; mv inter $fname";
    }
    # has job finished one ionic cycle so that CONTCAR is not blank?
    $update_poscar=FinishOneIonicCycle($curdir,$jobtype);
    if(!$update_poscar) {system "cd $curdir; cp POSCAR CONTCAR";}
  }
  # make sure ll_out etc. is removed
  system("cd $curdir; $Bin/vfin.pl $foldername > /dev/null");
  if(!$update_poscar){system "cd $curdir/$foldername;rm -f CONTCAR CONTCAR.con CONTCAR.xyz";} # LX: shall we delete this folder? 
  return;
}

sub GetInter {
  my $curdir=shift;
  my $fname="";
  my @inters=(); 
  opendir InDir, $curdir or die "cannot open $curdir\n";
  @inters=readdir InDir;
  $fname = grep /^inter\d/, @inters;
  $fname="inter".$fname;
  closedir InDir;
  return $fname;
}
# ---------------------------------------------------------------------------------------
# Calculate normal mode frequencies [($twonegatives,$frequencyfilepath)=DynMatrix($curpredir,$DynMatFile,$max)]
# ---------------------------------------------------------------------------------------
sub DynMatrix {
  my ($curpredir,$DynMatFile,$maximumnegative)=@_;
  my ($frequencyfilepath,$curpredirdyn);
  my $imaginery=0;
  my $twonegatives=0;
  my @line=();
  my $outzipped=0;
  my $zip = $ENV{'VTST_ZIP'} ;
  my $dummy="";
  if($zip eq ''){ $zip = 'gzip' ;}
  $curpredirdyn=$curpredir;
  if(-e "$curpredir/final"){ # a done cleaned up dynmat
    $curpredirdyn=$curpredir."/final";
  }else{
    $curpredirdyn=$curpredir;
  }
  $frequencyfilepath=$curpredirdyn."/".$DynMatFile;
  if(-e "$frequencyfilepath") {system "rm $frequencyfilepath";}
  if(-e "$curpredirdyn/OUTCAR.gz") {$outzipped=1; system "cd $curpredirdyn; gunzip OUTCAR.gz";}
  elsif(-e "$curpredirdyn/OUTCAR.bz2") {$outzipped=1; system "cd $curpredirdyn; bunzip2 OUTCAR.bz2";}
  $dummy=`cd $curpredirdyn; $Bin/dymmatrix.pl DISPLACECAR OUTCAR`;
  #print "dynmat dir: $curpredir"."\n";
  if(-e "$frequencyfilepath"){
    open(DYNMAT, "<$frequencyfilepath") || die "$frequencyfilepath cannot be opened.\n";
    while($line=<DYNMAT>) {
      $line=~s/^\s+//;
      @line=split(/\s+/,$line);
      if($line[3] == 1){$imaginery++;}
      if($imaginery>$maximumnegative){$twonegatives=1;}
    }
    close DYNMAT;
  }else{
    die "$frequencyfilepath was not successfully created by dynmatrix.pl\n";
  }
  if($outzipped) {system "cd $curpredir; $zip OUTCAR";}
  return ($twonegatives,$frequencyfilepath);
}

# --------------------------------------------------------------------------------------------------
# What if DynMatrix gives a bad dynmat? DealwithBadDynMatrix($curpredir)
# --------------------------------------------------------------------------------------------------
sub DealwithBadDynMatrix{
  my $curpredir=shift;
  my $nothing;
  $nothing="";
  #print "A bad dynmat was generated. We could either increase or decrease the finite displacement\n";
  #print "or we may to reconverge the saddle to a better resolution, etc.,\n";
  #print "but right now, we don't do anything.\n";
  # a standard prefactor 1E12 will be used if this process is still counted.
}

# --------------------------------------------------------------------------------------------------
# Get force and energy from an OUTCAR(.gz,.bzip2) before/after cleanup [($energy,$force)=GetEF($dir)]
# --------------------------------------------------------------------------------------------------
sub GetEF{
  my $dir=shift;
  my ($energy,$force);
  my @line=();
  my ($zip,$line,$outzipped);
  $zip = $ENV{'VTST_ZIP'} ;
  if($zip eq ''){ $zip = 'gzip'; }
  $energy=$force="na";
  $outzipped=0;
  if((-e "$dir/final/OUTCAR.gz")||(-e "$dir/final/OUTCAR.bz2")||(-e "$dir/final/OUTCAR")){
    if(-e "$dir/final/OUTCAR.gz"){$outzipped=1; system "cd $dir/final; gunzip OUTCAR.gz";}
    if(-e "$dir/final/OUTCAR.bz2"){$outzipped=1; system "cd $dir/final; bunzip2 OUTCAR.bz2";}
    chomp($line=`cd $dir/final; grep 'energy  without entropy=' OUTCAR | tail -1`);
    $energy=VaspEnergy_OUTCAR($line);
    chomp($line=`cd $dir/final; grep 'FORCES: max atom, RMS' OUTCAR | tail -1`);
    $force=VaspForce_OUTCAR($line);
    if($outzipped){system "cd $dir/final; $zip -9 OUTCAR";}
  }elsif((-e "$dir/OUTCAR")||(-e "$dir/OUTCAR.gz")||(-e "$dir/OUTCAR.bz2")){
    if(!(-e "$dir/OUTCAR")){
      if(-e "$dir/OUTCAR.gz"){$outzipped=1; system "cd $dir; gunzip OUTCAR.gz";}
      if(-e "$dir/OUTCAR.bz2"){$outzipped=1; system "cd $dir; bunzip2 OUTCAR.bz2";}
    }
    chomp($line=`cd $dir; grep 'energy  without entropy=' OUTCAR | tail -1`);
    $energy=VaspEnergy_OUTCAR($line);
    chomp($line=`cd $dir; grep 'FORCES: max atom, RMS' OUTCAR | tail -1`);
    $force=VaspForce_OUTCAR($line);
    if($outzipped){system "cd $dir; $zip -9 OUTCAR";}
  }
  return ($energy,$force); 
}

sub VaspEnergy_OUTCAR{
  my $line=shift;
  my $energy="na";
  my @lines=();
  $line=~s/^\s+//;
  if($line ne ""){
    @lines=split(/\s+/,$line);
    if($lines[6]=~/^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$/) {$energy=$lines[6];}
  }
  return $energy;
}

sub VaspForce_OUTCAR{
  my $line=shift;
  my $force="na";
  my @lines=();
  $line=~s/^\s+//;
  if($line ne ""){
    @lines=split(/\s+/,$line);
    if($lines[4]=~/^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$/) {$force=$lines[4];}
  }
  return $force;
}
# ------------------------------------------------------------------------
# Compare two POSCARs [$same=CompareTwoPOSCAR($curspdir,$curpsdir,$Rdiffmax)]
# ------------------------------------------------------------------------
sub CompareTwoPOSCAR{
  my ($curspdir,$curpsdir,$Rdiffmax)=@_;
  my ($coordinates1,$coordinates2,$total_atoms1,$total_atoms2);
  my ($basis,$lattice,$num_atoms,$selectiveflag,$selective,$description,$filetype);
  my ($same, $i, $j, $differ,$poscarfilename,$mag);
  #print " Reading POSCAR1:$curspdir\n";
  $poscarfilename=$curspdir;
  ($coordinates1,$basis,$lattice,$num_atoms,$total_atoms1,$selectiveflag,$selective,$description,$filetype)=read_poscar($poscarfilename);
  set_bc($coordinates1,$total_atoms1);
  #print " Reading POSCAR2:$curpsdir\n";
  $poscarfilename=$curpsdir;
  ($coordinates2,$basis,$lattice,$num_atoms,$total_atoms2,$selectiveflag,$selective,$description,$filetype)=read_poscar($poscarfilename);
  set_bc($coordinates2,$total_atoms2);
  if($total_atoms1 != $total_atoms2){die "In comparing two poscars, the numbers of atoms are not same: $curspdir $curpsdir\n";}
  dirkar($coordinates1,$basis,$lattice,$total_atoms1);
  dirkar($coordinates2,$basis,$lattice,$total_atoms2);
  ($differ,$mag)=minimum_image_dist($coordinates1,$coordinates2,$basis,$lattice,$total_atoms1);
  $same=1;
  $Rdiffmax=$Rdiffmax*$Rdiffmax;
  for($i=0;$i<$total_atoms1;$i++){ # atom by atom
    if($mag->[$i] > $Rdiffmax){
      $same=0;
      return $same;
    }
  }
  return $same;
}

# --------------------------------------------------------------------
# Check if two POSCARs are indistinguishable: return $same and a linker
# --------------------------------------------------------------------
sub indistinguishable{
  my ($poscar1,$poscar2,$Rdiffmax)=@_;
  my ($coordinates1,$coordinates2,$total_atoms1,$total_atoms2);
  my ($basis,$lattice,$num_atoms1,$num_atoms2,$selectiveflag,$selective,$description,$filetype);
  my ($mag,$same,$i,$j,$k,$m,$differ,$matched,$mark,$start,$end,%included);
  my %linker=();
  my @dummy=();
  #print " Reading POSCAR1:$poscar1\n";
  ($coordinates1,$basis,$lattice,$num_atoms1,$total_atoms1,$selectiveflag,$selective,$description,$filetype)=read_poscar($poscar1);
  set_bc($coordinates1,$total_atoms1);
  #print " Reading POSCAR2:$poscar2\n";
  ($coordinates2,$basis,$lattice,$num_atoms2,$total_atoms2,$selectiveflag,$selective,$description,$filetype)=read_poscar($poscar2);
  set_bc($coordinates2,$total_atoms2);
  if($total_atoms1 == $total_atoms2 && @$num_atoms1 == @$num_atoms2) {
    #print "@$num_atoms1\n";
    $mark=@$num_atoms1 - 1;
    $end=$total_atoms1;
    $Rdiffmax=$Rdiffmax*$Rdiffmax;
    if($end < 0) {die "In comparing two POSCARs, there are no atoms in either one.";}
    for($i=$mark;$i>=0;$i--) {#the most flexible atom groups are usually placed at the end of the POSCAR file.
      if($num_atoms1->[$i] == $num_atoms2->[$i]) {
        %included=();
        $start=$end - $num_atoms1->[$i];
        #if($i==0) {print "start at:$start, end at $end\n"};
        for($j=$start;$j<$end;$j++) { #locating a group of atoms
          #print "atom: $j\n";
          for($k=$start;$k<$end;$k++) {
            for($m=0;$m<3;$m++){
              $differ->[0][$m]=pbc($coordinates1->[$j][$m]-$coordinates2->[$k][$m]);
            }
            $differ=dirkar($differ,$basis,$lattice,1);
            $matched=1;
            $mag=0.0;
            for($m=0;$m<3;$m++){ $mag+=$differ->[0][$m]*$differ->[0][$m]; }
            if($mag > $Rdiffmax) {
              $matched=0;
              #print $mag."***";
            }
            #print "\n";
            if($matched) {
              if(!(exists($included{$k}))) {
                $linker{$j}=$k;
                $included{$k}=$j;
                last;
              }
            }
          } # end of for($K ...
          if(!(exists($linker{$j}))) {
            $same=0;
            return ($same,\%linker);
          }
        } # end of for($j ...
        $end=$start;
      }else{
        $same=0;
        return ($same,\%linker);
      } 
    } # end of for($mark... (finishing one group
    if((keys %linker) == $total_atoms1) { # if all atoms from both POSCARS uniquely matched up
      $same=1;
      #sub numerically {$a <=> $b};
      #for $i (sort numerically keys %linker) {print $i.":".$linker{$i}."\n";}
    }
  }else{
    #print "Warning in indistinguishable: the number(s) of atoms and/or types in two POSCARs are not same. Weird!\n";
    $same=0;
    return ($same,\%linker);
  }
  return ($same,\%linker);
}

# --------------------------------------------------------------------------------------
# Reshuffle one POSCAR based on the linker{poscar1}=poscar2 provided by indistinguishable
# --------------------------------------------------------------------------------------
sub ReshufflePOSCAR{
  my ($poscar1,$poscar2,$linker)=@_;
  my ($basis,$lattice,$num_atoms,$selectiveflag,$selective,$description,$total_atoms,$filetype);
  my ($i,$j,$k,$coordinates,$R,$Select);
  ($coordinates,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective,$description,$filetype)=read_poscar($poscar2);
  set_bc($coordinates,$total_atoms);
  if((keys %$linker) != $total_atoms) {
    die "Error in reshuffle: the size of linker is not same as the number of total atoms from $poscar2\n";
  }
  for($i=0;$i<$total_atoms;$i++) {
    for($j=0;$j<3;$j++) {
      $R->[$i][$j]=$coordinates->[$linker->{$i}][$j];
    }       
    $Select->[$i]=$selective->[$linker->{$i}];
  }           
  write_poscar($R,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$Select,$description,$poscar1,$filetype);
}               

# --------------------------------------------------------------------------------------
# Reshuffle MODECARs based on the linker{poscar1}=poscar2 provided by indistinguishable
# --------------------------------------------------------------------------------------
sub ReshuffleMODECAR{
  my ($modecar1,$modecar2,$linker)=@_;
  my ($i,$j,$k,$coordinates,$total_atoms,$R);
  ($coordinates,$total_atoms)=read_othercar($modecar2);
  #set_bc($coordinates,$total_atoms);
  if((keys %$linker) != $total_atoms) {
    die "Error in reshuffle: the size of linker is not same as the number of total atoms from $modecar2\n";
  }
  for($i=0;$i<$total_atoms;$i++) {
    for($j=0;$j<3;$j++) {
      $R->[$i][$j]=$coordinates->[$linker->{$i}][$j];
    }
  }
  write_othercar($R,$total_atoms,$modecar1);
}

# --------------------------------------------------------------------------------------
# Reshuffle POSCARs in a $quench folder based on the linker.
# --------------------------------------------------------------------------------------
sub ReshuffleMins {
  my ($mins,$linker)=@_;
  my $min1=$mins."/min1";
  my $min2=$mins."/min2";
  my $readme=$mins."/readme";
  my (@targets,$poscarfile);
  @targets=($min1."/POSCAR", $min2."/POSCAR");
  if(-e "$min1/final/POSCAR") {push @targets, $min1."/final/POSCAR";}
  if(-e "$min1/final/CONTCAR") {push @targets, $min1."/final/CONTCAR";}
  if(-e "$min2/final/POSCAR") {push @targets, $min2."/final/POSCAR";}
  if(-e "$min2/final/CONTCAR") {push @targets, $min2."/final/CONTCAR";}
  for $poscarfile (@targets) {
    ReshufflePOSCAR($poscarfile,$poscarfile,$linker);
  }     
  #update con and xyz files. I am not sure if we need to do this -- LX
  if(-e "$min1/final/POSCAR") {$poscarfile=`cd "$min1/final"; $Bin/pos2con.pl POSCAR; $Bin/con2xyz.pl POSCAR.con`;}
  if(-e "$min1/final/CONTCAR") {$poscarfile=`cd "$min1/final"; $Bin/pos2con.pl CONTCAR; $Bin/con2xyz.pl CONTCAR.con`;}
  if(-e "$min2/final/POSCAR") {$poscarfile=`cd "$min2/final"; $Bin/pos2con.pl POSCAR; $Bin/con2xyz.pl POSCAR.con`;}
  if(-e "$min2/final/CONTCAR") {$poscarfile=`cd "$min2/final"; $Bin/pos2con.pl CONTCAR; $Bin/con2xyz.pl CONTCAR.con`;}
  open(README, ">$readme");
  print README "Attention: POSCAR and CONTCAR files have been reshuffled, including those in final folder, if any.\n";
  close README;
}

# --------------------------------------------------------------------------------------
# Reshuffle POSCAR_sp and mins in a process folder based on the linker.
# --------------------------------------------------------------------------------------      
sub ReshuffleProcess{
  my ($prdir,$quench,$linker)=@_; 
  my $mins=$prdir."/".$quench;
  my $readme=$prdir."/readme";
  my (@targets,$POSCAR,$MODECAR);
  @targets=($prdir."/POSCAR_sp", $prdir."/POSCAR");
  if(-e "$prdir/final/CONTCAR") {push @targets, $prdir."/final/CONTCAR";}
  if(-e "$prdir/final/POSCAR") {push @targets, $prdir."/final/POSCAR";}
  if(-e "$prdir/final/CENTCAR") {push @targets, $prdir."/final/CENTCAR";}
  for $POSCAR (@targets){ 
    ReshufflePOSCAR($POSCAR,$POSCAR,$linker);
  }
  @targets=($prdir."/MODECAR");
  if(-e "$prdir/final/MODECAR") {push @targets, $prdir."/final/MODECAR";}
  if(-e "$prdir/final/NEWMODECAR") {push @targets, $prdir."/final/NEWMODECAR";}  
  for $MODECAR (@targets){
    ReshuffleMODECAR($MODECAR,$MODECAR,$linker);
  }
  #update con and xyz files. I am not sure if we need to do this -- LX
  if(-e "$prdir/final/POSCAR") {$POSCAR=`cd "$prdir/final"; $Bin/pos2con.pl POSCAR; $Bin/con2xyz.pl POSCAR.con`;} 
  if(-e "$prdir/final/CONTCAR") {$POSCAR=`cd "$prdir/final"; $Bin/pos2con.pl CONTCAR; $Bin/con2xyz.pl CONTCAR.con`;} 
  ReshuffleMins($mins,$linker);
  open(README, ">$readme");
  print README "Attention: POSCAR and MODECAR files have been reshuffled, including those in final folder, if any.\n";
  print README "However, dynmat folders, if any, are not reshuffled since there is no need to do it.\n";
  close README;
} 

# --------------------------------------------------------------------------------------
# Calculate the radiator for one poscar file
#($rad)=Radiator($start,$end,$coordinates,$total_atoms); # RT: Maybe rename this function
# --------------------------------------------------------------------------------------
sub Radiator{
  my ($start,$end,$coordinates,$total_atoms)=@_;
  my ($i,$j,$gcenter,$num,$rad);
  #calculate the center of the geometry
  $num=$end - $start;
  $gcenter->[0]=[(0.0,0.0,0.0)];
  #print "Center:@{$gcenter->[0]}\n";
  for($i=$start;$i<$end;$i++) {
    for($j=0;$j<3;$j++) {
      $gcenter->[0][$j]=$gcenter->[0][$j]+pbc($coordinates->[$i][$j]-$coordinates->[$start][$j]);
    }
  }
  if($num) {
    for($j=0;$j<3;$j++) {
      $gcenter->[0][$j]=$gcenter->[0][$j]/$num;
    }
  }
  for($j=0;$j<3;$j++) {
    $gcenter->[0][$j]=$gcenter->[0][$j]+$coordinates->[$start][$j];
  } # It looks redundant, but pbc is the problem, just as the same care we take in building NBE bands.
  #print "Center before set_bc:@{$gcenter->[0]}\n";
  set_bc($gcenter,1);
  #print "Center after set_bc:@{$gcenter->[0]}\n";
  # calculate the deviation for each atom
  for($i=$start;$i<$end;$i++) {
    for($j=0;$j<3;$j++) {
      $rad->[$i][$j]=pbc($coordinates->[$i][$j]-$gcenter->[0][$j]);
    }
  }  
  set_bc($rad,$total_atoms); # $rad cannot be an accumulator!
  #for($i=$start;$i<$end;$i++) {
  #  for($j=0;$j<3;$j++) {
  #    print "***".$rad->[$i][$j];
  #  }
  #  print "\n";
  #}
  return $rad;
}

# --------------------------------------------------------------------------------------
# $rad_dis_f=CalcRadialDistFunc($coordinates,$total_atoms,$basis,$lattice,$rad_dis_f,$j,$size_bins,$amin)
# --------------------------------------------------------------------------------------
sub CalcRadialDistFunc{
  my ($coordinates,$total_atoms,$basis,$lattice,$rad_dis_f,$j,$size_bins,$amin)=@_;
  my ($i,$m,$diff,$distance,$num_bin,$count);
  $count=int($amin/$size_bins);
  $rad_dis_f->{$j}=[()];
  for($i=0;$i<=$count;$i++) {
    push @{$rad_dis_f->{$j}}, 0;
  }
  for($i=0;$i<$total_atoms;$i++) {# let's focus on the distribution function for each type of atoms first
    if($i != $j) {
      for($m=0;$m<3;$m++) {
        $diff->[0][$m]=pbc($coordinates->[$i][$m] - $coordinates->[$j][$m]);
      }
      dirkar($diff,$basis,$lattice,1);
      $distance=magnitude($diff,1);
      $num_bin=int($distance/$size_bins);
      $rad_dis_f->{$j}[$num_bin]++;
    }
  }
  return $rad_dis_f;
}

# --------------------------------------------------------------------------------------
# Compare the radial distribution function for two atoms from different coordinates
# ($same,$rad_dis_f1,$rad_dis_f2)=CheckRadialDistFunc($coordinates1,$coordinates2,$total_atoms1,$basis,$lattice,$rad_dis_f1,$j,$rad_dis_f2,$k,$amin);
# --------------------------------------------------------------------------------------
sub CheckRadialDistFunc{
  my ($coordinates1,$coordinates2,$total_atoms1,$basis,$lattice,$rad_dis_f1,$j,$rad_dis_f2,$k,$amin)=@_;
  my $size_bins=0.2;
  my ($i,$num_bins,$same);
  if(!exists($rad_dis_f1->{$j})) {$rad_dis_f1=CalcRadialDistFunc($coordinates1,$total_atoms1,$basis,$lattice,$rad_dis_f1,$j,$size_bins,$amin);}
  if(!exists($rad_dis_f2->{$k})) {$rad_dis_f2=CalcRadialDistFunc($coordinates2,$total_atoms1,$basis,$lattice,$rad_dis_f2,$k,$size_bins,$amin);}
  $num_bins=@{$rad_dis_f1->{$j}};
  if($num_bins != @{$rad_dis_f2->{$k}}) {die "In CheckRadialDistFunc: fatal error, boxes are of different sizes\n";}
  $same=1;
  for($i=0;$i<$num_bins;$i++) {
    if($rad_dis_f1->{$j}[$i] != $rad_dis_f2->{$k}[$i]) {
      $same=0;
      last;
    }
  }
  return ($same,$rad_dis_f1,$rad_dis_f2);
}

# ------------------------------------------------------------------------
# Standard deviation of a list of values
# ------------------------------------------------------------------------
sub standard_deviation{
  my $listvalues=shift;
  my ($i,$average,$accumulator,$length,$similarity);
  $average=$accumulator=0.0;
  $length=@$listvalues;
  if(!$length) {die "In standard deviation, a zero-length list was found\n";}
  if($length==1) {return 1.0;}
  for $i (@$listvalues) {
    $average+=$i;
  }
  $average=$average/$length;
  for $i (@$listvalues) {
    $accumulator+=($i-$average)**2;
  }
  $similarity=sqrt($accumulator)/($length-1);
  return $similarity;
}

# ------------------------------------------------------------------------
# How parallel: $simularity=CalcParallel($coordinates1,$coordinates2,\%linker);
# ------------------------------------------------------------------------
sub CalcParallel{
  my ($coordinates1,$coordinates2,$linker)=@_;
  my $simularity=0;
  my ($j,$k,$m,$differ,$distance);
  my @dummy=();
  for $j (keys %$linker) {
    $k=$linker->{$j};
    for($m=0;$m<3;$m++) {
      $differ->[0][$m]=pbc($coordinates1->[$j][$m]-$coordinates2->[$k][$m]);
    }
    $distance=magnitude($differ, 1);
    push @dummy, $distance;
  }
  $similarity=standard_deviation(\@dummy);
  if($similarity < 1.0e-8) {$similarity=1.0e-8;}
  $similarity=1.0/$similarity;
  return $similarity;
}


# ------------------------------------------------------------------------
# Evaluate two POSCARs [($thatisit, similarity)=EvaluateTwoPOSCAR($POSCAR,$POSCAR_j,$Rdiffmax)]
# ------------------------------------------------------------------------
sub EvaluateTwoPOSCAR{
  my ($POSCAR,$POSCAR_j,$Rdiffmax,)=@_;
  my ($coordinates1,$coordinates2,$total_atoms1,$total_atoms2,$num_atoms1,$num_atoms2,$selective1,$selective2,$filetype);
  my ($basis,$lattice,$selectiveflag,$description,$frozen,$numfrozen,$vector1,$vector2,$vector3);
  my ($i,$j,$k,$m,$differ,$matched,$mark,$start,$end,%included,$rad_dis_f1,$rad_dis_f2);
  my %linker=();
  my @dummy=();
  my ($thatisit,$similarity,$samei,$checkedatoms,$dij_max);
  $thatisit="";
  $similarity=0;
  $numfrozen=0;
  $checkedatoms=0;
  #print " Reading POSCAR1:$POSCAR\n";
  $poscarfilename=$POSCAR;
  ($coordinates1,$basis,$lattice,$num_atoms1,$total_atoms1,$selectiveflag,$selective1,$description,$filetype)=read_poscar($poscarfilename);
  set_bc($coordinates1,$total_atoms1);
  #print " Reading POSCAR2:$curpsdir\n";
  $poscarfilename=$POSCAR_j;
  ($coordinates2,$basis,$lattice,$num_atoms2,$total_atoms2,$selectiveflag,$selective2,$description,$filetype)=read_poscar($poscarfilename);
  set_bc($coordinates2,$total_atoms2);
  if($total_atoms1 == $total_atoms2 && @$num_atoms1 == @$num_atoms2) {
    #print "@$num_atoms1\n";
    $mark=@$num_atoms1 - 1;
    $end=$total_atoms1;
    if($end < 0) {die "In comparing two POSCARs, there are no atoms in either one.";}
    for($j=0;$j<3;$j++){
      $vector1->[0][$j]=$basis->[$j][0];
      $vector2->[0][$j]=$basis->[$j][1];
      $vector3->[0][$j]=$basis->[$j][2];
    }
    $dummy[0]=magnitude($vector1,1);
    $dummy[1]=magnitude($vector2,1);
    $dummy[2]=magnitude($vector3,1);
    $dij_max=sqrt($dummy[0]*$dummy[0]+$dummy[1]*$dummy[1]+$dummy[2]*$dummy[2])/2.0;
    #print "dij_max=$dij_max\n";
    for($i=$mark;$i>=$mark;$i--) {#the most flexible atom groups are usually placed at the end of the POSCAR file.
      if($num_atoms1->[$i] == $num_atoms2->[$i]) {
        $checkedatoms=$checkedatoms+$num_atoms1->[$i];
        %included=();
        $start=$end - $num_atoms1->[$i];
        #if($i==0) {print "start at:$start, end at $end\n"};
        $rad1=Radiator($start,$end,$coordinates1,$total_atoms1);
        $rad2=Radiator($start,$end,$coordinates2,$total_atoms2);
        for($j=$start;$j<$end;$j++) { # locating a group of atoms
          #print "pair: $j\n";
          @dummy=split(/\s+/, $selective1->[$j]);
          $frozen=join("", @dummy);
          #print "$frozen\n";
          if(lc($frozen) eq "fff") { # skip the frozen atom. You make sure the frozen atoms are same through the simulation
            $numfrozen++;
            next;
          }
          for($k=$start;$k<$end;$k++) {
            for($m=0;$m<3;$m++){
              $differ->[0][$m]=pbc($rad1->[$j][$m]-$rad2->[$k][$m]);
            }
            $differ=dirkar($differ,$basis,$lattice,1);
            $matched=1;
            for($m=0;$m<3;$m++){
              #print $differ->[0][$m]."***";
              if(abs($differ->[0][$m]) > $Rdiffmax) {
                $matched=0;
                last;
              }
            }
            #print "\n";
            if($matched) {
              if(!(exists($included{$k}))) {
                ($same,$rad_dis_f1,$rad_dis_f2)=CheckRadialDistFunc($coordinates1,$coordinates2,$total_atoms1,$basis,$lattice,$rad_dis_f1,$j,$rad_dis_f2,$k,$dij_max);
                #print "@{$rad_dis_f1->{$j}}\n";
                #print "@{$rad_dis_f2->{$k}}\n";
                #print "$j and $k are same? $same\n";  
                if($same) {
                  $linker{$j}=$k;
                  $included{$k}=$j;
                  last;
                }
              }
            }
          } # end of for($K ...
          if(!(exists($linker{$j}))) {
            $similarity=0;
            return ($thatisit,$similarity,\%linker);
          }
        } # end of "for($j ..."
        $end=$start;
      }else{
        $similarity=0;
        return ($thatisit,$similarity,\%linker);
      } 
    } # end of for($mark... (finishing one group
    if((keys %linker) == ($checkedatoms-$numfrozen)) { 
      # if all free atoms from both POSCARS uniquely matched up. You must guarantee that all the frozen parts are same and consistent through the simulation!
      $thatisit=$POSCAR_j;
      # choose the one with the best parallel quality (this part can be upgraded to be more general if we use some transformation matrix)
      $similarity=CalcParallel($coordinates1,$coordinates2,\%linker);
      #sub numerically {$a <=> $b};
      #for $i (sort numerically keys %linker) {print $i.":".$linker{$i}."\n";}
    }
  }else{
    #print "Warning in indistinguishable: the number(s) of atoms and/or types in two POSCARs are not same. Weird!\n";
    $similarity=0;
    return ($thatisit,$similarity,\%linker);
  }
  return ($thatisit,$similarity,\%linker);
}

# --------------------------------------------------------------------
# Write state information to st.dat [WriteSt($curdir,$StFile,$st)]
# --------------------------------------------------------------------
sub WriteSt{
  my ($curdir,$StFile,$st)=@_;
  my $stfilename=$curdir."/".$StFile;
  my ($energy, $force, $dummy);
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

# -----------------------------------------------------------------------
# Read state information from st.dat [%st=ReadSt($curdir,$StFile)]
# -----------------------------------------------------------------------
sub ReadSt{
  my ($curdir,$StFile)=@_;
  my $stfilename=$curdir."/".$StFile;
  my %st=();
  my $line="";
  my @line=();
  my ($i, $j, $numprst);
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
    if($line eq "" || $line eq "\n"){ $i--; next; } # skip empty lines
    @line=split(/\s+/,$line);
    if($i==1){# the first line must be st and st energy etc.
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
        $st{lc($line[0])}[2] = $st{lc($line[0])}[2] + $st{$curdir}[0]; # convert barriers into energy.
      }
    }
  }
  if($st{"numprst"}[0] != $numprst){die "The numbers of st00xx states conflict in $stfilename\n";}
  close ST;
  if(exists($st{""})) {delete $st{""}};
  return(\%st);
}

# ---------------------------------------------------------------------------------------------------------
# relax an initial saddle point guess to prevent atoms from getting too close to each other
# spring_relaxation($coordinates,$basis,$lattice,$totalatoms,$selective);
# ---------------------------------------------------------------------------------------------------------
sub spring_relaxation{
  my ($R,$basis,$lattice,$totalatoms,$selective)=@_;
  my $epsilon=0.1;
  my $stepmax=10000;
  my $cutoff=0.2;
  my $drmax=0.2;
  my @line=();
  my ($i,$j,$k,$difference,$Rij,$lattice_vec,@BOX,@BOX_,@HALFBOX,@HALFBOX_,$ORTHOGONAL);
  my ($force,$fij,$phiij,$step,$goodsound,$converged,$frozen,$dr,$rcut);
  for($i=0;$i<3;$i++){
    for($j=0;$j<3;$j++){
      $lattice_vec->[$i][$j]=$basis->[$j][$i];
    }
  }
  $ORTHOGONAL=0;
  if($ARGV[0] == 1){$ORTHOGONAL=1};
  #----------------------------------
  #Some variables for ortogonal boxes
  if($ORTHOGONAL){
    @BOX=(abs($lattice_vec->[0][0]),abs($lattice_vec->[1][1]),abs($lattice_vec->[2][2]));
    @BOX_=(-1.0*$BOX[0], -1.0*$BOX[1], -1.0*$BOX[2]);
    @HALFBOX=($BOX[0]*0.5, $BOX[1]*0.5, $BOX[2]*0.5);
    @HALFBOX_=(-0.5*$BOX[0], -0.5*$BOX[1], -0.5*$BOX[2]);
  }
  $rcut=$cutoff*$cutoff;
  # figure out frozen flags for each degree of freedom
  for($i=0;$i<$totalatoms;$i++) {
    @line=split(/\s+/,$selective->[$i]);
    for($j=0;$j<3;$j++) {
      $frozen->[$i][$j]=$line[$j];
    }
  }
  # start the steepest descent loop
  $converged=0;
  for($step=0;$step<$stepmax;$step++) {
    if($converged) {last;}
    print "Step No. $step\n";
    for($i=0;$i<$totalatoms;$i++) {
      for($j=0;$j<3;$j++) {
        $force->[$i][$j]=0.0;
      }
    }
    $goodsound=1;
    for($i=0;$i<$totalatoms-1;$i++) {
      for($j=$i+1;$j<$totalatoms;$j++) {
        for($k=0;$k<3;$k++){
          $difference->[$k]=$R->[$i][$k]-$R->[$j][$k];
        }
        ($difference, $Rij)=FindMinimumImage($difference,$lattice_vec,\@BOX,\@HALFBOX,\@HALFBOX_);
        if($Rij < $rcut){
          $goodsound=0;
          $Rij=sqrt($Rij);
          ($fij,$phiij)=linear_repulsion_pot($Rij,$cutoff);
          #($fij,$phiij)=exp_repulsion_pot($Rij,$cutoff);
          $fij=$fij/$Rij;
          for($k=0;$k<3;$k++){
            $phiij=$fij*$difference->[$k];
            $force->[$i][$k]-=$phiij;
            $force->[$j][$k]+=$phiij;
          }
        }
      }
    }
    if($goodsound) {
      $converged=1;
      next;
    }
    for($i=0;$i<$totalatoms;$i++) {
      for($j=0;$j<3;$j++) {
        if(lc($frozen->[$i][$j]) eq "f") {next;}
        $dr=$epsilon*$force->[$i][$j];
        if(abs($dr) > $drmax) {$dr=$drmax*$dr/abs($dr);}
        #print "$i: $j -- dr: $dr\n"; 
        $R->[$i][$j]+=$dr;
      }
    }
    if($ORTHOGONAL){
      for($i=0;$i<$totalatoms;$i++){
        for($j=0;$j<3;$j++){
          while($R->[$i][$j] > $BOX[$j]){ $R->[$i][$j]-=$BOX[$j]; }
          while($R->[$i][$j] < $BOX_[$j]){ $R->[$i][$j]+=$BOX[$j]; }
        }
      }
    }else{
      $R=kardir($R,$basis,$lattice,$total_atoms);
      set_bc($R,$total_atoms);
      $R=dirkar($R,$basis,$lattice,$total_atoms);
    }
  }
  if($converged == 1) {
    print "fully relaxed after $step (1: trivial) steps\n";
  }else {
    print "Not fully relaxed after $stepmax steps, but we can;t wait. move on\n";
  }
}

sub linear_repulsion_pot{ # E=kx-b, k< 0
  my ($r,$cutoff)=@_;
  my $slope=-0.3;
  my ($energy, $fij);
  my $intercept=$slope*$cutoff;
  $energy=$slope*$r-$intercept;
  $fij=$slope;
  return ($fij,$energy);
} 

#----------------------------------------------------------------------
# FindMinimumImage: (INPUT) A vector (pointer) (OUTPT) the minimum image vector
#----------------------------------------------------------------------
sub FindMinimumImage{
  my ($dr_car,$lattice_vec,$ORTHOGONAL,$BOX,$HALFBOX,$HALFBOX_)=@_;
  my ($dsqmin,$done,$found,$v1,$v2,$v3,$drt_car,$dsq);
  my ($i,$j,$k);
  if($ORTHOGONAL){
   for($i=0;$i<@$dr_car;$i++){
     if($dr_car->[$i]<$HALFBOX_->[$i]){ $dr_car->[$i]+=$BOX->[$i]; }
     if($dr_car->[$i]>$HALFBOX->[$i]){ $dr_car->[$i]-=$BOX->[$i]; }
   }
   $dsqmin=VecDot($dr_car,$dr_car);
  }else{
   $dsqmin=VecDot($dr_car,$dr_car);
   $done=1;
   while($done){
    $found=1;
    for($i=-1;$i<2;$i++){
      $v1=VecMultiply($lattice_vec->[0],$i);
      for($j=-1;$j<2;$j++){
         $v2=VecMultiply($lattice_vec->[1],$j);
         for($k=-1;$k<2;$k++){
           $v3=VecMultiply($lattice_vec->[2],$k);
           $drt_car=VecSum($dr_car,$v1,$v2,$v3);
           $dsq=VecDot($drt_car,$drt_car);
           if($dsq<$dsqmin){
             $dr_car=$drt_car;
             $dsqmin=$dsq;
             $found=0;
           }
         }
      }
    }
    if($found){ $done=0 ; }
   }
  }
  #print sqrt($dsqmin)."\n";
  return ($dr_car,$dsqmin);
}

#----------------------------------------------------------------------
# VecDot: (INPUT) two vectors (pointers) (OUTPUT) a dot product
#----------------------------------------------------------------------
sub VecDot{
  my ($vec1, $vec2)=@_;
  my ($i, $dot);
  #check dimensions before calling this sub
  $dot=0.0;
  for($i=0;$i<@$vec1;$i++){
    $dot+=$vec1->[$i]*$vec2->[$i];
  }
  return $dot;
}

#----------------------------------------------------------------------
# VecSum: (INPUT) 3-d vectors (pointers) (OUTPUT) a sum vector (pointer)
#----------------------------------------------------------------------
sub VecSum{
  my ($dimension, $i, $j, @sum);
  $dimension=@_;
  @sum=(0.0,0.0,0.0);
  for($i=0;$i<$dimension;$i++){
    for($j=0;$j<3;$j++){
      $sum[$j]+=$_[$i]->[$j];
    }
  }
  return \@sum;
}

#----------------------------------------------------------------------
# VecMultiply: (INPUT) one vector (pointer) (OUTPUT) a new vector (pointer)
#----------------------------------------------------------------------
sub VecMultiply{
  my ($vec, $fac)=@_;
  my ($i, $vecout);
  for($i=0;$i<@$vec;$i++){
    $vecout->[$i]=$vec->[$i]*$fac;
  }
  return $vecout;
}

#----------------------------------------------------------------------
# ($differ,$mag)=minimum_image_dist($R1,$R2,$basis,$lattice,$total_atoms);
#----------------------------------------------------------------------
sub minimum_image_dist{
  my ($R1,$R2,$basis,$lattice,$total_atoms)=@_;
  my ($i,$j,$lattice_vec,@BOX,@HALFBOX,@HALFBOX_,$ORTHOGONAL);
  my ($differ,$mag);
  for($i=0;$i<3;$i++){
    for($j=0;$j<3;$j++){
      $lattice_vec->[$i][$j]=$basis->[$j][$i];
    }
  }
  $ORTHOGONAL=0;
  if($ARGV[0] == 1){$ORTHOGONAL=1};
  #----------------------------------
  #Some variables for ortogonal boxes
  if($ORTHOGONAL){
    @BOX=(abs($lattice_vec->[0][0]),abs($lattice_vec->[1][1]),abs($lattice_vec->[2][2]));
    @HALFBOX=($BOX[0]*0.5, $BOX[1]*0.5, $BOX[2]*0.5);
    @HALFBOX_=(-0.5*$BOX[0], -0.5*$BOX[1], -0.5*$BOX[2]);
  }
  for($i=0;$i<$total_atoms;$i++) {
    for ($j=0; $j<3; $j++) {
      $differ->[$i][$j]=$R1->[$i][$j]-$R2->[$i][$j];
    }
    ($differ->[$i],$mag->[$i])=FindMinimumImage($differ->[$i],$lattice_vec,$ORTHOGONAL,\@BOX,@HALFBOX,@HALFBOX_);
  }
  return ($differ,$mag);
}

# ---------------------------------------------------------------------------------------------------------
# ($recovered,$whereabout,$st) = GrabQuenchedSP_OneSt($target,$curdir,$prdir,$quench,$dynmat,$turnoff_dynmat,
#               $Equivalency,$SearchesAlgo,$Ediffmax,$Rdiffmax,$DynMatFile,$stfile,$st);
# ---------------------------------------------------------------------------------------------------------
sub GrabQuenchedSP_OneSt{
  my ($target,$curdir,$prdir,$quench,$dynmat,$turnoff_dynmat,$Equivalency,$SearchesAlgo,$Ediffmax,$Rdiffmax,$DynMatFile,$stfile,$st)=@_;
  my $curspdir=$curdir."/".$prdir;
  my ($energy,$numprst,$j,$status,$quality,$curpsdir1,$curpsdir2,$same,$recovered,$whereabout,$targetpr,$targetst,$stinfo,$minj);
  my %minsinfo=();
  my ($linker,$yes);
  $energy=$st->{$prdir}[2];
  $targetst=substr($target, -6);
  #print "$target: $targetst\n";
  $recovered=0;
  $stinfo=ReadSt($target,$stfile);
  $numprst=$stinfo->{"numprst"}[0];
  for($j=1;$j<=$numprst;$j++){
    $targetpr=DirName("pr", $j);
    $whereabout="";
    $yes=0;
    # just compare it with any good or bad saddle which may have been quenched
    $status=lc($stinfo->{$targetpr}[0]);
    $quality=lc($stinfo->{$targetpr}[1]);
    if($quality eq "good" || (($quality eq "bad") && ($status eq "done"))){
      if($status eq "done" || $status eq "dynmat") { # as long as a then-good saddle was quenched.
        $curpsdir1=$target."/".$targetpr."/POSCAR_sp";
        $curpsdir2=$curspdir."/POSCAR_sp";
        if(abs($energy - $stinfo->{$targetpr}[2]) < $Ediffmax) {
          $same=CompareTwoPOSCAR($curpsdir2,$curpsdir1,$Rdiffmax); #compare two POSCAR_sp files.
          if($Equivalency && !$same) {
            ($same,$linker)=indistinguishable($curpsdir2,$curpsdir1,$Rdiffmax);
             $yes=$same;
          } 
          if($same) {
            #print "It is same as $curpsdir1\n";
            $whereabout=$curpsdir1;
            if((-e "$target/$targetpr/$quench/min1/final/CONTCAR") && (-e "$target/$targetpr/$quench/min2/final/CONTCAR") && ((-e "$target/$targetpr/$quench/min1/final/OUTCAR.gz") || (-e "$target/$targetpr/$quench/min1/final/OUTCAR.bz2") || (-e "$target/$targetpr/$quench/min1/final/OUTCAR")) && ((-e "$target/$targetpr/$quench/min2/final/OUTCAR.gz") || (-e "$target/$targetpr/$quench/min2/final/OUTCAR.bz2") || (-e "$target/$targetpr/$quench/min2/final/OUTCAR"))) {
              system "cp -r $target/$targetpr/$quench $curspdir/";
              if($Equivalency && $yes) {
                #print "In $curspdir: reshuffle of mins POSCARs from $target/$targetpr/$quench just happened.\n";
                ReshuffleMins($curspdir."/".$quench,$linker);
              }
              if($quality eq "good"){
                $minj=substr($stinfo->{$targetpr}[3], -4);
                $minsinfo{$stinfo->{$targetpr}[3]}=$stinfo->{$targetpr}[4];
                if($minj eq "min2") {
                  $minsinfo{"$quench/min1"}=$stinfo->{$targetst}[0];
                }elsif($minj eq "min1") {
                  $minsinfo{"$quench/min2"}=$stinfo->{$targetst}[0];
                }else{
                  #print "Something wrong with the mins information at $target. Ignore this one, but better go check it.\n";
                  next;
                }
              }else{ #it was a bad saddle
                $minsinfo{"$quench/min1"}=$stinfo->{$targetpr}[3];
                $minsinfo{"$quench/min2"}=$stinfo->{$targetpr}[4];
              }
              $st=EvalSP($curdir,$prdir,$quench,$SearchesAlgo,\%minsinfo,$Ediffmax,$Rdiffmax,$st);
              if($st->{$prdir}[1] eq "good") {
                if($turnoff_dynmat) {
                  $st->{$prdir}[0]="done";
                }else{
                  if(-e  "$target/$targetpr/$dynmat/final/$DynMatFile") {
                    system "cp -r $target/$targetpr/$dynmat $curspdir/";
                    $st->{$prdir}[0]="done";
                  }else{
                    $st->{$prdir}[0]="dynmat";
                  }
                }
              }else{
                $st->{$prdir}[0]="done";
              } #end of the "if($st->{$prdir}[1] eq "good")" block
              $recovered = 1;
              last;
            }else{
              #print "Notice:  $curpsdir1 is broken. Skip this one and check others.\n";
            }
          }#end of the "if($same)" block  
        }# end of the "abs($energy - $stinfo->{$targetpr}[2]) < $Ediffmax" block
      }
    }# end of the "if($quality eq "good")" block
  }#end of the "for" loop
  return ($recovered,$whereabout,$st);
}
# ---------------------------------------------------------------------------------------------------------
# Check if the new saddle has been quenched in the previous processes, or previous states or any state in the pool
# ($recycled,$st)=GrabQuenchedSP($curdir,$prdir,$SearchesAlgo,$quench,$dynmat,$stfile,$turnoff_dynmat,
#                $Equivalency,$stpool,$Ediffmax,$Rdiffmax,$DynMatFile,$stfile,$stdirs,$st);
# ---------------------------------------------------------------------------------------------------------
sub GrabQuenchedSP{
  my ($curdir,$prdir,$SearchesAlgo,$quench,$dynmat,$stfile,$turnoff_dynmat,$Equivalency,$stpool,$Ediffmax,$Rdiffmax,$DynMatFile,$stfile,$stdirs,$st)=@_;
  my ($filepath,$stpoolsize,$POSCAR_sp,$i,$j);
  my ($curspdir,$target,$recovered,$whereabout);
  my ($numst,$headpath);
  my @inmaindir=();
  my @stdir=();
  $curspdir=$curdir."/".$prdir;
  $numst=@$stdirs;
  $POSCAR_sp=$curspdir."/POSCAR_sp";
  $recovered=0;
  $whereabout="";
  for($j=($numst-2);$j>=0;$j--) {
    $target=$stdirs->[$j];
    #print "recovered= $recovered\n";
    ($recovered,$whereabout,$st) = GrabQuenchedSP_OneSt($target,$curdir,$prdir,$quench,$dynmat,$turnoff_dynmat,$Equivalency,$SearchesAlgo,$Ediffmax,$Rdiffmax,$DynMatFile,$stfile,$st);
    if($recovered) {return ($recovered,$whereabout,$st);}
  }
  $stpoolsize=@$stpool;
  for($i=0;$i<$stpoolsize;$i++) {
    $headpath=$stpool->[$i];
    opendir MAINSTDIR, $headpath or die "can't open this dir!";
    @inmaindir=readdir MAINSTDIR;
    @stdir = grep /^st\d{4}$/, @inmaindir;
    closedir MAINSTDIR;
    $numst=@stdir;
    for($j=0;$j<$numst;$j++){
      $target=$headpath."/".$stdir[$j];
      ($recovered,$whereabout,$st) = GrabQuenchedSP_OneSt($target,$curdir,$prdir,$quench,$dynmat,$turnoff_dynmat,$Equivalency,$SearchesAlgo,$Ediffmax,$Rdiffmax,$DynMatFile,$stfile,$st);
      if($recovered) {return ($recovered,$whereabout,$st);}
    }
  }
  return ($recovered,$whereabout,$st);
}

# ---------------------------------------------------------------------------------------------------------
# Look back one pr000x folder 
# ($recovered,$included,$stprs,$st) = LookBackOneSt($target,$stdir,$previousst,$energy,$quench,$recovered,$included,$Equivalency,$Ediffmax,$Rdiffmax,$stfile,$stprs,$st)
# ---------------------------------------------------------------------------------------------------------
sub LookBackOneSt{
  my ($target,$stdir,$previousst,$energy,$quench,$recovered,$included,$Equivalency,$Ediffmax,$Rdiffmax,$stfile,$stprs,$st)=@_;
  my ($pr000i,$POSCAR_i,$POSCAR,$same,$repeat,$repeated,$newprdir,$filepath,$sp1,$sp2,$stinfo,$i,$k,$temp,$linker);
  my $numprstinfo=0;
  my $yes=0;
  #print "current st: $stdir  target st: $target recovered: $recovered quench: $quench stfile $stfile energy: $energy\n";
  print "$stdir LookBackOneSt $target to recover good saddles (if any). Before this, $recovered saddles have been recovered.\n";
  $POSCAR=$stdir."/POSCAR";
  $stinfo=ReadSt($target,$stfile);
  if(!(exists($stinfo->{"numprst"}))) {
    print "Warning in $target, $stfile is broken, no numprst. So,set it to zero.\n";
    $stinfo->{"numprst"}=[(0)];
  }
  $numprstinfo=$stinfo->{"numprst"}[0];
  for($i=1;$i<=$numprstinfo;$i++) {
    $pr000i=DirName("pr",$i);
    $temp="";
    if((lc($stinfo->{$pr000i}[1]) eq "good") && (lc($stinfo->{$pr000i}[0]) eq "done" || lc($stinfo->{$pr000i}[0]) eq "dynmat")) {
      if(abs($energy-$stinfo->{$pr000i}[4]) < $Ediffmax) {
        $POSCAR_i=$target."/".$pr000i."/".$stinfo->{$pr000i}[3]."/POSCAR";
        $same=CompareTwoPOSCAR($POSCAR_i,$POSCAR,$Rdiffmax);
        if($Equivalency && !$same) {
          ($same,$linker)=indistinguishable($POSCAR,$POSCAR_i,$Rdiffmax);
          $yes=$same;
        }
        if($same) {
          if($previousst eq $target) {$included->{$pr000i}="included";} # mark this process as included.
          $repeated=0;
          for($k=0;$k<@$stprs;$k++) {
            if(abs($stinfo->{$pr000i}[2] - $st->{$stprs->[$k]}[2]) < $Ediffmax) {
              $sp1=$stdir."/".$stprs->[$k]."/POSCAR_sp";
              $sp2=$target."/".$pr000i."/POSCAR_sp";
              if($Equivalency && $yes) {
                $temp=$target."/".$pr000i."/POSCAR_sp_reshuffled";
                ReshufflePOSCAR($temp,$sp2,$linker);
                $sp2=$temp;
              }
              $repeat=CompareTwoPOSCAR($sp1,$sp2,$Rdiffmax);
              if($temp){system "rm $temp";}
              if($repeat) {
                $repeated=1;
                last;
              }  
            }
          }
          if($repeated) {next;}
          $recovered++;
          $filepath=$target."/".$pr000i;
          $newprdir=DirName("pr",$recovered);
          system "cp -r $filepath $stdir/$newprdir; cd $stdir; cd $newprdir; cd $quench; mv min1 mininter; mv min2 min1; mv mininter min2";
          if($Equivalency && $yes) {ReshuffleProcess($stdir."/".$newprdir,$quench,$linker);}
          $st->{$newprdir}=$stinfo->{$pr000i};
          $st->{$newprdir}[4]=$stinfo->{$target}[0];
          push @{$stprs}, $newprdir;
        }
      }
    }elsif((lc($stinfo->{$pr000i}[1]) eq "bad") && (lc($stinfo->{$pr000i}[0]) eq "done")){
      if(!(-e $target."/".$pr000i."/".$quench)) {next;} # this bad saddle was not quenched; probably a negative barrier when it was filtered.
      if(abs($energy-$stinfo->{$pr000i}[4]) < $Ediffmax) {
        $POSCAR_i=$target."/".$pr000i."/".$quench."/min2/POSCAR";
        $same=CompareTwoPOSCAR($POSCAR_i,$POSCAR,$Rdiffmax);
        if($Equivalency && !$same) {
          ($same,$linker)=indistinguishable($POSCAR,$POSCAR_i,$Rdiffmax);
          $yes=$same;
        }
        if($same) {
          if($previousst eq $target) {$included->{$pr000i}="included";} # mark this process as included.
          $repeated=0;
          for($k=0;$k<@$stprs;$k++) {
            if(abs($stinfo->{$pr000i}[2] - $st->{$stprs->[$k]}[2]) < $Ediffmax) {
              $sp1=$stdir."/".$stprs->[$k]."/POSCAR_sp";
              $sp2=$target."/".$pr000i."/POSCAR_sp";
              if($Equivalency && $yes) {
                $temp=$target."/".$pr000i."/POSCAR_sp_reshuffled";
                ReshufflePOSCAR($temp,$sp2,$linker);
                $sp2=$temp;
              }
              $repeat=CompareTwoPOSCAR($sp1,$sp2,$Rdiffmax);
              if($temp){system "rm $temp";}
              if($repeat) {
                $repeated=1;
                last;
              }  
            }
          }
          if($repeated) {next;}
          $recovered++;
          $filepath=$target."/".$pr000i;
          $newprdir=DirName("pr",$recovered);
          system "cp -r $filepath $stdir/$newprdir";
          if($Equivalency && $yes) {ReshuffleProcess($stdir."/".$newprdir,$quench,$linker);}
          $st->{$newprdir}=$stinfo->{$pr000i};
          $st->{$newprdir}[1]="good";
          $st->{$newprdir}[4]=$st->{$newprdir}[3];
          $st->{$newprdir}[3]=$quench."/min1";
          push @{$stprs}, $newprdir;
          next;
        }
      }
      if(abs($energy-$stinfo->{$pr000i}[3]) < $Ediffmax) {
        $POSCAR_i=$target."/".$pr000i."/".$quench."/min1/POSCAR";
        $same=CompareTwoPOSCAR($POSCAR_i,$POSCAR,$Rdiffmax);
        if($Equivalency && !$same) {
          ($same,$linker)=indistinguishable($POSCAR,$POSCAR_i,$Rdiffmax);
          $yes=$same;
        }
        if($same) {
          if($previousst eq $target) {$included->{$pr000i}="included";}
          $repeated=0;
          for($k=0;$k<@$stprs;$k++) {
            if(abs($stinfo->{$pr000i}[2] - $st->{$stprs->[$k]}[2]) < $Ediffmax) {
              $sp1=$stdir."/".$stprs->[$k]."/POSCAR_sp";
              $sp2=$target."/".$pr000i."/POSCAR_sp";
              if($Equivalency && $yes) {
                $temp=$target."/".$pr000i."/POSCAR_sp_reshuffled";
                ReshufflePOSCAR($temp,$sp2,$linker);
                $sp2=$temp;
              } 
              $repeat=CompareTwoPOSCAR($sp1,$sp2,$Rdiffmax);
              if($temp){system "rm $temp";}
              if($repeat) {
                $repeated=1;
                last;
              }  
            }
          }
          if($repeated) {next;}
          $recovered++;
          $filepath=$target."/".$pr000i;
          $newprdir=DirName("pr",$recovered);
          system "cp -r $filepath $stdir/$newprdir";
          if($Equivalency && $yes) {ReshuffleProcess($stdir."/".$newprdir,$quench,$linker);}
          $st->{$newprdir}=$stinfo->{$pr000i};
          $st->{$newprdir}[1]="good";
          $st->{$newprdir}[3]=$quench."/min2";
          push @{$stprs}, $newprdir;
        }
      }        
    }else{
      next;
    }
  }
  return ($recovered,$included,$stprs,$st);
}

# ---------------------------------------------------------------------------------------------------------
# Make use of good and bad saddles at the previous steps and in the state pool. 
# ($recovered,$included,$stprs,$jobs,$st)= UseGoodBad($stdir,$curstep,$numsaddles,$SearchesAlgo,$quench,
# $stpool,$stenergyfile,$stfile,$recovered,$included,$Equivalency,$Ediffmax,$Rdiffmax,$stprs,$akmc_step,$jobs,$st)
# ---------------------------------------------------------------------------------------------------------
sub UseGoodBad{
  my $stdir=shift;
  my $curstep=shift;
  my $numsaddles=shift;
  my $SearchesAlgo=shift;
  my $quench=shift;
  my $stpool=shift;
  my $stenergyfile=shift;
  my $stfile=shift;
  my $recovered=shift;
  my $included=shift;
  my $Equivalency=shift;
  my $Ediffmax=shift;
  my $Rdiffmax=shift;
  my $stprs=shift;
  my $akmc_step=shift;
  my $jobs=shift;
  my $st=shift;
  my @dummy=();
  my ($stenergy,$stinfo,$filepath,$stpoolsize,$energy,$same,$POSCAR,$POSCAR_i,$curst,$numst,@totalst,$i,$j,$k);
  my ($newprdir,$pr000i,$target,$sp1, $sp2,$repeated,$repeat,$previousst,$totalsteps,$headpath);
  $energy=$akmc_step->{$curstep}[3];
  $POSCAR=$stdir."/POSCAR";
  $stpoolsize=@$stpool;
  $stenergy=ReadStEnergyFile($stenergyfile);
  @totalst=(reverse sort keys %$stenergy);
  $numst=@totalst;
  if($numst > 0){
    $curst=$totalst[0]; # it has been sorted reversely
    # print $curst."  ".$numst."\n";
    ($curst eq DirName("st",$numst)) || die " Error in UseGoodBad: The state directories must be listed in order as st00xx\n";
  }
  @dummy=(sort keys %$akmc_step);
  $totalsteps=@dummy;
  $totalsteps=$totalsteps-1;
  if($totalsteps > 0) {
    $previousst=$akmc_step->{$dummy[$totalsteps-1]}[2];
  }else{
    $previousst="none";
  }
  print "previous st: $previousst\n";
  for($j=0;$j<$numst;$j++) {
    $target=$totalst[$j];
    #print "recovered= $recovered\n";
    ($recovered,$included,$stprs,$st) = LookBackOneSt($target,$stdir,$previousst,$energy,$quench,$recovered,$included,$Equivalency,$Ediffmax,$Rdiffmax,$stfile,$stprs,$st);
    print "recovered= $recovered\n";
  }
  for($i=0;$i<$stpoolsize;$i++) {
    $headpath=$stpool->[$i];
    $filepath=$headpath."/".$stenergyfile;
    $stenergy=ReadStEnergyFile($filepath);
    @totalst=(reverse sort keys %$stenergy);
    $numst=@totalst;
    for($j=0;$j<$numst;$j++){
      $target=$headpath."/".$totalst[$j];
      ($recovered,$included,$stprs,$st) = LookBackOneSt($target,$stdir,$previousst,$energy,$quench,$recovered,$included,$Equivalency,$Ediffmax,$Rdiffmax,$stfile,$stprs,$st);
    }
  }
  return ($recovered,$included,$stprs,$jobs,$st);
}

# ---------------------------------------------------------------------------------------------------------
# Make use of good saddles at the previous step.
# ($recovered,$included,$stprs,$jobs,$st)=UsePreviousGood($stdir,$stdirpre,$curstep,$numsaddles,$SearchesAlgo,$quench,$Rdiffmax,$stfile,$recovered,$included,$stprs,$akmc_step,$jobs,$st)
# ---------------------------------------------------------------------------------------------------------
sub UsePreviousGood{
  my $stdir=shift;
  my $stdirpre=shift;
  my $curstep=shift;
  my $numsaddles=shift;
  my $SearchesAlgo=shift;
  my $quench=shift;
  my $Rdiffmax=shift;
  my $stfile=shift;
  my $recovered=shift;
  my $included=shift;
  my $stprs=shift;
  my $akmc_step=shift;
  my $jobs=shift;
  my $st=shift;
  my $displacement_max=0.25;
  my ($stinfo,$POSCAR,$POSCAR_f,$POSCAR_sp,$numprst,$i,$j,$k,$same,$r_mag, $rdiffmax_sq,$differ,$dummy);
  my ($newprdir,$pr000i,$POSCAR_01,$POSCAR_02, $POSCAR_i,$MODECAR);
  my ($coordinates1,$coordinates2,$coordinates3,$norm,$total_atoms1,$total_atoms2,$total_atoms3,$total_atoms4);
  my ($basis,$lattice,$num_atoms,$selectiveflag,$selective,$description,$filetype);
  my @rij=();
  my $suffocator="";
  $rdiffmax_sq=$Rdiffmax*$Rdiffmax;
  if(lc($SearchesAlgo) eq "dimer") {
    # read two POSCARs and calculate the magnitude of the displacement for each atom
    $POSCAR=$stdirpre."/POSCAR";
    ($coordinates1,$basis,$lattice,$num_atoms,$total_atoms1,$selectiveflag,$selective,$description,$filetype)=read_poscar($POSCAR);
    set_bc($coordinates1,$total_atoms1);
    $POSCAR_f=$stdir."/POSCAR";
    ($coordinates2,$basis,$lattice,$num_atoms,$total_atoms2,$selectiveflag,$selective,$description,$filetype)=read_poscar($POSCAR_f);
    set_bc($coordinates2,$total_atoms2);
    if($total_atoms1 != $total_atoms2) {die "$POSCAR and $POSCAR_f : different number of atoms.\n";}
    $coordinates1=dirkar($coordinates1,$basis,$lattice,$total_atoms1);
    $coordinates2=dirkar($coordinates2,$basis,$lattice,$total_atoms2);
    ($differ,$dummy)=minimum_image_dist($coordinates1,$coordinates2,$basis,$lattice,$total_atoms1);
    for($j=0;$j<$total_atoms1;$j++) {
      push @rij, sqrt($dummy);
    }
    # read in each of good saddles in the previous state and reset it based on the rij between initial and the chosen final state
    $stinfo=ReadSt($stdirpre,$stfile);
    if(!(exists($stinfo->{"numprst"}))) {
      #print "Warning in $stdirpre, $stfile is broken, no numprst. So, set it to be zero.\n";
      $stinfo->{"numprst"}=[(0)];
    }
    # print "numprst in $stdirpre: ".$stinfo->{"numprst"}[0]."\n";
    $numprst=$stinfo->{"numprst"}[0];
    for($i=1;$i<=$numprst;$i++) {
      $pr000i=DirName("pr",$i);
      if((lc($stinfo->{$pr000i}[1]) eq "good") && (lc($stinfo->{$pr000i}[0]) eq "done")) {
        if($included->{$pr000i}) {
          print "Process $pr000i has been ".$included->{$pr000i}.".\n";
          next;
        }
        $recycled=$stdirpre."/".$pr000i;
        $POSCAR_sp=$recycled."/POSCAR_sp";
        $MODECAR=$recycled."/MODECAR";
        ($coordinates3,$basis,$lattice,$num_atoms,$total_atoms3,$selectiveflag,$selective,$description,$filetype)=read_poscar($POSCAR_sp);
        ($norm, $total_atoms4)=read_othercar($MODECAR);
        set_bc($coordinates3,$total_atoms3);
        $coordinates3=dirkar($coordinates3,$basis,$lattice,$total_atoms3);
        if($total_atoms3 != $total_atoms1 || $total_atoms3 != $total_atoms4) {
          die "$POSCAR_sp, $MODECAR and $POSCAR : different number of atoms.\n";
        }
        for($j=0;$j<$total_atoms3;$j++) {
          if($rij[$j] >= $displacement_max) {
            for($k=0;$k<3;$k++){
              $coordinates3->[$j][$k]=$coordinates2->[$j][$k];
              $norm->[$j][$k]=0.0;
            }
          }
        }
        $r_mag=magnitude($norm,$total_atoms3);
        ##print "in UsePreviousGood: r_mag=***$r_mag***\n";
        if($r_mag==0){# all the atoms are shifted, this initial guess is a minimum state
          next;
        }
        # check if this sp search initial config is same as the minimum of the current st state
        ($differ,$dummy)=minimum_image_dist($coordinates3,$coordinates2,$basis,$lattice,$total_atoms3);
        $same=1;
        for($j=0;$j<$total_atoms3;$j++){  
           if($dummy > $rdiffmax_sq) {
             $same=0;
             last;
           }       
        }       
        #print "Same=$same\n";
        if($same) {next;} #this sp initial configuration guess is same as the final state; skip it.
        # prevent atoms from getting too close to each other
        #print "relaxing this initial guess ... \n";
        spring_relaxation($coordinates3,$basis,$lattice,$total_atoms3,$selective);
        #print "... done.\n";
        # ok, recycle this one.
        $recovered++;
        $newprdir=DirName("pr",$recovered);
        $suffocator=`cd $stdir; $Bin/diminit.pl $newprdir 0 0.1 0.1 8`; #just to create the folder pr000x
        $r_mag=1.0/$r_mag;
        $norm=vmult($norm,$r_mag,$total_atoms3); # normalize the recycled mode vector
        $POSCAR_01=$stdir."/".$newprdir."/POSCAR";
        $POSCAR_02=$stdir."/".$newprdir."/MODECAR";
        $coordinates3=kardir($coordinates3,$basis,$lattice,$total_atoms3);
        write_poscar($coordinates3,$basis,$lattice,$num_atoms,$total_atoms3,$selectiveflag,$selective,$description,$POSCAR_01,$filetype); 
        write_othercar($norm, $total_atoms3, $POSCAR_02);
        $included->{$pr000i}="recycled";
        $st->{$newprdir}=[("search","pending","na","na","na")];
      } # end of the if block
    } # end of the for loop
  }elsif(lc($SearchesAlgo) eq "lanczos") {
    die "UsePreviousGood is not ready for Lanczos yet.\n";
  }else{
    die "Unrecognizable saddle point finding algorithm $SearchesAlgo\n";
  }
  return ($recovered,$included,$stprs,$jobs,$st);
}

# ---------------------------------------------------------------------------------------------------------
# Make use of good saddles at the previous step.
# ($recovered,$included,$stprs,$jobs,$st)=UseVectorShift($stdir,$stdirpre,$curstep,$numsaddles,$SearchesAlgo,$quench,$stfile,$recovered,$included,$linker,$stprs,$akmc_step,$jobs,$st)
# ---------------------------------------------------------------------------------------------------------
sub UseVectorShift{
  my $stdir=shift;
  my $stdirpre=shift;
  my $curstep=shift;
  my $numsaddles=shift;
  my $SearchesAlgo=shift;
  my $quench=shift;
  my $stfile=shift;
  my $recovered=shift;
  my $included=shift;
  my $linker=shift;
  my $stprs=shift;
  my $akmc_step=shift;
  my $jobs=shift;
  my $st=shift;
  my ($stinfo,$POSCAR,$POSCAR_f,$POSCAR_sp,$numprst,$i,$j,$k,$r_mag,$dummy);
  my ($newprdir,$pr000i,$POSCAR_01,$POSCAR_02);
  my ($coordinates1,$coordinates2,$coordinates3,$total_atoms1,$total_atoms2,$total_atoms3,$dr);
  my ($basis,$lattice,$num_atoms,$selectiveflag,$selective,$description,$filetype);
  my $suffocator="";
  if(lc($SearchesAlgo) eq "dimer") {
    # read two POSCARs and calculate the magnitude of the displacement for each atom
    $POSCAR=$stdirpre."/POSCAR";
    ($coordinates1,$basis,$lattice,$num_atoms,$total_atoms1,$selectiveflag,$selective,$description,$filetype)=read_poscar($POSCAR);
    set_bc($coordinates1,$total_atoms1);
    $POSCAR_f=$stdir."/POSCAR";
    ($coordinates2,$basis,$lattice,$num_atoms,$total_atoms2,$selectiveflag,$selective,$description,$filetype)=read_poscar($POSCAR_f);
    set_bc($coordinates2,$total_atoms2);
    if($total_atoms1 != $total_atoms2) {die "$POSCAR and $POSCAR_f : different number of atoms.\n";}
    $coordinates1=dirkar($coordinates1,$basis,$lattice,$total_atoms1);
    $coordinates2=dirkar($coordinates2,$basis,$lattice,$total_atoms2);
    # read in each of good saddles in the previous state and reset it based on the rij between initial and the chosen final state
    $stinfo=ReadSt($stdirpre,$stfile);
    if(!(exists($stinfo->{"numprst"}))) {
      #print "Warning in $stdirpre, $stfile is broken, no numprst. Set it to zero.\n";
      $stinfo->{"numprst"}=[(0)];
    }
    $numprst=$stinfo->{"numprst"}[0];
    for($i=1;$i<=$numprst;$i++) {
      $pr000i=DirName("pr",$i);
      if((lc($stinfo->{$pr000i}[1]) eq "good") && (lc($stinfo->{$pr000i}[0]) eq "done")) {
        if($included->{$pr000i}) {
          print "Process $pr000i has been ".$included->{$pr000i}.".\n";
          next;
        }
        $recycled=$stdirpre."/".$pr000i;
        $POSCAR_sp=$recycled."/POSCAR_sp";
        ($coordinates3,$basis,$lattice,$num_atoms,$total_atoms3,$selectiveflag,$selective,$description,$filetype)=read_poscar($POSCAR_sp);
        set_bc($coordinates3,$total_atoms3);      
        if($total_atoms3 != $total_atoms1) {die "in UseVectorShift: $POSCAR_sp and $POSCAR : different number of atoms.\n";}
        $coordinates3=dirkar($coordinates3,$basis,$lattice,$total_atoms3);
        ($dr,$dummy)=minimum_image_dist($coordinates3,$coordinates1,$basis,$lattice,$total_atoms3);
        $coordinates3=vsum($coordinates2,$dr,$total_atoms2);
        $coordinates3=kardir($coordinates3,$basis,$lattice,$total_atoms3);
        set_bc($coordinates3,$total_atoms3);
        $coordinates3=dirkar($coordinates3,$basis,$lattice,$total_atoms3);
        # prevent atoms from getting too close to each other
        #print "relaxing this initial guess ... ...\n";
        spring_relaxation($coordinates3,$basis,$lattice,$total_atoms3,$selective);
        #print "... ... done\n";
        
        # ok, recyle this one.
        $recovered++;
        $newprdir=DirName("pr",$recovered);
        $suffocator=`cd $stdir; $Bin/diminit.pl $newprdir 0 0.1 0.1 8`;
        ($dr,$dummy)=minimum_image_dist($coordinates3,$coordinates2,$basis,$lattice,$total_atoms3);
        $r_mag=0.0;
        foreach $j (@$dummy){
          $r_mag+=$j;
        } 
        $r_mag=sqrt($r_mag);
        #print "in UseVectorShift: r_mag=***$r_mag***\n";
        $r_mag=1.0/$r_mag;
        $dr=vmult($dr,$r_mag,$total_atoms3);
        $POSCAR_01=$stdir."/".$newprdir."/POSCAR";
        $POSCAR_02=$stdir."/".$newprdir."/MODECAR";
        $coordinates3=kardir($coordinates3,$basis,$lattice,$total_atoms3);
        write_poscar($coordinates3,$basis,$lattice,$num_atoms,$total_atoms3,$selectiveflag,$selective,$description,$POSCAR_01,$filetype); 
        write_othercar($dr,$total_atoms3,$POSCAR_02);
        $included->{$pr000i}="recycled";
        $st->{$newprdir}=[("search","pending","na","na","na")];
      }# end of the if block
    } # end of the for loop
  }elsif(lc($SearchesAlgo) eq "lanczos") {
    die "Usevectorshift is not ready for lanczos yet.\n";
  }else{
    die "Unrecognizable sp algorithm $SearchesAlgo\n";
  }
  return ($recovered,$included,$stprs,$jobs,$st);
}

# ---------------------------------------------------------------------------------------------------------
# Check if there is similar atom groups in the previous state
# $stprdir=SurfaceVecShift($curstep,$stpool,$stenergyfile,$Ediffmax,$Rdiffmax,$akmc_step);
# ---------------------------------------------------------------------------------------------------------
sub SurfaceVecShift{
  my ($curstep,$stpool,$stenergyfile,$Ediffmax,$Rdiffmax,$akmc_step)=@_;
  my ($energy,$POSCAR,$POSCAR_j,$stpoolsize,$stenergy,@totalst,$numst,$curst,$i,$j,$similarity,$sim,$stprdir,$linker);
  my $thatisit="";
  my $finallinker="";
  $similarity=0;
  $stprdir=$akmc_step->{$curstep-1}[2]; #There is a previous st at this stage
  $energy=$akmc_step->{$curstep}[3];
  $POSCAR="POSCAR";
  $stpoolsize=@$stpool;
  $stenergy=ReadStEnergyFile($stenergyfile);
  @totalst=(sort keys %$stenergy);
  $numst=@totalst;
  if($numst > 0){
    $curst=$totalst[$numst-1];
    #print $curst. "  ".$numst."\n";
    ($curst eq DirName("st",$numst)) || die "Error in CheckStpool: The state directories must be listed in order as st00xx\n";
  }
  for $j (keys %$stenergy){
    if (abs($energy-$stenergy->{$j}) < $Ediffmax) {
      $POSCAR_j=$j."/POSCAR";
      #print "Evaluate: $POSCAR_j and $POSCAR\n";
      ($thatisit,$sim,$linker)=EvaluateTwoPOSCAR($POSCAR,$POSCAR_j,$Rdiffmax);
      #print "thatisit:***$thatisit***sim:***$sim***\n";
      if($thatisit && ($sim > $similarity)) {
        $stprdir=$j;
        $similarity=$sim;
        $finallinker=$linker;
      }
    }
  }
  for ($i=0;$i<$stpoolsize;$i++){
    $filepath=$stpool->[$i]."/".$stenergyfile;
    $stenergy=ReadStEnergyFile($filepath);
    for $j (keys %$stenergy){
      if (abs($energy-$stenergy->{$j}) < $Ediffmax) {
        $POSCAR_j=$stpool->[$i]."/".$j."/POSCAR";
        ($thatisit,$sim,$linker)=EvaluateTwoPOSCAR($POSCAR,$POSCAR_j,$Rdiffmax);
        if($thatisit && ($sim > $similarity)) {
          $stprdir=$stpool->[$i]."/".$j;
          $similarity=$sim;
          $finallinker=$linker;
        }
      }
    }
  }
  return ($stprdir,$finallinker);
}

# --------------------------------------------------------------------------
# KdbAddPr($curspdir)
# --------------------------------------------------------------------------
sub KdbAddPr{
  my $curspdir=shift;
  my $dummy="";
  $dummy=`$Bin/kdbaddvpr.pl $curspdir`;
}

# ---------------------------------------------------------------------------------------------------------
# $numcon=CreateKdbMatches($stdir, $poscar,$KDBcutoff)
# ---------------------------------------------------------------------------------------------------------
sub CreateKdbMatches{
  my ($stdir, $poscar, $KDBcutoff)=@_;
  my $numcon="";
  my @dummy=();
  chomp($numcon=`cd $stdir; $Bin/kdbquery.pl $poscar $KDBcutoff`);
  @dummy=split(/\s+/, $numcon);
  $numcon=$dummy[-1];
  return $numcon;
}

# ---------------------------------------------------------------------------------------------------------
# ($recovered)=UseKdbMatches($stdir,$numcon,$SearchesAlgo,$recovered);
# ---------------------------------------------------------------------------------------------------------
sub UseKdbMatches{
  my ($stdir,$numcon,$SearchesAlgo,$recovered)=@_;
  my ($i, $prdir, $dummy, $targetSP, $targetMode,@InDIR,@num);
  my $kdbmatches="kdbmatches";
  my $sp_header="SADDLE_";
  my $mode_header="MODE_";
  if(lc($SearchesAlgo) eq "dimer"){
   if(!($numcon=~/^\d+$/)){
     print "Warning: In UseKdbMatches: the number of saddle guesses $numcon is not a digital number.\n";
       $dummy=$stdir."/".$kdbmatches;
       opendir INKDBDIR, $dummy or die "can't open this $dummy dir!";
       @InDir=readdir INKDBDIR;
       @num = grep /^SADDLE_\d+$/, @InDir;
       $numcon=@num;
       @num = grep /^MODE_\d+$/, @InDir;
       $i=@num;
       if($numcon != $i){ $numcon=(($numcon-$i) < 0) ? $numcon : $i; }
       closedir INKDBDIR;
       print $numcon."\n";
   }
   for($i=0; $i<$numcon; $i++){
     $targetSP=$kdbmatches."/".$sp_header.$i;
     $targetMode=$kdbmatches."/".$mode_header.$i;
     if((-e "$stdir/$targetSP") && (-e "$stdir/$targetMode")){
      $recovered++;
      $prdir=DirName("pr", $recovered);
      print "cp $targetSP $prdir/POSCAR;cp $targetMode $prdir/MODECAR\n";
      $dummy=`cd $stdir;mkdir $prdir;cp KPOINTS POTCAR $prdir;cp $targetSP $prdir/POSCAR;cp $targetMode $prdir/MODECAR`;
      #if($dummy eq ""){ system("cd $stdir; rm -r $kdbmatches"); } # do not keep $kdbmatches
     }
   } #end of the for loop
  }elsif(lc($SearchesAlgo) eq "lanczos") {
    die "UseKdbMatches: is not ready for Lanczos yet.\n";
  }else{
    die "Unrecognizable saddle point finding algorithm $SearchesAlgo\n";
  }
  return $recovered;
}

# ---------------------------------------------------------------------------------------------------------
# Setup saddle points searches.  ($jobs,$st)=PreSPSearches($stdir,$curstep,$numsaddles,$SearchesAlgo,$quench,$stpool,
# $PrRecycleAll,$PrRecycle,$UseKDB,$KDBcutoff,$PrRecycleShift,$SurfaceRecShift,$stenergyfile,$stfile,$Equivalency,$Ediffmax,
# $Rdiffmax,$DisplaceAlgo,$DisplaceRange,$NN_rcut,$MaxCoordNum,$akmc_step,$jobs,$st)
# ---------------------------------------------------------------------------------------------------------
sub PreSPSearches{
  my $stdir=shift;
  my $curstep=shift;
  my $numsaddles=shift;
  my $SearchesAlgo=shift;
  my $quench=shift;
  my $stpool=shift;
  my $PrRecycleAll=shift;
  my $PrRecycle=shift;
  my $UseKDB=shift;
  my $KDBcutoff=shift;
  my $PrRecycleShift=shift;
  my $SurfaceRecShift=shift;
  my $stenergyfile=shift;
  my $stfile=shift;
  my $Equivalency=shift;
  my $Ediffmax=shift;
  my $Rdiffmax=shift;
  my $DisplaceAlgo=shift;
  my $DisplaceRange=shift;
  my $NN_rcut=shift;
  my $MaxCoordNum=shift;
  my $akmc_step=shift;
  my $jobs=shift;
  my $st=shift;
  my @dummy=();
  my @stprs=();
  my %included=();
  my ($stdirpre,$prestep,$recovered,$reused,$sp2bcreated,$stprs,$included);
  my ($spprocessed,$i,$prdir,$prstart,$spcreated,$linker,$numcon);
  $recovered=0;
  $reused=0;
  $spprocessed=0;
  $stprs=\@stprs;
  $included=\%included;
  $linker="";
  if($stdir ne $akmc_step->{$curstep}[0]) {die "Error: current st dir $stdir is not same as $akmc_step->{$curstep}[0] in the akmc_step file\n";}
  if(-e "$stdir/INCAR_sp") {system "cd $stdir; cp INCAR_sp INCAR";}
  if((-e "$stdir/INCAR") && (-e "$stdir/KPOINTS") && (-e "$stdir/POTCAR") && (-e "$stdir/POSCAR") && (-e "$stdir/DISPLACECAR_sp")){
    if ($SearchesAlgo == "dimer") {
      if($PrRecycleAll) {# check good and bad saddles to see if the "final" state is same as current initial state in the current step
        print "Preparing saddle point searches by retrieving previous process information.\n";
        $spprocessed=$recovered;
        ($recovered,$included,$stprs,$jobs,$st)=UseGoodBad($stdir,$curstep,$numsaddles,$SearchesAlgo,$quench,$stpool,$stenergyfile,$stfile,$recovered,$included,$Equivalency,$Ediffmax,$Rdiffmax,$stprs,$akmc_step,$jobs,$st);
        $spprocessed=$recovered-$spprocessed;
        $reused+=$spprocessed;
        print "In the current step $curstep, state $stdir, we just retrieved $spprocessed good saddle points from previous step (if any) and the state pool\n";
        print "===========================================================\n";
      }
      if($UseKDB){
        print "Create kdbmatches from kdbdatbase ... ...\n";
        $numcon=CreateKdbMatches($stdir,"POSCAR",$KDBcutoff);
        $spprocessed=$recovered;
        print "Use kdbmatches to prepare saddle point searches... ...\n";
        ($recovered)=UseKdbMatches($stdir,$numcon,$SearchesAlgo,$recovered);
        $spprocessed=$recovered-$spprocessed;
        print "In the current step $curstep, state $stdir, we just prepared $spprocessed saddle point searches from kdbmatches in $stdir\n";
        $reused+=$spprocessed;
      }
      if($PrRecycle) {# make a better initial guess by comparing individual atom shifts between good saddles and the current config
        print "Preparing saddle point searches by recycling good saddles from the previous state (if any).\n";
        sub numerically { $a <=> $b}
        @dummy=(sort numerically keys %$akmc_step);
        $prestep=@dummy;
        $prestep=$prestep-1;
        if($prestep > 0) {
          $stdirpre=$akmc_step->{$dummy[$prestep-1]}[2];
          print "Previous state: $stdirpre\n";
          $spprocessed=$recovered;
          ($recovered,$included,$stprs,$jobs,$st)=UsePreviousGood($stdir,$stdirpre,$curstep,$numsaddles,$SearchesAlgo,$quench,$Rdiffmax,$stfile,$recovered,$included,$stprs,$akmc_step,$jobs,$st);
          $spprocessed=$recovered-$spprocessed;
          $reused+=$spprocessed;
          print "In the current step $curstep, state $stdir, $spprocessed good saddle points were recycled from state $stdirpre.\n";
        }
        print "===========================================================\n";
      }
      if($PrRecycleShift) {# add each vector connecting good saddles and their own initial states to the current initial state
        print "Preparing saddle point searches by vector-shifting good saddles from the selected (default:previous) state (if any).\n";
        sub numerically { $a <=> $b}
        @dummy=(sort numerically keys %$akmc_step);
        $prestep=@dummy;
        $prestep=$prestep-1;
        if($prestep > 0) {
          if($SurfaceRecShift) { # compare configurations and try to select an equivalent or similar state
            ($stdirpre,$linker)=SurfaceVecShift($curstep,$stpool,$stenergyfile,$Ediffmax,$Rdiffmax,$akmc_step);
          }else{
            $stdirpre=$akmc_step->{$dummy[$prestep-1]}[2];
          }
          print "Previous state: $stdirpre\n";
          $spprocessed=$recovered;
          ($recovered,$included,$stprs,$jobs,$st)=UseVectorShift($stdir,$stdirpre,$curstep,$numsaddles,$SearchesAlgo,$quench,$stfile,$recovered,$included,$linker,$stprs,$akmc_step,$jobs,$st);
          $spprocessed=$recovered-$spprocessed;
          $reused+=$spprocessed;
          print "In the current step $curstep, state $stdir, we just vector-shifted $spprocessed good saddle points from state $stdirpre.\n";
        }
        print "===========================================================\n";
      }
      $spcreated=$recovered;
      if($reused < $numsaddles) {
        $sp2bcreated=$numsaddles - $reused;
        $prstart=$recovered+1;
        $spcreated=$recovered+$sp2bcreated;
        for($i=$prstart;$i<=$spcreated;$i++){
          print "Preparing saddle point searches using random displacements\n";
          $prdir=DirName("pr",$i);
          system "cd $stdir; $Bin/diminit.pl $prdir $DisplaceAlgo $DisplaceRange $NN_rcut $MaxCoordNum > /dev/null";
        }
        print "===========================================================\n";
      }
    }elsif ($SearchesAlgo eq "lanczos") {
      #system "cd $stdir; $Bin/laninit.pl $numsaddles > /dev/null";
      die "Saddle point searching script $Bin/laninit.pl is not ready yet.\n";
    }
  }else{
     die "Errors in PreSPSearches: Files in $stdir for setting up saddle point searches are missing\n";
  }
  return ($spcreated,$jobs,$st);
}

# -----------------------------------------------------------
# Check if a stopped SP job is converged [$converged=CheckSP($curspdir,$SearchesAlgo,$energy,$force);]
# -----------------------------------------------------------
sub CheckSP{
  my ($curspdir,$SearchesAlgo,$energy,$force)=@_;
  my $converged=0;
  if(lc($SearchesAlgo) eq "dimer"){
    $converged=CheckMins($curspdir,$energy,$force);
  }elsif(lc($SearchesAlgo)=="lanczos"){
    die "Checking Lanzcos jobs is not ready yet.\n;"
    # $converged=CheckLanczos($curspdir,$energy,$force);
  }else{
    die "Unrecognizable SearchesAlgo \n;"
  }
  return $converged;
}

# ----------------------------------------------------------------------------------
# Get energy and force of the center of the dimer [($energy,$force)=dimef($curspdir)]
# ----------------------------------------------------------------------------------
sub dimef{
  my $curspdir=shift;
  my ($energy,$force);
  my ($i,$j,$gotforce,$lastone,@results,$dummy,@dummy,$outzipped);
  my $zip = $ENV{'VTST_ZIP'} ;
  if($zip eq ''){ $zip = 'gzip' ; }
  $outzipped=0;
  if(!(-e "$curspdir/OUTCAR")){
    if(-e "$curspdir/OUTCAR.gz"){ $outzipped=1; system "cd $curspdir; gunzip OUTCAR.gz";}
    if(-e "$curspdir/OUTCAR.bz2"){ $outzipped=1; system "cd $curspdir; bunzip2 OUTCAR.bz2";}
  }
  $energy=$force="na";
  if(!(-e "$curspdir/OUTCAR")){
    #print "Warning in dimef: there is no OUTCAR in $curspdir.\n";
    return ($energy, $force);
  }
  $dummy=`cd $curspdir; grep -E "energy  without entropy=|Dimer: Central Point|FORCES:" OUTCAR`;
  @results=split /\n/, $dummy;
  #@results=`cd $curspdir; grep -f $Bin/amy OUTCAR`;
  $lastone=@results-1;
  $gotforce=0;
  for($i=$lastone;$i>=0;$i--){ # look backward
    $dummy=$results[$i];
    if($gotforce){
      if($dummy=~/^\s+energy  without entropy/){ # two whitespace before "without entropy"
        $energy=VaspEnergy_OUTCAR($dummy);
        last;
      }
    }elsif($dummy=~/^\s+Dimer: Central Point/){
      $gotforce=1;
      for($j=$i+1;$j<=$lastone;$j++){
        $dummy=$results[$j];
        if($dummy=~/^\s+FORCES\:/){
          $force=VaspForce_OUTCAR($dummy);
          last;
        }
      }
    } # end of the first inner for loop
  }
  if($outzipped){system "$zip $curspdir/OUTCAR";}
  return ($energy,$force);
}

# ----------------------------------------------------------------------------------
sub dimef0{
  my $curspdir=shift;
  my ($energy,$force);
  my ($dummye,@dummye,$dummyf,@dummyf,$outzipped);
  my $zip = $ENV{'VTST_ZIP'} ;
  if($zip eq ''){ $zip = 'gzip' ; }
  $outzipped=0;
  if(-e "$curspdir/OUTCAR.gz") {
      system "cd $curspdir; gunzip OUTCAR.gz";
      $outzipped=1;
  }elsif(-e "$curspdir/OUTCAR.bz2"){
      system "cd $curspdir; bunzip2 OUTCAR.bz2";
      $outzipped=1;
  }
  $energy=$force="na";
  if(!(-e "$curspdir/OUTCAR")){
    #print "Warning in dimef: there is no OUTCAR in $curspdir.\n";
    return ($energy, $force);
  }
  
  $dummye=`cd $curspdir; grep -E "energy without entropy|Dimer: Central Point" OUTCAR | sed -n '/energy without entropy/{N;/\\n.*Dimer: Central Point/p;}' | tail -2 | head -1`;
  $dummye=~s/^\s+//;
  @dummye=split(/\s+/,$dummye);
  $energy=$dummye[7];

  $dummyf=`cd $curspdir; grep -E "FORCES:|Dimer: Central Point" OUTCAR | sed -n '/Dimer: Central Point/{n;p;}' | tail -1`;
  $dummyf=~s/^\s+//;
  @dummyf=split(/\s+/,$dummyf);
  $force=$dummyf[4];

  if($outzipped){system "$zip $curspdir/OUTCAR";}
  return ($energy,$force);
}

# ----------------------------------------------------------------------------------
# Get the curvature for a saddle which is either converged or being running
# ----------------------------------------------------------------------------------
sub GetCurvature{
  my ($curspdir,$SearchesAlgo)=@_;
  my $curvature="na";
  my @line=();
  my $line=();
  if(lc($SearchesAlgo) eq "dimer"){
    if(-e "$curspdir/DIMCAR"){
      $line=`tail -1 $curspdir/DIMCAR`;
      $line=~s/^\s+//;
      if($line ne ""){ # start of "amy" 
        @line=split(/\s+/,$line);
        if($line[4]=~/^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$/){
          $curvature=$line[4];
        }elsif($line[4] eq "---"){
          $line=`tail -2 $curspdir/DIMCAR | head -1`;
          $line=~s/^\s+//;
          if($line ne ""){
            @line=split(/\s+/,$line);
            if($line[4]=~/^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$/){
              $curvature=$line[4];
            }
          }
        }
      } # end of "amy"
    }
  }elsif(lc($SearchesAlgo) eq "lanczos"){
   die "GetCurvature Script for Lanczox is not ready yest.\n";
  }else{
   die "In $curspdir, unrecoganizable $SearchesAlgo.\n";
  }
  return $curvature;
}

# ----------------------------------------------------------------------------------
# Get force and energy from a saddle [($energy,$force,$curvature)=GetSPEF($curspdir,$SearchesAlgo)]
# Only numerical values will be returned. All the other cases are assigned "na";
# ----------------------------------------------------------------------------------
sub GetSPEF{
  my ($curspdir,$SearchesAlgo)=@_;
  my ($energy,$force, $curvature);
  my @line=();
  my $line="";
  if(lc($SearchesAlgo) eq "dimer"){
    if (-e "$curspdir/final"){#cleanedup
      $curvature=GetCurvature($curspdir."/final",$SearchesAlgo);
      ($energy,$force)=dimef($curspdir."/final");
      if($energy eq "na" || $force eq "na" || $curvature eq "na"){
        #print "energy=$energy force=$force curvature=$curvature\n";
        die "Error In $curspdir, the run is supposed to be finished, but we have no energy/force/curvature. Check it.\n";
      }
    }else{ # not cleaned up
      $curvature=GetCurvature($curspdir,$SearchesAlgo);
      ($energy,$force)=dimef($curspdir);
      if($energy eq "na" || $force eq "na" || $curvature eq "na"){
        #print "dimer has no energy and force yet. Please wait a while if it is still running.\n";
      }
    }
  }elsif(lc($SearchesAlgo) eq "lanczos"){
   die "GetSPEF Script for Lanczos is not ready yet.\n";
  }else{
   die "In $curspdir, unrecognizable $SearchesAlgo.\n";
  } 
  #print "energy=$energy force=$force curvature=$curvature \n";
  return ($energy,$force,$curvature);
}

# ----------------------------------------------------
# Start a quench job [MakeSPQuench($curspdir,$SearchesAlgo,$SimR)]
# ----------------------------------------------------
sub MakeSPQuench{
  my ($curspdir,$SearchesAlgo,$SimR)=@_;
  my $idle=0;
  if (-e "$curspdir/mins") {
    print "WARNING: $curspdir/mins already exists. We will not make a new mins folder.Check to make sure it is what you want.\n";
    $idle=1;
    #system "cd $curspdir; rm -r mins";
    return $idle;
  }
  if (-e "$curspdir/final"){#cleanedup
    if(lc($SearchesAlgo) eq "dimer"){
      if(-e "$curspdir/MODECAR") {
        system "cd $curspdir; $Bin/dimmins.pl POSCAR MODECAR $SimR 0 > /dev/null";
      }else{
        die "In $curspdir, there is no MODECAR in the supposedly converged dimer.\n";
      }
    }elsif(lc($SearchesAlgo) eq "lanczos"){
       #system "cd $curspdir; $Bin/lanzmins.pl POSCAR MODECAR $SimR 0 > /dev/null";
       die "Lanczos quench setup script is not ready yet.\n";
    }else{
       die "In $curspdir, unrecognizable SearchesAlgo.\n";
    }
  }else{#not cleaned-up
    if(-e "$curspdir/NEWMODECAR") {
      system "cd $curspdir; $Bin/dimmins.pl CONTCAR NEWMODECAR $SimR 0 > /dev/null";
    }else{
      die "In $curspdir, there is no NEWMODECAR in the supposedly converged dimer.\n";
    }
  }
  return $idle;
}

# -----------------------------------------------------
# Get a saddle configuration  [$POSCAR_sp=GetSP($curspdir,$SearchesAlgo)]
# -----------------------------------------------------
sub GetSP{
  my ($curspdir,$SearchesAlgo)=@_;
  my $POSCAR_sp="";
  if(lc($SearchesAlgo) eq "dimer"){
    if(-e "$curspdir/final"){# cleaned up
      system "cd $curspdir; cp POSCAR POSCAR_sp";
    }elsif(-e "$curspdir/CENTCAR"){
      system "cd $curspdir; cp CENTCAR POSCAR_sp";
    }else{
      die "Error in GetSP: no CENTCAR in a finished run (not cleaned-up) folder\n ";
    }
  }elsif(lc($SearchesAlgo) eq "lanczos"){
    die "Lanczos saddle-grabbing script is not ready yet.\n";
  }else{
    die "In $curspdir, unrecognizable SearchesAlgo.\n";
  }
  return $POSCAR_sp;
}

# ------------------------------------------
# Write jobs [WriteJobs($jobfilename,%jobs)]
# ------------------------------------------
sub WriteJobs{
  my ($jobfilename,$jobs)=@_;
  my $header="";
  my $j;
  $header="joblocation"."  "."JobStatus"."  "."JobId"."  "."JobType"."  "."Energy(eV)"."  "."Force(eV/A)"."  "."Curvature";
  open(JOB,">$jobfilename");
  print JOB $header."\n";
  for $j (sort keys %$jobs){
    print JOB $j."  "."@{$jobs->{$j}}"."\n";
  }
  close JOB; 
}

# ---------------------------------------------------------
# Start a minimization job [($numrunjobs,$jobs)=StartMins($curdir,$numrunjobs,$MaxJobs,$jobs)]
# ---------------------------------------------------------
sub StartMins{
  my ($curdir,$numrunjobs,$MaxJobs,$jobs)=@_;
  my $jobtype="minimization";
  ($numrunjobs,$jobs)=HandleJobSubmission($numrunjobs, $MaxJobs, $curdir, $jobtype, $jobs);
  return ($numrunjobs,$jobs);
}

# -----------------------------------------------------------
# Handle the initialminimization job [($numrunjobs,$jobs)=HandleMins($curdir,$numrunjobs,$MaxJobs,$jobs)]
# -----------------------------------------------------------
sub HandleMins{
   my ($curdir,$numrunjobs,$MaxJobs,$jobs)=@_;
   my ($jobstatus,$jobtype,$converged);
   my ($energy, $force);
   # curdir is the directory that contains min
   if(!(exists($jobs->{$curdir}))){die "There is one $curdir directory without any job information in the jobfile.\n";}
   $jobstatus=lc($jobs->{$curdir}[0]);
   $jobtype=lc($jobs->{$curdir}[2]);
   if($jobstatus eq "running"){
     #print "Job $jobid is still running.\n"
     ($energy,$force)=GetEF($curdir);
     if($energy ne "na"){$jobs->{$curdir}[3]=$energy;}
     if($force ne "na"){$jobs->{$curdir}[4]=$force;}
   }elsif($jobstatus eq "queue" || $jobstatus eq "submittted" ){
     #print "Job in $curdir is still  $jobstatus. Check the queue and make sure it's still alive. We quit.\n";
   }elsif($jobstatus eq "stopped"){ # check if the job is converged.
     ($energy,$force)=GetEF($curdir);
     $converged=CheckMins($curdir,$energy,$force);
     if($energy ne "na"){$jobs->{$curdir}[3]=$energy;}
     if($force ne "na"){$jobs->{$curdir}[4]=$force;}
     if($converged){
       CleanupJobs($curdir,$converged,$jobtype);
       $jobs->{$curdir}[0]="completed";
     }else{
       # check and re-submit the job.
       #if($energy eq "na" || $force eq "na"){
       #   print "job in $curdir never started or at least OUTCAR is not out yet. We'll retry, but better check it!\n";
       #}
       CleanupJobs($curdir,$converged,$jobtype);
       ($numrunjobs,$jobs)=HandleJobSubmission($numrunjobs, $MaxJobs, $curdir, $jobtype, $jobs);
     }
   }elsif($jobstatus eq "2bsubmitted"){
     ($numrunjobs,$jobs)=HandleJobSubmission($numrunjobs, $MaxJobs, $curdir, $jobtype, $jobs);
   }elsif($jobstatus eq "completed"){
     # check and cleanup if necessary
     if(!(-e "$curdir/final")) {CleanupJobs($curdir,1,$jobtype);}
   }elsif($jobstatus eq "killed") {
     #print "Warning: in HandleMins, the job in $curdir has been killed (probably by the user). Make sure that is what you want.\n";
   }else{
     die "Error in $curdir: unrecognizable job status\n";
   }
   return ($numrunjobs,$jobs);
}

# ---------------------------------------------------------------------------------------------------------------
# Handle a dynmat job [($numrunjobs,$jobs,$comment)=HandleDynmat($curpredir,$DynMatFile,$numrunjobs,$MaxJobs,$jobs)]
# ---------------------------------------------------------------------------------------------------------------
 sub HandleDynmat{
   my ($curpredir,$DynMatFile,$numrunjobs,$MaxJobs,$jobs)=@_;
   my ($jobstatuspre,$jobtype);
   my ($converged,@fuss,$dummy,$comment);
   my $maximumnegatives=1;
   $comment="";
   if (-e "$curpredir"){
      if(!(exists($jobs->{$curpredir}))){die "There is one $curpredir without any job information in the jobfile.\n";}
      $jobstatuspre=lc($jobs->{$curpredir}[0]);
      $jobtypepre=lc($jobs->{$curpredir}[2]);
      if(($jobstatuspre eq "running") || ($jobstatuspre eq "queue") || ($jobstatuspre eq "submittted")){
        $dummy=0; # just eat a line and kill the time.
      }elsif($jobstatuspre eq "stopped"){
         $converged=CheckDynmat($curpredir);
         if($converged) {
            @fuss=split(/\//, $curpredir); # st000x/dynmat and st000x/pr000x/dynmat
            if(@fuss == 2){$maximumnegatives=0;}  # dynmat for the initial state
            ($twonegatives,$dynmatfilename)=DynMatrix($curpredir,$DynMatFile,$maximumnegatives);
            if($twonegatives){# bad dynmat e.g., reduce the displacement in the DISPLACECAR or refine the saddle etc.
              CleanupJobs($curpredir,0,$jobtypepre);
              $jobs->{$curpredir}[0]="killed";
              $comment="BadDynmat";
              # fix this rescue steps later on
              #DealwithBadDynMatrix($curpredir);
              #($numrunjobs,$jobs)=HandleJobSubmission($numrunjobs, $MaxJobs, $curpredir, $jobtypepre, $jobs);
              #if($jobs->{$curpredir}[0] eq "killed") {$comment=$jobs->{$curpredir}[1];}
            }else{
              CleanupJobs($curpredir,$converged,$jobtypepre);
              $jobs->{$curpredir}[0]="completed";
            }
         }else{
            # check and re-submit the dynmat job.
            #if (!(-e "$curpredir/OUTCAR")){print "job in $curpredir never started. will retry, but better check it!";} 
            CleanupJobs($curpredir,$converged,$jobtypepre); #It does not copy CONTCAR to POSCAR.
            ($numrunjobs,$jobs)=HandleJobSubmission($numrunjobs, $MaxJobs, $curpredir, $jobtypepre, $jobs);
            if($jobs->{$curpredir}[0] eq "killed") {$comment=$jobs->{$curpredir}[1];}
         }
      }elsif($jobstatuspre eq "completed"){
        #print "It's weird we get to this step. You probably manually set jobstatus for the $curpredir.\n";
        ($twonegatives,$dynmatfilename)=DynMatrix($curpredir,$DynMatFile,$maximumnegatives);
        if($twonegatives){#e.g., reduce the displacement in the DISPLACECAR
          CleanupJobs($curpredir,0,$jobtypepre);
          $jobs->{$curpredir}[0]="killed";
          $comment="BadDynmat";
          # fix this rescue steps later on
          #DealwithBadDynMatrix($curpredir);
          #($numrunjobs,$jobs)=HandleJobSubmission($numrunjobs, $MaxJobs, $curpredir, $jobtypepre, $jobs);
          #if($jobs->{$curpredir}[0] eq "killed") {$comment=$jobs->{$curpredir}[1];}
        }else{
          if (-e "$curpredir/OUTCAR") {CleanupJobs($curpredir, 1, $jobtypepre);}
        }
      }elsif($jobstatuspre eq "2bsubmitted") {
        ($numrunjobs,$jobs)=HandleJobSubmission($numrunjobs, $MaxJobs, $curpredir, $jobtypepre, $jobs);
        if($jobs->{$curpredir}[0] eq "killed") {$comment=$jobs->{$curpredir}[1];}
      }elsif($jobstatuspre eq "killed"){
        #print "Warning from HandleDynmat: the job in $curpredir has been killed (probably by the user). Make sure that is what you want.\n"
      }else{
        die "Unrecognizable job status in $curpredir.\n";
      }
   }else{
      #print "Something weird:  no dynmat folder($curpredir).\n";
      ($numrunjobs,$jobs)=StartDynmat($curpredir,$numrunjobs,$MaxJobs,$jobs);
      if($jobs->{$curpredir}[0] eq "killed") {$comment=$jobs->{$curpredir}[1];}
   }
   return ($numrunjobs, $jobs, $comment);   
}

# -------------------------------------------------------------------------------------------------------
# Start a new st00xx [(%jobs,%st)=StartSt($curdir,$MaxJobs,$NumSearches,$SearchesAlgo,$quench,$dynmat,$turnoff_dynmat,
# $curstep,$stpool,$PrRecycleAll,$PrRecycle,$UseKDB,$KDBcutoff,$PrRecycleShift,$SurfaceRecShift,$stenergyfile,$stfile,$Equivalency,
# $Ediffmax,$Rdiffmax,$DisplaceAlgo,$DisplaceRange,$NN_rcut,$MaxCoordNum,%akmc_step,%jobs,%st)]
# -------------------------------------------------------------------------------------------------------
sub StartSt{
  my $curdir=shift;
  my $MaxJobs=shift; 
  my $NumSearches=shift;
  my $SearchesAlgo=shift;
  my $quench=shift;
  my $dynmat=shift;
  my $turnoff_dynmat=shift;
  my $curstep=shift;
  my $stpool=shift;
  my $PrRecycleAll=shift;
  my $PrRecycle=shift;
  my $UseKDB=shift;
  my $KDBcutoff=shift;
  my $PrRecycleShift=shift;
  my $SurfaceRecShift=shift;
  my $stenergyfile=shift;
  my $stfile=shift;
  my $Equivalency=shift;
  my $Ediffmax=shift;
  my $Rdiffmax=shift;
  my $DisplaceAlgo=shift; 
  my $DisplaceRange=shift; 
  my $NN_rcut=shift;
  my $MaxCoordNum=shift;
  my $akmc_step=shift;
  my $jobs=shift;
  my $st=shift;
  my ($curdyndir,$curprdir,$prdir,$curspdyndir,$curspdir);
  my ($jobtype,$i,$sp2bcreated,$submittedjobs,$spcreated);
  my $curstenergy=0;
  $curstenergy=$akmc_step->{$curstep}[3];
  $submittedjobs=0;
  if(!$turnoff_dynmat) {
    # start the dynmat job for the initial state configuration
    $curdyndir=$curdir."/".$dynmat;
    ($submittedjobs,$jobs)=StartDynmat($curdyndir,$submittedjobs,$MaxJobs,$jobs);
    $st->{"dynmat"}=[()];
    push @{$st->{"dynmat"}},$jobs->{$curdyndir}[0]; # submitted or 2bsubmitted or killed
  }else{
    $st->{"dynmat"}=[()];
    push @{$st->{"dynmat"}},"done";
  }
  # prepare the saddle point searching jobs.
  #$sp2bcreated=(($MaxJobs-$submittedjobs) < $NumSearches) ? $MaxJobs-$submittedjobs : $NumSearches;
  $sp2bcreated=(($MaxJobs-$submittedjobs) > 0) ? $MaxJobs-$submittedjobs : 0;
  ($spcreated,$jobs,$st)=PreSPSearches($curdir,$curstep,$sp2bcreated,$SearchesAlgo,$quench,$stpool,$PrRecycleAll,$PrRecycle,$UseKDB,$KDBcutoff,$PrRecycleShift,$SurfaceRecShift,$stenergyfile,$stfile,$Equivalency,$Ediffmax,$Rdiffmax,$DisplaceAlgo,$DisplaceRange,$NN_rcut,$MaxCoordNum,$akmc_step,$jobs,$st);
  print "Total number of saddle points prepared (including recovered good ones and recycled guesses): $spcreated\n";
  # start the saddle point searching jobs.
  for($i=1;$i<=$spcreated;$i++){
    $prdir=DirName("pr",$i);
    $curspdir=$curdir."/".$prdir;
    $curspdyndir=$curspdir."/".$dynmat;
    if(!(-e $curspdir)) {die "Saddle point search folder $prdir is missing in $curdir\n";}
    if(exists($st->{$prdir})) { # it is either included or recycled
       if(lc($st->{$prdir}[0]) eq "done") { # it was a done state in the recovered st folder w/o dynmat calculation
         if(lc($st->{$prdir}[1]) ne "good") {die "In StartSt: an included $prdir sp is not good. Dead wrong!\n";}
         if($turnoff_dynmat) { # for the current run
           next;
         }else{
           if(-e $curspdyndir) {
             if(-e "$curspdyndir/final") {
               next;
             }else{ # weird, it was marked as done but dynmat was at least not cleaned up. 
               # set the jobs hash to make sure it is to be checked for convergence
               system "cp $curdir/INCAR_* $curspdir; cp $curdir/akmc_dynmat.sub $curspdir;cp $curdir/INCAR_* $curspdyndir; cp $curdir/akmc_dynmat.sub $curspdyndir";
               $jobs->{$curspdyndir}=[("running","jobid_marker","dynmat","na","na")]; # there should be no job called "jobid_marker" in the queue
               $st->{$prdir}[0]="dynmat";
             }
           }else{ # dynmat was not calculated
             system "cp $curdir/INCAR_* $curspdir; cp $curdir/akmc_dynmat.sub $curspdir";
             ($submittedjobs,$jobs)=StartDynmat($curspdyndir,$submittedjobs,$MaxJobs,$jobs);
             $st->{$prdir}[0]="dynmat";
             if(lc($jobs->{$curspdyndir}[0]) eq "killed"){
               $st->{$prdir}[0]="killed";
               $st->{$prdir}[1]=$jobs->{$curspdyndir}[1]."_dynmat";
             }
           }
         }
       }elsif(lc($st->{$prdir}[0]) eq "dynmat") { # it was not a done state in the recovered st folder w/o dynmat calculation
         if(lc($st->{$prdir}[1]) ne "good") {die "In StartSt: an included $prdir sp is not good. Dead wrong! \n";}
         if($turnoff_dynmat) { # for the current run
           $st->{$prdir}[0]="done";
         }else{
           if(-e $curspdyndir) {
             if(-e "$curspdyndir/final") { #weird, this usually means the status must be "done". anyway, it doesn't hurt.
               next;
             }else{
               # set the jobs hash to make sure it is to be checked for convergence
               system "cp $curdir/INCAR_* $curspdir; cp $curdir/akmc_dynmat.sub $curspdir; cp $curdir/INCAR_* $curspdyndir; cp $curdir/akmc_dynmat.sub $curspdyndir";
               $jobs->{$curspdyndir}=[("running","jobid_marker","dynmat","na","na")]; #there should be no job called "jobid_marker" in the queue
               $st->{$prdir}[0]="dynmat";
             }
           }else{ # weird, dynmat was not calculated
             system "cp $curdir/INCAR_* $curspdir; cp $curdir/akmc_dynmat.sub $curspdir";
             ($submittedjobs,$jobs)=StartDynmat($curspdyndir,$submittedjobs,$MaxJobs,$jobs);
             $st->{$prdir}[0]="dynmat";
             if(lc($jobs->{$curspdyndir}[0]) eq "killed"){
               $st->{$prdir}[0]="killed";
               $st->{$prdir}[1]=$jobs->{$curspdyndir}[1]."_dynmat";
             }
           }
         }
       }elsif(lc($st->{$prdir}[0]) eq "search") { # it was recycled
         # copy INCAR_min, INCAR_dynmat and INCAR_sp into pr000x since diminit.pl copy INCAR and akmc.sub for dimer only.
         system "cp $curdir/INCAR_* $curspdir; cp $curdir/akmc_*.sub $curspdir";
         $jobtype=$SearchesAlgo;
         # $jobs->{$curspdir}=[($jobstatus, $jobid, $jobtype, "na", "na", "na")];
         ($submittedjobs,$jobs)=HandleJobSubmission($submittedjobs, $MaxJobs, $curspdir, $jobtype, $jobs);
         $st->{$prdir}=[("search","pending","na","na","na")];
         if($jobs->{$curspdir}[0] eq "killed") {
           $st->{$prdir}[0]="killed";
           $st->{$prdir}[1]=$jobs->{$curspdir}[1]."_search";
         }
       }else{
         die "In StartSt, if St is preset, it must be either done, dynmat, or search.\n";
       }
    }else{ # a random search
       #copy INCAR_min, INCAR_dynmat and INCAR_sp into pr000x since diminit.pl copy INCAR and akmc.sub for dimer only.
       system "cp $curdir/INCAR_* $curspdir; cp $curdir/akmc_*.sub $curspdir";
       $jobtype=$SearchesAlgo;
       # $jobs->{$curspdir}=[($jobstatus, $jobid, $jobtype, "na", "na","na")];
       ($submittedjobs,$jobs)=HandleJobSubmission($submittedjobs, $MaxJobs, $curspdir, $jobtype, $jobs);
       $st->{$prdir}=[("search","pending","na","na","na")];
       if($jobs->{$curspdir}[0] eq "killed") {
         $st->{$prdir}[0]="killed";
         $st->{$prdir}[1]=$jobs->{$curspdir}[1]."_search";
       }
    }
  }
  @{$st->{$curdir}}=($curstenergy);
  @{$st->{"status"}}=("running");
  @{$st->{"NumSearchesLeft"}}=($NumSearches);
  @{$st->{"numprst"}}=($spcreated);
  @{$st->{"prinfolist"}}=("process","status","quality","barrier","final/(min1(eV))","final/min2(eV)");
  return ($jobs,$st);
}

# ----------------------------------------------------------------------------------------------------
# Start a single SP search[($numrunjobs,%jobs,%st)=StartOneSP($curdir,$prdir,$SearchesAlgo,
# $numrunjobs,$MaxJobs,$UseKDB,$DisplaceAlgo,$DisplaceRange,$NN_rcut,$MaxCoordNum,%jobs,%st)] 
# ----------------------------------------------------------------------------------------------------
sub StartOneSP{
  my ($curdir,$prdir,$SearchesAlgo,$numrunjobs,$MaxJobs,$UseKDB,$DisplaceAlgo,$DisplaceRange,$NN_rcut,$MaxCoordNum,$jobs,$st)=@_;
  my $curspdir=$curdir."/".$prdir;
  my $jobtype="";
  if(lc($SearchesAlgo) eq "dimer") {
    system "cd $curdir; $Bin/diminit.pl $prdir $DisplaceAlgo $DisplaceRange $NN_rcut $MaxCoordNum > /dev/null";
    #copy INCAR_min, INCAR_dynmat and INCAR_sp into pr000x since diminit.pl copy INCAR and akmc.sub  for dimer only.
    system "cp $curdir/INCAR_* $curspdir/; cp $curdir/akmc_*.sub $curspdir/";
    $jobtype=lc($SearchesAlgo);
    # $jobs->{$curspdir}=[($jobstatus, $jobid, $jobtype, "na", "na", "na")];
    ($numrunjobs,$jobs)=HandleJobSubmission($numrunjobs, $MaxJobs, $curspdir, $jobtype, $jobs);
    $st->{$prdir}=[("search","pending","na","na","na")];
    if($jobs->{$curspdir}[0] eq "killed") {
      $st->{$prdir}[0]="killed";
      $st->{$prdir}[1]=$jobs->{$curspdir}[1]."_search";
    }
    @{$st->{"status"}}=("running");
  }elsif(lc($SearchesAlgo) eq "lanczos") {
    die "no StartOneSP for lanczoyet.\n";
  }
  return ($numrunjobs,$jobs,$st);
}

# ----------------------------------------------------------------------------------------------------
# Check if a saddle is repeated. [($numsamesaddle,$st)=CheckRepeatedSP($curdir,$prdir,$numprst,$SearchesAlgo,$st)] 
# ----------------------------------------------------------------------------------------------------
sub CheckRepeatedSP{
  my ($curdir,$prdir,$numprst,$SearchesAlgo,$st)=@_;
  my @repeat_prdir=();
  my ($i, $j, $curspdir, $curprdir, $dummy, $numsamesaddle);
  $numsamesaddle=1; # count prdir itself
  $dummy="repeat*".$prdir;
  for($i=1;$i<=$numprst;$i++){
    $curprdir=DirName("pr", $i);
    if(($st->{$curprdir}[1] eq $dummy) && (lc($st->{$curprdir}[0]) eq "done")){
      #the sp repeats a done saddle;
      #$st->{$curprdir}[1]=$dummy."*".$st->{$prdir}[1]; # such as "repeat*pr0003*good" 
      $numsamesaddle++;
    }
  }
  return ($numsamesaddle,$st);
}

# --------------------------------------------------------------------------------------------------
# Start a quench job [($numrunjobs,$jobs,$quenchstatus2,$comment)=
# StartQuench($curspdir,$SearchesAlgo,$numrunjobs,$MaxJobs,$SimR,$jobs)]
# --------------------------------------------------------------------------------------------------
sub StartQuench{
  my ($curspdir,$SearchesAlgo,$numrunjobs,$MaxJobs,$SimR,$jobs)=@_;
  my ($j,$curquendir,$quenchstatus2,$comment,$idle);
  my $quenjobtype="minimization";
  $quenchstatus2=$comment="";
  if(!(exists($jobs->{$curspdir}))) {die "In $curspdir, no job information in the jobs hash.\n"}
  $idle=MakeSPQuench($curspdir, $SearchesAlgo, $SimR); #st0001/pr0003/mins/min1/ min2/ etc.
  if(!(-e "$curspdir/mins")) {die "In $curspdir, the MakeSPQuench did not work successfully. Check it.\n";}
  for($j=1;$j<3;$j++){
    $curquendir=$curspdir."/mins/min".$j;
    if($idle){ # mins have been created before, we register it to the jobs hash, but we would do nothing.
      if(!exists($jobs->{$curquendir})){
        $jobs->{$curquendir}=[("stopped","jobid_marker", $quenjobtype,"na", "na")];
      }
      $comment="nojobssubmitted";
      $quenchstatus2="2bchecked";
      next;
    }
    #cp INCAR_min into quench folders
    system "cp $curspdir/INCAR_min $curquendir/; cp $curspdir/akmc_min.sub $curquendir/";
    ($numrunjobs,$jobs)=HandleJobSubmission($numrunjobs, $MaxJobs, $curquendir, $quenjobtype, $jobs);
    if(lc($jobs->{$curquendir}[0]) eq "killed"){
      $quenchstatus2="killed";
      $comment=$jobs->{$curquendir}[1];
    }
  }
  return ($numrunjobs,$jobs,$quenchstatus2,$comment);
}

# ---------------------------------------------------------------------------------------------------------------------------
# Handle a saddle point search job 
# [($numrunjobs,$jobs,$comment)=HandleSP($curdir,$prdir,$SearchesAlgo,$SPEnergyMax,$numrunjobs,$MaxJobs,$nowait,$jobs)]
# ---------------------------------------------------------------------------------------------------------------------------
sub HandleSP{
  my ($curdir,$prdir,$SearchesAlgo,$SPEnergyMax,$numrunjobs,$MaxJobs,$nowait,$jobs)=@_;
  my ($spstatus,$spjobtype,$converged,$comment);
  my ($energy, $force, $curspdir,$curvature, $POSCAR_sp, $succeed, $dummy);
  $comment="";
  $curspdir=$curdir."/".$prdir;
  if(!(exists($jobs->{$curspdir}))) {die "In $curspdir, no job information in the jobs hash.\n"}
  $spstatus=lc($jobs->{$curspdir}[0]);
  $spjobtype=lc($jobs->{$curspdir}[2]);
  if($spstatus eq "running"){
    ($energy,$force,$curvature)=GetSPEF($curspdir,$SearchesAlgo);
    #print "energy=$energy force=$force curvature=$curvature\n";
    if($energy ne "na"){ $jobs->{$curspdir}[3]=$energy; $UPDATEPROCESS{"energy"}=$energy; }
    if($force ne "na"){ $jobs->{$curspdir}[4]=$force; $UPDATEPROCESS{"force"}=$force; }
    if($curvature ne "na"){ $jobs->{$curspdir}[5]=$curvature; $UPDATEPROCESS{"curvature"}=$curvature; }
    if($energy ne "na"  && $energy > $SPEnergyMax){ # a simpler check
      # it's going crazy or heading to a high energy saddle at best.
      KillJobs($curspdir, $jobs->{$curspdir}[1]);
      $jobs->{$curspdir}[0]="killed";
      $comment="SaddleEnergyTooHigh";
      CleanupJobs($curspdir,0,$SearchesAlgo);
      $numrunjobs--;
    }elsif($energy eq "na" || $force eq "na" || $curvature eq "na"){
      #print "dimer just started. no energy/force/curvature yet. Please wait if it is still running.\n";
      $comment="dimer just started, wait.";
    }
  }elsif($spstatus eq "submitted" || $spstatus eq "queue"){
    $dummy=0; # a trivial line
    if(!$nowait && $spstatus eq "queue"){
      # Halt searches in queue until further notice
      KillJobs($curspdir, $jobs->{$curspdir}[1]);
      $jobs->{$curspdir}[0]="2bsubmitted";
      $comment="wait for quench/dynmat\n";
    }
  }elsif($spstatus eq "stopped"){
    ($energy,$force,$curvature)=GetSPEF($curspdir,$SearchesAlgo);
    $converged=CheckSP($curspdir,$SearchesAlgo,$energy,$force);
    #print "energy=$energy force=$force curvature=$curvature\n";
    if($energy ne "na"){ $jobs->{$curspdir}[3]=$energy; $UPDATEPROCESS{"energy"}=$energy; }
    if($force ne "na"){ $jobs->{$curspdir}[4]=$force; $UPDATEPROCESS{"force"}=$force; }
    if($curvature ne "na"){ $jobs->{$curspdir}[5]=$curvature; $UPDATEPROCESS{"curvature"}=$curvature;}
    if($converged){
      if($curvature < 0.0) {
        #print "$curspdir is converged\n";
        $POSCAR_sp=GetSP($curspdir,$SearchesAlgo);
        $jobs->{$curspdir}[0]="completed";
        CleanupJobs($curspdir,$converged,$SearchesAlgo);
      }else{#bad curvature
        $jobs->{$curspdir}[0]="killed";
        $comment="PositiveCurvature";
        CleanupJobs($curspdir,0,$SearchesAlgo);
      }
    }else{
      if($energy eq "na" || $force eq "na"){
        #print "job in $curspdir could never start. We'll retry, but better check it!\n";
        CleanupJobs($curspdir,$converged,$SearchesAlgo);
        if($nowait){
          ($numrunjobs,$jobs)=HandleJobSubmission($numrunjobs, $MaxJobs, $curspdir, $spjobtype, $jobs);
          if($jobs->{$curspdir}[0] eq "killed") {
             $comment=$jobs->{$curspdir}[1];
          }
        }else{
          $jobs->{$curspdir}[0]="2bsubmitted";
        }
      }else{
        if($energy > $SPEnergyMax){
           #it's heading to a high energy saddle at least; don't resubmit this job
           $jobs->{$curspdir}[0]="killed";
           $comment="SaddleEnergyTooHigh";
           CleanupJobs($curspdir,0,$SearchesAlgo);
        }else{
           #not so bad, so resubmit the job
           CleanupJobs($curspdir,$converged,$SearchesAlgo);
           if($nowait){
            ($numrunjobs,$jobs)=HandleJobSubmission($numrunjobs, $MaxJobs, $curspdir, $spjobtype, $jobs);
            if($jobs->{$curspdir}[0] eq "killed") {$comment=$jobs->{$curspdir}[1];}
           }else{
            $jobs->{$curspdir}[0]="2bsubmitted";
           }
        }
      }
    } # end of if($converged) block 
  }elsif($spstatus eq "completed"){
    # it's weird that a completed and converged search was not followed by a quench. it could be some manual operation.
    # we will clean up the folder and prepare for starting a quench job anyway.
    ($energy,$force,$curvature)=GetSPEF($curspdir,$SearchesAlgo);
    #print "energy=$energy force=$force curvature=$curvature\n";
    if($energy ne "na"){ $jobs->{$curspdir}[3]=$energy; $UPDATEPROCESS{"energy"}=$energy; }
    if($force ne "na"){ $jobs->{$curspdir}[4]=$force; $UPDATEPROCESS{"force"}=$force; }
    if($curvature ne "na"){ $jobs->{$curspdir}[5]=$curvature; $UPDATEPROCESS{"curvature"}=$curvature; 
    }else{
      die "The saddle search in $curspdir must be converged if you MANUALLY marked the job in jobs.dat as 'completed'!\n";
    }
    if($curvature < 0.0) {
      $POSCAR_sp=GetSP($curspdir,$SearchesAlgo);
      $jobs->{$curspdir}[0]="completed";
      CleanupJobs($curspdir,1,$SearchesAlgo);
    }else{#bad curvature
      $jobs->{$curspdir}[0]="killed";
      $comment="PositiveCurvature";
      CleanupJobs($curspdir,0,$SearchesAlgo);        
    }
  }elsif($spstatus eq "2bsubmitted"){
    if($nowait){
      ($numrunjobs,$jobs)=HandleJobSubmission($numrunjobs, $MaxJobs, $curspdir, $spjobtype, $jobs);
      if($jobs->{$curspdir}[0] eq "killed") {$comment=$jobs->{$curspdir}[1];}
    }
  }elsif($spstatus eq "killed"){
    #print "Warning from HandleSP: the job in $curspdir has been killed(probably by the owner). Make sure that is what you want.\n"; 
    if($jobs->{$curspdir}[1]=~/^\d+$/){ # it has a jobid. warning: this digital job id may not be platform-independent
      $comment="UnknownReason";
    }else{ # killed with a reason given by owner or code
      $comment=$jobs->{$curspdir}[1];
    }
  }else{
    die "Unrecognizable jobs status in $curspdir.\n";
  }
  return ($numrunjobs, $jobs,$comment);
} 
     
# ----------------------------------------------------------------------------
# Check if a quench job is stuck at the saddle: 
# ($converged, $tryagain)=DidQuenchDoIt($curspquendir,$j,$SimR,$Ediffmax,$spenergy,$energy,$force)
# ----------------------------------------------------------------------------
sub DidQuenchDoIt{
  my ($curspquendir,$j,$SimR,$Ediffmax,$spenergy,$energy,$force)=@_;
  my $converged=1;
  my $tryagain=0;
  my ($dr,$minj,$curminj,@dummy,$curspdir);
  $minj="min".$j;
  $curminj=$curspquendir."/".$minj;
  if(-e "$curminj/final") {
    # print "It has been cleaned up with the good sign. We take it as well converged.\n";
    return ($converged, $tryagain);
  }
  #chomp($curspdir=`cd $curspquendir; cd ..; pwd`);
  #($spenergy,$force,$curvature)=GetSPEF($curspdir,"dimer");
  #print "spenergy=$spenergy\n";
  if(!($spenergy=~/^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$/)) { # since it is read in from $StFile
    die "In DidQuenchDoIt: No saddle energy/barrier was loaded in; check it.\n";
  }
  if($energy eq "na") { # freshly obtained
    die "Error in DidQuenchDoIt: no energy in $curminj. This shouldn't happen. Something is messed up.\n";
  }
  if(abs($spenergy-$energy) < $Ediffmax) {
    chomp($dr=`cd $curspquendir; dist.pl ../POSCAR_sp $minj/CONTCAR | tail -1`);
    if($dr < $Rdiffmax || (abs($dr-$SimR)<$Rdiffmax && $SimR <= 0.1)) {
      if(!(-e "$curspquendir/$minj/POSCAR_2SimR")) {
        #print "Warning in $curspquendir/$minj, quench didn't do it since the push was too small. We will try to push more.\n";
        $converged=0;
        $tryagain=1;
      #}else{
      # print "Warning in $curspquendir/$minj, quench probably didn't do it and we tried to push more. Have to take it. Check it manually.\n";
      }
    }
  }
  return ($converged,$tryagain);
}

# ------------------------------------------------------------------------------------------------------------------
# Handle a quench job [($numrunjobs,$quenchstatus,$jobs,$comment)=HandleQuench($curspquendir,$SearchesAlgo,$numrunjobs,$MaxJobs,
#           $SimR,$Ediffmax,$Rdiffmax,$spenergy,$jobs)]
# ------------------------------------------------------------------------------------------------------------------
sub HandleQuench{
  my ($curspquendir,$SearchesAlgo,$numrunjobs,$MaxJobs,$SimR,$Ediffmax,$Rdiffmax,$spenergy,$jobs)=@_;
  my ($quenchstatus,$quenjobtype,$quenjobstatus,$converged,$comment);
  my ($energy, $force, $converged, $minj,$minj_short);
  my ($i,$j,$dummy,$tryagain,$TwoSimR);
  my $deadmark=0;
  $comment="";
  if(lc($SearchesAlgo) eq "dimer"){
    $i=0;
    for($j=1;$j<3;$j++){# two images
      $minj_short="min".$j;
      $minj=$curspquendir."/".$minj_short;
      if(!(exists($jobs->{$minj}))) {die "In $minj, no job information in the jobs hash.\n";}
      $quenjobstatus=lc($jobs->{$minj}[0]);
      $quenjobtype=lc($jobs->{$minj}[2]);
      if($quenjobstatus eq "running"){
        ($energy,$force)=GetEF($minj);
        if($energy ne "na") { $jobs->{$minj}[3]=$energy; $UPDATEPROCESS{"energy"}=$energy;}
        if($force ne "na") { $jobs->{$minj}[4]=$force; $UPDATEPROCESS{"force"}=$force;}
        $UPDATEPROCESS{"comment"}=$UPDATEPROCESS{"comment"}."$minj_short(running)";
      }elsif($quenjobstatus eq "submitted" || $quenjobstatus eq "queue"){
        $UPDATEPROCESS{"comment"}=$UPDATEPROCESS{"comment"}."$minj_short(queue)";
      }elsif($quenjobstatus eq "stopped"){
        $tryagain=0;
        ($energy,$force)=GetEF($minj);
        $converged=CheckMins($minj,$energy,$force);
        if($energy ne "na") { $jobs->{$minj}[3]=$energy;}
        if($force ne "na") { $jobs->{$minj}[4]=$force;}
        if($converged){($converged,$tryagain)=DidQuenchDoIt($curspquendir,$j,$SimR,$Ediffmax,$spenergy,$energy,$force);}
        if($converged){
          CleanupJobs($minj,$converged,$quenjobtype);
          $jobs->{$minj}[0]="completed";
          $i++;
        }else{
          if($energy ne "na"){$UPDATEPROCESS{"energy"}=$energy;}
          if($force ne "na") {$UPDATEPROCESS{"force"}=$force;}
          CleanupJobs($minj,$converged,$quenjobtype);
          if($tryagain) {
            $TwoSimR=2.0*$SimR;
            system "cd $curspquendir/..; $Bin/dimmins.pl POSCAR MODECAR $TwoSimR $minj_short;";
            system "cd $minj; cp POSCAR POSCAR_2SimR;";
            $UPDATEPROCESS{"comment"}=$UPDATEPROCESS{"comment"}."$minj_short(re-quenched)";
          }else{$UPDATEPROCESS{"comment"}=$UPDATEPROCESS{"comment"}."$minj_short(continued)";}
          ($numrunjobs,$jobs)=HandleJobSubmission($numrunjobs, $MaxJobs, $minj, $quenjobtype, $jobs);
          if($jobs->{$minj}[0] eq "killed") {
            $comment=$jobs->{$minj}[1];
            $UPDATEPROCESS{"comment"}=$UPDATEPROCESS{"comment"}."(killed_$comment)";
          }
        }  
      }elsif($quenjobstatus eq "completed"){
        # if it's not cleaned up, we will clean up the folder and start a quench job anyway.
        $converged=1;
        ($energy,$force)=GetEF($minj);
        if($energy ne "na") { $jobs->{$minj}[3]=$energy;}
        if($force ne "na") { $jobs->{$minj}[4]=$force;}
        if($converged){($converged,$tryagain)=DidQuenchDoIt($curspquendir,$j,$SimR,$Ediffmax,$spenergy,$energy,$force);}
        if($converged){
          $i++;
          CleanupJobs($minj,$converged,$quenjobtype);
        }else{
          if($energy ne "na"){$UPDATEPROCESS{"energy"}=$energy;}
          if($force ne "na") {$UPDATEPROCESS{"force"}=$force;}
          CleanupJobs($minj,$converged,$quenjobtype);
          if($tryagain) {
            $TwoSimR=2.0*$SimR;
            system "cd $curspquendir/..; $Bin/dimmins.pl POSCAR MODECAR $TwoSimR $minj_short;";
            system "cd $minj; cp POSCAR POSCAR_2SimR;";
            $UPDATEPROCESS{"comment"}=$UPDATEPROCESS{"comment"}."$minj_short(re-quenched)";
          }else{$UPDATEPROCESS{"comment"}=$UPDATEPROCESS{"comment"}."$minj_short(continued)";}
          ($numrunjobs,$jobs)=HandleJobSubmission($numrunjobs, $MaxJobs, $minj, $quenjobtype, $jobs);
          if($jobs->{$minj}[0] eq "killed") {
            $comment=$jobs->{$minj}[1];
            $UPDATEPROCESS{"comment"}=$UPDATEPROCESS{"comment"}."(killed_$comment)";
          }
        }
      }elsif($quenjobstatus eq "2bsubmitted"){
        ($numrunjobs,$jobs)=HandleJobSubmission($numrunjobs, $MaxJobs, $minj, $quenjobtype, $jobs);
        if($jobs->{$minj}[0] eq "killed") {
          $comment=$jobs->{$minj}[1];
          $UPDATEPROCESS{"comment"}=$UPDATEPROCESS{"comment"}."$minj_short(killed_$comment)";
        }
      }elsif($quenjobstatus eq "killed"){
        $deadmark=1;
        $comment=$minj_short."Failed";
        #print "Warning from HandleQuench: the job in $minj has been killed (probably by the user). Make sure that is what you want.\n";
      }else{
        die "Unrecognizable jobs status in $minj.\n";
      }
    } # end of the "for" loop.
  }elsif(lc($SearchesAlgo) eq "lanczos"){
    die "Error in HanndleQuench, LANCZOS is not integrated in this akmc code yet.\n";
  }else{
    die "Error in HanddleQuench, unrecognizable saddle point algorithm.\n";
  }
  if ($i == 2) {
    $quenchstatus="completed";
    if(abs($jobs->{$curspquendir."/min1"}[3]-$jobs->{$curspquendir."/min2"}[3]) < $Ediffmax){
      $dummy=CompareTwoPOSCAR($curspquendir."/min1/POSCAR",$curspquendir."/min2/POSCAR",$Rdiffmax);
      if($dummy){ # trouble: both ends are same;
        $j=$curspquendir."_1sttrial"; 
        if(-e $j){ # we tried but it still ends up same.
          $quenchstatus="killed";
          $comment="min1=min2";
        }else{
          $quenchstatus="pushmore";
          system("mv $curspquendir $j");
          $comment="failed: same endpoints";
        }
      }
    }
    $UPDATEPROCESS{"energy"}=$jobs->{$minj}[3];
    $UPDATEPROCESS{"force"}=$jobs->{$minj}[4];
  }else{
    $quenchstatus="running";
  }
  if($deadmark) {
    $quenchstatus="killed";
    $UPDATEPROCESS{"comment"}="QuenchDoomed:$comment";
  }
  return ($numrunjobs,$quenchstatus,$jobs,$comment);
}

# -----------------------------------------------------------------------
# Evaluate saddle points [%st=EvalSP($curdir,$prdir,$quench,$SearchesAlgo,$minsinfo,$Ediffmax,$Rdiffmax,$st)]
# -----------------------------------------------------------------------
sub EvalSP{
  my ($curdir,$prdir,$quench,$SearchesAlgo,$minsinfo,$Ediffmax,$Rdiffmax,$st)=@_;
  my $minj="";
  my ($i,$j,$minj,$minj_short);
  my $initial="";
  my ($POSCARPath, $minenergy, $minPOSCAR,$final);
  my ($final_energy,$min1energy,$min2energy);
  $energy=$st->{$curdir}[0];
  $POSCARPath="$curdir/POSCAR";
  if(lc($SearchesAlgo) eq "dimer"){
    $initial="";
    $final="";
    for($j=1;$j<3;$j++){
      $minj_short=$quench."/min".$j;
      $minj=$curdir."/".$prdir."/".$minj_short;
      $minPOSCAR=$minj."/POSCAR";
      if(exists($minsinfo->{$minj_short})) {
        $minenergy=$minsinfo->{$minj_short};
      }else{
        die "In EvalSP of $curdir: No mins information passed in.\n";
      }
      if($j == 1) {
         $min1energy=$minenergy;
      }else{
         $min2energy=$minenergy;
      }
      if(abs($energy - $minenergy) < $Ediffmax){
         $same=CompareTwoPOSCAR($POSCARPath,$minPOSCAR,$Rdiffmax); # compare two POSCAR files.
         if($same) {
           if(!$initial){
             $initial=$minj;
           }else{
             # both minima end up on the same initial minimum, but it should have been excluded.
             print "Quench in $prdir lead to same minima. This should not happened. Go check it.\n";
             $final=$minj_short;
             $final_energy=$minenergy;
           }
         }else{
           $final=$minj_short;
           $final_energy=$minenergy;
         }
      }else{
         $final=$minj_short;
         $final_energy=$minenergy;
      }      
    } # end of the "for" loop.
    if($initial) {
      $st->{$prdir}[1]="good";
      $st->{$prdir}[3]=$final;
      $st->{$prdir}[4]=$final_energy;
    }else{
      $st->{$prdir}[1]="bad";
      $st->{$prdir}[3]=$min1energy;
      $st->{$prdir}[4]=$min2energy;
    }
  }elsif(lc($SearchesAlgo) eq "lanczos"){
    die "Error in HanndleQuench, LANCZOS is not integrated in this akmc code yet.\n";
  }else{
    die "Error in HanddleQuench, unrecognizable saddle point algorithm.\n";
  }
  return $st;
}

# -------------------------------------------------------------------------------------
# Start a postquench dynmat job [($numrunjobs,$jobs)=StartDynmat($curspdyndir,$numrunjobs,$MaxJobs,$jobs);]
# -------------------------------------------------------------------------------------
sub StartDynmat{
  my ($curspdyndir,$numrunjobs,$MaxJobs,$jobs)=@_;
  my $jobtype="dynmat";
  my $curspdir="";
  my $idle;
  if(-e "$curspdyndir"){ # dynmat already exists, do nothing.
    if(!exists($jobs->{$curspdyndir})){
     $jobs->{$curspdyndir}=[("stopped", "5xg42hd52d2v423fjgf",$jobtype,"na", "na")];
    }
    return ($numrunjobs,$jobs);
  }
  MakeNewSt($curspdyndir);
  chomp($curspdir=`cd $curspdyndir; cd ..; pwd`);
  #print "$curspdir\n";
  if(-e "$curspdir/POSCAR_sp"){ # this is a dynmat for a process.
    system "cp $curspdir/POSCAR_sp $curspdyndir/POSCAR";
  }elsif(-e "$curspdir/POSCAR"){
    system "cp $curspdir/POSCAR $curspdyndir/POSCAR";
  }
  ($numrunjobs,$jobs)=HandleJobSubmission($numrunjobs, $MaxJobs, $curspdyndir, $jobtype, $jobs);
  return ($numrunjobs,$jobs);
}


# -------------------------------------------------------------------------------------
# Evalaute a DONE process again:$st=ReEvaluateProcess($curdir,$prdir,$quench,$SearchesAlgo,$Ediffmax,$Rdiffmax,$jobs,$st)  # LX: this subroutine is not being used.
# -------------------------------------------------------------------------------------
sub ReEvaluateProcess{
  my ($curdir,$prdir,$quench,$SearchesAlgo,$turnoff_dynmat,$Ediffmax,$Rdiffmax,$jobs,$st)=@_;
  my ($curspquendir,$curspdir);
  my ($reassessed_succeed,$energy_re,$force_re,$curvature_re,%minsinfo);
  $curspdir=$curdir."/".$prdir;
  $curspquendir=$curspdir."/".$quench;
  #print "$prdir is to be re-evaluated: update the barrier, etc information.\n";
  #get the barrier again
  $reassessed_succeed=0;
  if(exists($jobs->{$curspdir}) && $jobs->{$curspdir}[3]=~/^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$/){
    $st->{$prdir}[2]=$jobs->{$curspdir}[3];
  }else{
    ($st->{$prdir}[2],$force_re,$curvature_re)=GetSPEF($curspdir,$SearchesAlgo);
  }
  if($st->{$prdir}[2]=~/^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$/){
    #get the quench information again.
    if(exists($jobs->{$curspquendir."/min1"}) && exists($jobs->{$curspquendir."/min2"})){
      %minsinfo=("$quench/min1"=>$jobs->{$curspquendir."/min1"}[3],"$quench/min2"=>$jobs->{$curspquendir."/min2"}[3]);
      $reassessed_succeed=1;
    }elsif(exists($jobs->{$curspquendir."/min1"})){
      ($energy_re,$force_re)=GetEF($curspquendir."/min2");
      if($energy_re=~/^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$/){
        %minsinfo=("$quench/min1"=>$jobs->{$curspquendir."/min1"}[3], "$quench/min2"=>$energy_re);
        $reassessed_succeed=1;
      }else{
        print "Re-evaluation of $quench/min1 shows no digital barriers. We skip the rest of re-evaluation.\n";
      }
    }elsif(exists($jobs->{$curspquendir."/min2"})){
        ($energy_re,$force_re)=GetEF($curspquendir."/min1");
        if($energy_re=~/^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$/){
          %minsinfo=("$quench/min2"=>$jobs->{$curspquendir."/min2"}[3], "$quench/min1"=>$energy_re);
          $reassessed_succeed=1;
        }else{
          print "Re-evaluation of $quench/min2 shows no digital barriers. We skip the rest of re-evaluation.\n";
        }
    }else{
        ($energy_re,$force_re)=GetEF($curspquendir."/min1");
        if($energy_re=~/^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$/){
          %minsinfo=("$quench/min1"=>$energy_re);
          ($energy_re,$force_re)=GetEF($curspquendir."/min2");
          if($energy_re=~/^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$/){
            $minsinfo{"$quench/min2"}=$energy_re;
            $reassessed_succeed=1;
          }else{
            print "Re-evaluation of $quench/min2 shows no digital barriers. We skip the rest of re-evaluation.\n";
          }
        }else{
          print "Re-evaluation of $quench/min1 shows no digital barriers. We skip the rest of re-evaluation.\n";
        }
    }
  }else{
    print "Re-evaluation of $prdir shows no digital barriers. We skip the rest of re-evaluation.\n";
  }
  if($reassessed_succeed){
    $st=EvalSP($curdir,$prdir,$quench,$SearchesAlgo,\%minsinfo,$Ediffmax,$Rdiffmax,$st);
    if($st->{$prdir}[1] eq "good"){
      if($turnoff_dynmat){$st->{$prdir}[0]="done";
      }else{$st->{$prdir}[0]="dynmat";}
    }else{
      $st->{$prdir}[0]="done";
    }
  }
  return $st;
}

# ------------------------------------------------------------------------------------------------------
# Filter saddle points after convergence 
# [%st=FilterSP($curdir,$prdir,$numprst,$BarrierMax,$SearchesAlgo,$quench,$dynmat,$stfile,$turnoff_dynmat,
#         $Equivalency,$stpool,$GrabQuenched,$Ediffmax,$Rdiffmax,$DynMatFile,$stdirs,$st)]
# ------------------------------------------------------------------------------------------------------
sub FilterSP{
  my $curdir=shift;
  my $prdir=shift;
  my $numprst=shift;
  my $BarrierMax=shift;
  my $SearchesAlgo=shift;
  my $quench=shift;
  my $dynmat=shift;
  my $stfile=shift;
  my $turnoff_dynmat=shift;
  my $Equivalency=shift;
  my $stpool=shift;
  my $GrabQuenched=shift;
  my $Ediffmax=shift;
  my $Rdiffmax=shift;
  my $DynMatFile=shift;
  my $stdirs=shift;
  my $st=shift;
  my $curspdir=$curdir."/".$prdir;
  my $barrier=0.0;
  my $j=0;
  my ($target, $spenergy,$status, $quality,$same, $curpsdir1, $curpsdir2,$whereabout);
  my $recycled=0;
  $spenergy=$st->{$prdir}[2];
  if($spenergy eq "na" || $spenergy eq "") {die "In $curspdir, the so-called converged saddle point run has no force or energy info.\n";}
  $barrier=$spenergy - $st->{$curdir}[0];
  if(($barrier >0.0) && ($barrier <= $BarrierMax)){
    $st->{$prdir}[1]="promising";
    for($j=1;$j<=$numprst;$j++){
      $target=DirName("pr", $j);
      if($target ne $prdir){
        #just compare it with anyone not marked as repeated, no matter whether the target saddle be good, high energy or bad,
        $status=lc($st->{$target}[0]);
        $quality=lc($st->{$target}[1]);
        if($quality eq "good" || $quality eq "bad" || $quality eq "promising" || $quality eq "highenergy"){
           $curpsdir1=$curdir."/".$target."/POSCAR_sp";
           $curpsdir2=$curspdir."/POSCAR_sp";
           if(abs($spenergy - $st->{$target}[2]) < $Ediffmax){
              $same=CompareTwoPOSCAR($curpsdir1,$curpsdir2,$Rdiffmax); # compare two POSCAR_sp files.
              if($same) {
                $st->{$prdir}[1]="repeat"."*".$target;
                last;
              }
           }
        }
      }
    } # end of the "for" loop
    if($st->{$prdir}[1] eq "promising") { # found no repeat in the current state folder, how about others?
      $j=!((-e "$curspdir/$quench/min1") && (-e "$curspdir/$quench/min2")); # make sure no qunech folders yet.
      if($GrabQuenched && $j) { # check if the saddle has been done in previous states
        ($recycled,$whereabout,$st)=GrabQuenchedSP($curdir,$prdir,$SearchesAlgo,$quench,$dynmat,$stfile,$turnoff_dynmat,$Equivalency,$stpool,$Ediffmax,$Rdiffmax,$DynMatFile,$stfile,$stdirs,$st); # after finishing debug, we will remove variable $whereabout.
        #if($recycled) {print "The saddle had been quenched at $whereabout, ";}
      }
    }
  }elsif($barrier > $BarrierMax){
    $st->{$prdir}[1]="highenergy";
  }else{
    #print "In $curspdir, a negative barrier of $barrier. mark it as bad\n";
    $st->{$prdir}[1]="bad";
  }
  return ($recycled,$st);
}

# -----------------------------------------------------------------------------
# Calculate the prefactor [$prefactor=Calc_Prefactor($dynmat_init,$dynmat_sp)]
# -----------------------------------------------------------------------------
sub Calc_Prefactor{
  my ($dynmat_init,$dynmat_sp)=@_;
  my $prefactor=0.0;
  my $line="";
  my @line=();
  my $last=100;
  if(!(-e $dynmat_init) || !(-e $dynmat_sp)) {
    print "Warning : $dynmat_init or $dynmat_sp does not exist. We will assume a prefactor of 1e12 1/s.\n";
    $prefactor=1.0e12;
    return $prefactor;
  }
  chomp($line=`$Bin/dymprefactor.pl $dynmat_init $dynmat_sp | tail -1;`);
  $line=~s/^\s+//;
  @line=split(/\s+/,$line);
  $last=@line;
  if($line[$last-1] eq "frequency") {
    print "Warning of a failed prefactor calculation: \n";
    print "@line"."\n";
    print "We will replace it with the standard prefactor of 1e12 1/s.\n";
    $prefactor=1.0e12;
  }else{
    $prefactor=$line[$last-1]*1.0e12 ; # Terraherz
  }
  return $prefactor;
}

# ---------------------------------------------------------------------------
# Calculate the rate constant [$rate=Calc_Rate($barrier, $prefactor, $Temperature)]
# ---------------------------------------------------------------------------
sub Calc_Rate{
  my ($barrier,$prefactor,$temperature)=@_;
  my $rateconst=0.0;
  if ($temperature==0) {die "In Calc_Rate: temperature is zero.\n";}
  $rateconst=$prefactor*exp(-1.0*$barrier*11604.5/$temperature);
  return $rateconst;
}

# ---------------------------------------------------------------------------------------------------------
# Categorize saddle points ($unique_sp,$total_sp,$good_saddles_inprocess,$degeneracy,$energy)=CountSP($curdir,$SearchesAlgo,$NumKT,$Temperature,$st)
# ---------------------------------------------------------------------------------------------------------
sub CountSP{
  my ($curdir,$SearchesAlgo,$NumKT,$Temperature,$st)=@_;
  my $MinSP="";
  my %energy=();
  my %degeneracy=();
  my @dummy=();
  my ($unique_sp,$total_sp,$good_saddles_inprocess,$belowTenKT, $belowTenKT_total,$i,$prdir,$numprst,$numsamesaddle,$tenKT,$spmin,$spminKT);
  $unique_sp=$total_sp=$good_saddles_inprocess=$belowTenKT=$belowTenKT_total=0;
  $numprst=$st->{"numprst"}[0];
  for($i=1;$i<=$numprst;$i++) {
    $prdir=DirName("pr",$i);
    # analyze good saddle points
    if(lc($st->{$prdir}[1]) eq "good") {
      $good_saddles_inprocess++; # including done good saddles and those in dynmat steps
      if(lc($st->{$prdir}[0]) eq "done") {
        ($numsamesaddle,$st)=CheckRepeatedSP($curdir,$prdir,$numprst,$SearchesAlgo,$st);
        $degeneracy{$prdir}=[($numsamesaddle)];
        $energy{$prdir}=$st->{$prdir}[2];
        $unique_sp++;
        $total_sp=$total_sp + $numsamesaddle;
      }
    }
  }
  @dummy=(keys %energy); # both are at the same order
  if(@dummy== 0) {
    #print "no good saddle has been found yet.\n";
    return ($unique_sp, $total_sp, $good_saddles_inprocess, $belowTenKT,$belowTenKT_total,$MinSP,\%degeneracy, \%energy);
  }else{
    $MinSP=$dummy[0];
    $spmin=$energy{$MinSP};
    for($i=1;$i<@dummy;$i++) {
      $prdir=$dummy[$i];
      if($spmin > $energy{$prdir}) {
        $MinSP=$prdir;
        $spmin=$energy{$MinSP};
      }
    }
    $tenKT=$NumKT*$Temperature/11604.5; # e.g., eV for 10KT
    $spminKT=$spmin+$tenKT;
    for($i=0;$i<@dummy;$i++) {
      $prdir=$dummy[$i];
      if($energy{$prdir} <= $spminKT) {
        $belowTenKT++;
        $belowTenKT_total=$belowTenKT_total+$degeneracy{$prdir}[0];
      }
    }  
    return ($unique_sp,$total_sp,$good_saddles_inprocess,$belowTenKT,$belowTenKT_total,$MinSP,\%degeneracy,\%energy);
  }
}  

# -----------------------------------------------------------------------------------------------------------------
# $goodandsound=EvaluateLowEnergySPPopulation($curdir,$NumSearches,$SearchesAlgo,$NumKT,$Temperature,$population,$st)
# -----------------------------------------------------------------------------------------------------------------
sub EvaluateLowEnergySPPopulation{
  my ($curdir,$NumSearches,$SearchesAlgo,$NumKT,$Temperature,$population,$st)=@_;
  my $goodandsound=0;
  my $prdir="";
  my $en="";
  my ($minimum,$percent_pop,$numprst);
  my @ener=();
  my ($unique_sp,$total_sp,$good_saddles_inprocess,$belowTenKT,$belowTenKT_total,$MinSP,$degeneracy,$energy);
  $unique_sp=$total_sp=$good_saddles_inprocess=$belowTenKT=0;
  $numprst=$st->{"numprst"}[0];
  print "counting good saddle points that have already been found ... ...\n";
  ($unique_sp,$total_sp,$good_saddles_inprocess,$belowTenKT,$belowTenKT_total,$MinSP,$degeneracy,$energy)=CountSP($curdir,$SearchesAlgo,$NumKT,$Temperature,$st);
  print "$unique_sp unique saddle points among total $total_sp good saddle points found.\n";
  if(exists($degeneracy{""})) {delete $degeneracy{""};}
  if($unique_sp >= $NumSearches) {
   if($population <= 0) {
     $goodandsound=1;
     return $goodandsound;
   }else{
     for $prdir (sort keys %$degeneracy) {
       push @{$degeneracy->{$prdir}}, $degeneracy->{$prdir}[0]/$total_sp ;
     }
   
     @ener=(sort {$energy->{$a} <=> $energy->{$b}} keys %$energy); 
     $minimum=$energy->{$ener[0]};
     $percent_pop=0.0;
     for $en (@ener) {
       if (($energy->{$en} - $minimum) < 0.01) {
          # equivalent low energy saddle points probably due to symmetry
          $percent_pop=$percent_pop + $degeneracy->{$en}[1];
       }
       # how do we do the population stuff? now, 20% for the lowest energy saddle
     }
     if($percent_pop >= $population) {
       print "Good low energy saddle population: target: $population and reality: $percent_pop \n";
       $goodandsound=1;
       return $goodandsound;
     }
   } # end of if population
  } # end of the first if
  # now, either we don't have enough good saddles or a good low energy saddle population 
  return $goodandsound;
}

# -----------------------------------------------------------------------------------------------------------------
# ($NumSearchesLeft,$goodandsound)=EvaluateLowEnergySPTenKT($curdir,$NumSearches,$SearchesAlgo,$NumKT,$Temperature,$population,$st_old,$st)
# -----------------------------------------------------------------------------------------------------------------
sub EvaluateLowEnergySPTenKT{
  my ($curdir,$NumSearches,$SearchesAlgo,$NumKT,$Temperature,$population,$st_old,$st)=@_;
  my ($NumSearchesLeft,$numprst);
  my $goodandsound=0;
  my ($unique_sp,$total_sp,$numsamesaddle,$good_saddles_inprocess,$belowTenKT,$belowTenKT_total,$MinSP,$degeneracy,$energy);
  my ($unique_sp0,$total_sp0,$numsamesaddle0,$good_saddles_inprocess0,$belowTenKT,$belowTenKT_total0,$MinSP0,$degeneracy0,$energy0);
  my $totalgoodbelowtenKT=0;
  # count how many unique good saddle points were found before the update of the st information.
   print "===========================================================\n";
   print "Counting good saddle points already found before updating the state.\n";
   ($unique_sp0,$total_sp0,$good_saddles_inprocess0,$belowTenKT0,$belowTenKT_total0,$MinSP0,$degeneracy0,$energy0)=CountSP($curdir,$SearchesAlgo,$NumKT,$Temperature,$st_old);
   print "$unique_sp0 unique saddle points among total $total_sp0 good saddle points found.\n";
   print "$belowTenKT0 unique saddle points below $NumKT KT+SPMIN among total $belowTenKT_total0 good saddle points below $NumKT KT+SPMIN\n";
   if(exists($degeneracy0{""})) {delete $degeneracy0{""};}
  # count how many unique good saddle points have been found.
   print "===========================================================\n";
   print "Counting good saddle points already found.\n";
   ($unique_sp,$total_sp,$good_saddles_inprocess,$belowTenKT,$belowTenKT_total,$MinSP,$degeneracy,$energy)=CountSP($curdir,$SearchesAlgo,$NumKT,$Temperature,$st);
   print "$unique_sp unique saddle points among total $total_sp good saddle points found.\n";
   print "$belowTenKT unique saddle points below $NumKT KT+SPMIN among total $belowTenKT_total good saddle points below $NumKT KT+SPMIN\n";
   if(exists($degeneracy{""})) {delete $degeneracy{""};}
   $NumSearchesLeft=$st->{"NumSearchesLeft"}[0];
   $numprst=$st->{"numprst"}[0];
   if($total_sp0 < $total_sp) { # found more good saddle points
     if($total_sp0 == 0) { # just found the first spmin       
       $NumSearchesLeft=$NumSearches;
     }else { # in the middle of the run
       print "The lowest saddle: $MinSP0 (before) and $MinSP (now)\n";
       if($MinSP0 eq $MinSP) { # the lowest energy saddle is same
          print "Unique saddles below $NumKT KT+SPMIN: $belowTenKT0 (before) and $belowTenKT (now)\n";
         if($belowTenKT0 == $belowTenKT) { # found no new saddle between spmin and spmin+nKT
           $totalgoodbelowtenKT=$belowTenKT_total-$belowTenKT_total0;
           if($totalgoodbelowtenKT >= $NumSearchesLeft) {
             $NumSearchesLeft=0;
           }elsif($totalgoodbelowtenKT >= 0) {
             $NumSearchesLeft=$NumSearchesLeft-$totalgoodbelowtenKT;
             if($NumSearchesLeft<0){ $NumSearchesLeft=0;}
           }else{
             print "Warning: we found fewer saddles below SPMIN+10KT wst the same spmin?\n";
             print "Number of searches left could be wrong. Go check it!\n";
           }
         }elsif($belowTenKT0 < $belowTenKT) {
           $NumSearchesLeft=$NumSearches;
         }else{
           die "Error: belowTenKT0 > belowTenKT($belowTenKT0 \> $belowTenKT)\n";
         }
       }else{ # the lowest energy saddle is not same
         $NumSearchesLeft=$NumSearches;
       }
     }
   }elsif($total_sp0 == $total_sp){ # find no more good saddles
     if($total_sp0 == 0) { # at the very beginning
       $NumSearchesLeft=$NumSearches;
     }
   }else{
     die "Error: total_sp0 > total_sp($total_sp0 \> $total_sp)\n";
   }
   if($NumSearchesLeft==0) {
     $goodandsound=1;
   }else{
     $goodandsound=0;
   } 
   return ($NumSearchesLeft,$goodandsound);
}

# -----------------------------------------------------------------------------------------------------------------
# $goodandsound=EvaluateSPPopulation($curdir,$NumSearches,$SearchesAlgo,$NumKT,$Temperature,$st)
# -----------------------------------------------------------------------------------------------------------------
sub EvaluateSPPopulation{
  my ($curdir,$NumSearches,$SearchesAlgo,$NumKT,$Temperature,$st)=@_;
  my $goodandsound=0;
  my ($unique_sp,$total_sp,$numsamesaddle,$good_saddles_inprocess,$belowTenKT,$belowTenKT_total,$MinSP,$degeneracy,$energy);
  # count how many unique good saddle points have been found.
   print "Counting good saddle points already found.\n";
   ($unique_sp,$total_sp,$good_saddles_inprocess,$belowTenKT,$belowTenKT_total,$MinSP,$degeneracy,$energy)=CountSP($curdir,$SearchesAlgo,$NumKT,$Temperature,$st);
   print "$unique_sp unique saddle points among total $total_sp good saddle points found.\n";
   if($total_sp >= $NumSearches) { # found more good saddle points
     $goodandsound=1;
   } 
   return $goodandsound;
}

# ----------------------------------------------------------------------------------------------------------------------------
# Statistics of the saddle points distribution
# ($unique_sp,$total_sp,$good_saddles_inprocess,$belowTenKT,$belowTenKT_total,$MinSaddle,$degeneracy,$sp_energy,$st)=
# SaddleStat($curdir,$prdir,$SearchesAlgo,$ConvergenceAlgo,$NumKT,$Temperature,$NumSearches,$unique_sp,$total_sp,$good_saddles_inprocess,$belowTenKT,
# $belowTenKT_total,$MinSaddle,$degeneracy,$sp_energy,$st);
# ----------------------------------------------------------------------------------------------------------------------------
sub SaddleStat{
  my $curdir=shift;
  my $prdir=shift;
  my $SearchesAlgo=shift;
  my $ConvergenceAlgo=shift;
  my $NumKT=shift;
  my $Temperature=shift;
  my $NumSearches=shift;
  my $unique_sp=shift;
  my $total_sp=shift;
  my $good_saddles_inprocess=shift;
  my $belowTenKT=shift;
  my $belowTenKT_total=shift;
  my $MinSaddle=shift;
  my $MaxSaddle=shift;
  my $degeneracy=shift;
  my $sp_energy=shift;
  my $st=shift;
  my $NumSearchesLeft=$st->{"NumSearchesLeft"}[0];
  #print "Saddle analysis (Algo: $ConvergenceAlgo)... ... ";
  if($ConvergenceAlgo == 1) { # Count low energy good saddles only, including repeated.
    if(lc($st->{$prdir}[1]) eq "good") {# a unique good saddle
      if(lc($st->{$prdir}[0]) eq "done") { # a done good saddle
        if($MinSaddle) {
          if(($st->{$prdir}[2] <= $MaxSaddle)) { # a new unique saddle
            #print "$prdir is a new unique saddle below MinSaddle+NumKT.....\n";
            $NumSearchesLeft=$NumSearches; #reset the magi cnumber
            if($st->{$prdir}[2] < $st->{$MinSaddle}[2]) { #a new low, reset the whole thing
              #print "hmmmm,even lower than MinSaddle. Reset the NumSearchesLeft number to $NumSearches... ...\n";
              ($unique_sp,$total_sp,$good_saddles_inprocess,$belowTenKT,$belowTenKT_total,$MinSaddle,$degeneracy,$sp_energy)=CountSP($curdir,$SearchesAlgo,$NumKT,$Temperature,$st);
              $MaxSaddle=$st->{$MinSaddle}[2]+$NumKT*$Temperature/11604.5;
              if($degeneracy->{$prdir}[0]>1) { # default is 1
                #print "However, $degeneracy->{$prdir}[0] processes(including isself) repeat it. So, cut the NumSearchesLeft number ...";
                $NumSearchesLeft=$NumSearchesLeft-$degeneracy->{$prdir}[0]+1; # NumSearchesLeft is reset; itself doesn't contribute to the repeatedness
                if($NumSearchesLeft < 0) {$NumSearchesLeft=0};
              }
              #print "done";
              $UPDATEPROCESS{"comment"}="good\&a new lowest saddle";
            }else{ # in between minsp and maxsp
              #print "it is still above Minsaddle ... ...";
              ($prdirdegeneracy,$st)=CheckRepeatedSP($curdir,$prdir,$numprst,$SearchesAlgo,$st);
              $degeneracy->{$prdir}=[($prdirdegeneracy)];
              if($prdirdegeneracy > 1) {#default is 1
                #print "However, $prdirdegeneracy processes(including istelf) repeat it. So, cut the NumSearchesLeft number ...";
                $NumSearchesLeft=$NumSearchesLeft-$prdirdegeneracy+1; #NumSearchesLeft is reset; itself doesn't contribute to the repeatedness
                if($NumSearchesLeft < 0) {$NumSearchesLeft=0};
                $UPDATEPROCESS{"comment"}="good\&counted\&repeated $prdirdegeneracy times";
              }
              $sp_energy->{$prdir}=$st->{$prdir}[2];
              $unique_sp++;
              $total_sp=$total_sp+$prdirdegeneracy;
              $belowTenKT++;
              $belowTenKT_total=$belowTenKT_total+$prdirdegeneracy;
              #print "done";
            }
            #print "\n";
          }else{ # an high energy unique saddle
            #print "$prdir is a new unique saddle ABOVE MinSaddle+NumKT, so it does not contribute ...";
            ($prdirdegeneracy,$st)=CheckRepeatedSP($curdir,$prdir,$numprst,$SearchesAlgo,$st);
            $degeneracy->{$prdir}=[($prdirdegeneracy)];
            $sp_energy->{$prdir}=$st->{$prdir}[2];
            $unique_sp++;
            $total_sp=$total_sp+$prdirdegeneracy;
            #print "done\n";
            $UPDATEPROCESS{"comment"}="good\&uncounted";
          }
        }else{
          #print "This is a first good saddle ever found.";
          $NumSearchesLeft=$NumSearches;
          $MinSaddle=$prdir;
          $MaxSaddle=$st->{$MinSaddle}[2]+$NumKT*$Temperature/11604.5;
          ($prdirdegeneracy,$st)=CheckRepeatedSP($curdir,$prdir,$numprst,$SearchesAlgo,$st);
          $degeneracy->{$prdir}=[($prdirdegeneracy)];
          if($prdirdegeneracy > 1) {
            #print "However, $prdirdegeneracy processes(including itself) repeat it. So, cut the NumSearchesLeft number ...";
            $NumSearchesLeft=$NumSearchesLeft-$prdirdegeneracy+1;# NumSearchesLeft is reset; itself doesn't contribute to the repeatedness
            if($NumSearchesLeft < 0) {$NumSearchesLeft=0};
          }
          $sp_energy->{$prdir}=$st->{$prdir}[2];
          $unique_sp++;
          $total_sp=$total_sp+$prdirdegeneracy;
          $belowTenKT++;
          $belowTenKT_total=$belowTenKT_total+$prdirdegeneracy;
          #print "... done\n";
          $UPDATEPROCESS{"comment"}="the first good saddle";
        } # end of if(MinSaddle)
      }else{#Good but not done
        if(lc($st->{$prdir}[0]) eq "dynmat") {
          #print "A done good saddle is coming. Better wait for this process to be done or reduce the MaxJobs ...";
        }elsif(lc($st->{$prdir}[0]) eq "killed"){
          #print "A good saddle but the process is killed(thrown way). It could be caused by a fatally failed dynmat calculation(fix it manually)...";
        }else{
          #print "Warning in $prdir, funky process status $st->{$prdir}[0] ...";
        }
        #print "done\n";
        $UPDATEPROCESS{"comment"}="a good saddle is coming. better wait.";
      }
    }elsif($st->{$prdir}[1]=~/repeat/i){# a repeated saddle (just done)
      if(lc($st->{$prdir}[0]) eq "done") {
        $prdirrepeated=substr($st->{$prdir}[1],7);
        if((lc($st->{$prdirrepeated}[0]) eq "done") && (lc($st->{$prdirrepeated}[1]) eq "good")) { # it repeats a done old unique saddle
          $degeneracy->{$prdirrepeated}[0]++;
          $total_sp++;
          if(($st->{$prdirrepeated}[2] <= $MaxSaddle)) { # it should be above minsp
            #print "repeat $prdirrepeated which has been counted. ...";
            $belowTenKT_total++;
            $NumSearchesLeft--;
            if($NumSearchesLeft<0){ $NumSearchesLeft=0; }
            $UPDATEPROCESS{"comment"}="repeat $prdirrepeated(good\&counted)";
          }else{
            #print "repeat $prdirrepeated which is good, but over minsaddle+numKT and thus not counted ...";
            $UPDATEPROCESS{"comment"}="repeat $prdirrepeated(good\&uncounted)";
          }
        }else{
          #print "repeat $prdirrepeated which is ".$st->{$prdirrepeated}[0]." and ".$st->{$prdirrepeated}[1].", so it has NOT been counted ...";
          $UPDATEPROCESS{"comment"}="repeat $prdirrepeated(".$st->{$prdirrepeated}[0]."\&".$st->{$prdirrepeated}[1].")";
        }
      }else{
        if(lc($st->{$prdir}[0]) eq "killed") {
           #print "A repeated saddle is thrown way(killed) for some reason which you must know(go check it)...";
           if($UPDATEPROCESS{"comment"} eq "---"){$UPDATEPROCESS{"comment"}="repeated\&thrown away.check it!";
           }else{$UPDATEPROCESS{"comment"}=$UPDATEPROCESS{"comment"}." repeated\&thrown away.check it!";}
        }else{
           #print "Warning in $prdir, funky process status $st->{$prdir}[0] ...";
           $UPDATEPROCESS{"comment"}="Warning, funky process status $st->{$prdir}[0]";
        }
      }
      #print "done\n";
    }elsif(lc($st->{$prdir}[1]) eq "promising") {
      #print "a $st->{$prdir}[1] saddle. It has not contributed yet.\n";
    }else{
      #print "a $st->{$prdir}[1] saddle etc. it doesn't contribute.\n";
    } # end of ConvergenceAlgo==1 if block
  }elsif($ConvergenceAlgo ==0){ # Count the total good saddles only, including repeated.
    if(lc($st->{$prdir}[1]) eq "good") {# a unique good  saddle
      if(lc($st->{$prdir}[0]) eq "done") { # a done good saddle
        #print "$prdir is a new unique saddle ... ... ";
        $UPDATEPROCESS{"comment"}="good\&new";
        ($prdirdegeneracy,$st)=CheckRepeatedSP($curdir,$prdir,$numprst,$SearchesAlgo,$st);
        $degeneracy->{$prdir}=[($prdirdegeneracy)];
        if($prdirdegeneracy > 1) {#default is 1
          #print "has been repeated $prdirdegeneracy times(including itself). \n";
          $NumSearchesLeft=$NumSearchesLeft-$prdirdegeneracy;
          if($NumSearchesLeft < 0) {$NumSearchesLeft=0};
          $UPDATEPROCESS{"comment"}=$UPDATEPROCESS{"comment"}."repeated $prdirdegeneracy times";    
        }
        $unique_sp++;
        $total_sp=$total_sp+$prdirdegeneracy;
      }else{#Good but not done, it must be doing dynmat or be killed
        if((lc($st->{$prdir}[0]) ne "killed") && (lc($st->{$prdir}[0]) ne "dynmat")) {
          #print "Warning in $prdir, funky process status $st->{$prdir}[0]!\n";
          $UPDATEPROCESS{"comment"}="Warning: funky process status $st->{$prdir}[0]!";
        }
      }
    }elsif($st->{$prdir}[1]=~/repeat/i){# a repeated saddle (just done)
      if(lc($st->{$prdir}[0]) eq "done") {
        $prdirrepeated=substr($st->{$prdir}[1],7);
        if((lc($st->{$prdirrepeated}[0]) eq "done") && (lc($st->{$prdirrepeated}[1]) eq "good")) { # it repeats a done good unique saddle
          #print "$prdir repeats $prdirrepeated which is done and good, and thus counted.\n";
          $NumSearchesLeft--;
          if($NumSearchesLeft<0){ $NumSearchesLeft=0; }
          $total_sp++;
          $UPDATEPROCESS{"comment"}="repeat $prdirrepeated(good\&counted)";
        }else{
          #print "$prdir repeats $prdirrepeated which is a bad or killed saddle, or just not done yet. NOT counted.\n";
          $UPDATEPROCESS{"comment"}="repeat $prdirrepeated(bad\&uncounted)";
        }
      }else{
        #print "Warning in $prdir, funky process status $st->{$prdir}[0]!\n";
        $UPDATEPROCESS{"comment"}="Warning, funky process status $st->{$prdir}[0]!";
      }
    }else{
      #print "a $st->{$prdir}[1] saddle etc. it contributes nothing for the time being.\n";
    } # end of ConvergenceAlgo==0 if block    
  }else{
    die "Error in SaddleStat ; a wrong ConvergenceAlgo: $ConvergenceAlgo\n";
  }# end of the outmost if block
  $st->{"NumSearchesLeft"}[0]=$NumSearchesLeft;
  return($unique_sp,$total_sp,$good_saddles_inprocess,$belowTenKT,$belowTenKT_total,$MinSaddle,$MaxSaddle,$degeneracy,$sp_energy,$st);
}

# ---------------------------------------------------------------------------------------------------------
# Update the st info after going through each folder in st00xx [($numprfin,$numprrun,$numrunjobs,$jobs,$st)=
# CheckSt($curdir,$quench,$dynmat,$numprfin,$numprrun,$numrunjobs,$MaxJobs,$NumSearches,$UseKDB,$SearchesAlgo,$ConvergenceAlgo,
# $population,$NumKT,$Temperature,$DisplaceAlgo,$DisplaceRange,$NN_rcut,$MaxCoordNum,$nowait,$jobs,$st,$st_old)]
# ---------------------------------------------------------------------------------------------------------
sub CheckSt{
  my $curdir=shift;
  my $quench=shift;
  my $dynmat=shift;
  my $numprfin=shift;
  my $numprrun=shift;
  my $numrunjobs=shift;
  my $MaxJobs=shift;
  my $NumSearches=shift;
  my $UseKDB=shift;
  my $SearchesAlgo=shift;
  my $ConvergenceAlgo=shift;
  my $Population=shift;
  my $NumKT=shift;
  my $Temperature=shift;
  my $DisplaceAlgo=shift;
  my $DisplaceRange=shift;
  my $NN_rcut=shift;
  my $MaxCoordNum=shift;
  my $nowait=shift;
  my $jobs=shift;
  my $st=shift;
  my $st_old=shift;
  my $enoughsearches=0;
  my ($i,$prdir,$numallowedjobs,$numprst,$numprst_new,$targetsp,$targetmin,$targetdynmat,$success,$NumSearchesLeft,$NumPrsIn);
  $numprst=$st->{"numprst"}[0];
  #evaluate the saddle population
  if($st->{"NumSearchesLeft"}[0] == 0) { # no need to look at processes of the state.
    $enoughsearches=1; 
  }else{
    SWITCH: {
      if($ConvergenceAlgo == 0) {last SWITCH;} # updated at each step, so it is already taken care of.
      if($ConvergenceAlgo == 1) {last SWITCH;} # updated at each step, so it is already taken care of.
      if($ConvergenceAlgo == 2) {($NumSearchesLeft,$enoughsearches)=EvaluateLowEnergySPTenKT($curdir,$NumSearches,$SearchesAlgo,$NumKT,$Temperature,$population,$st_old,$st);
                            $st->{"NumSearchesLeft"}[0]=$NumSearchesLeft; #just to update st hash, not for any use in this subroutine.
                            last SWITCH;}
      if($ConvergenceAlgo == 3) {$enoughsearches=EvaluateLowEnergySPPopulation($curdir,$NumSearches,$SearchesAlgo,$NumKT,$Temperature,$population,$st);
                            last SWITCH;}
    }
  }
  if($enoughsearches) {
   # check and mark the NumSearchesLeft number showing no more processes are needed.
   if($st->{"NumSearchesLeft"}[0] != 0) {$st->{"NumSearchesLeft"}[0] == 0;} # now it is only for ConvergenceAlgo 3
   #check if it is required to wait for finishing quench and dynmat jobs for promising saddle points
   if(!$nowait){
    $NumPrsIn=AnalyzeSt($st);
    if($NumPrsIn->{"quench"} == 0 && $NumPrsIn->{"dynmat"} == 0){
      $nowait=1; # we are done.
    }else{
      print "It looks like enough saddle points have been found, but there are a few unknown saddles.\n";
      print "We have to halt searches until the current quench and dynmat processes are finished.\n";
      print "number of prs in quench: ".$NumPrsIn->{"quench"}.", number of prs in dynmat: ".$NumPrsIn->{"dynmat"}."\n";
    }
   }
   if($nowait) {
    # kill all the running processes
    print "We are satisfied with the saddle points found. Kill all the running saddle processes.\n";
    for($i=1;$i<=$numprst;$i++) {
      $prdir=DirName("pr",$i);
      $targetsp=$curdir."/".$prdir;
      if(lc($st->{$prdir}[0]) ne "done" && lc($st->{$prdir}[0]) ne "killed") {
        if (lc($st->{$prdir}[0]) eq "search") {
          if($jobs->{$targetsp}[0] eq "running" || $jobs->{$targetsp}[0] eq "queue" || $jobs->{$targetsp}[0] eq "submitted") {
            $success=KillJobs($targetsp, $jobs->{$targetsp}[1]);
            #print "tried killing a job in $targetsp jobid:".$jobs->{$targetsp}[1]."\n";
            #if(!$success) {print "can not kill job ".$jobs->{$targetsp}[1]. "in $targetsp.\n";}
          }
          if(-e "$targetsp/OUTCAR") { CleanupJobs($targetsp,0,$SearchesAlgo); }
        }elsif(lc($st->{$prdir}[0]) eq "quench") {
          $targetmin=$targetsp."/".$quench."/min1";
          if($jobs->{$targetmin}[0] eq "running" || $jobs->{$targetmin}[0] eq "queue" || $jobs->{$targetmin}[0] eq "submitted") {
            $success=KillJobs($targetmin, $jobs->{$targetmin}[1]);
            #print "tried killing a job in $targetmin jobid:".$jobs->{$targetmin}[1]."\n";
            #if(!$success) {print "can not kill job ".$jobs->{$targetmin}[1]. "in $targetmin.\n";}
          }    
          if(-e "$targetmin/OUTCAR") {CleanupJobs($targetmin, 0, "minimization");}
          $targetmin=$targetsp."/".$quench."/min2";
          if($jobs->{$targetmin}[0] eq "running" || $jobs->{$targetmin}[0] eq "queue" || $jobs->{$targetmin}[0] eq "submitted") {
            $success=KillJobs($targetmin, $jobs->{$targetmin}[1]);
            #print "tried killing a job in $targetmin jobid:".$jobs->{$targetmin}[1]."\n";
            #if(!$success) {print "can not kill job ".$jobs->{$targetmin}[1]. "in $targetmin.\n";}
          }
          if(-e "$targetmin/OUTCAR") {CleanupJobs($targetmin, 0, "minimization");}
        }elsif(lc($st->{$prdir}[0]) eq "dynmat") {
          $targetdynmat=$targetsp."/".$dynmat;
          if($jobs->{$targetdynmat}[0] eq "running" || $jobs->{$targetdynmat}[0] eq "queue" || $jobs->{$targetdynmat}[0] eq "submitted") {
            $success=KillJobs($targetdynmat, $jobs->{$targetdynmat}[1]);
            #print "tried killing a job in $targetdynmat jobid:".$jobs->{$targetdynmat}[1]."\n";
            #if(!$success) {print "can not kill job ".$jobs->{$targetdynmat}[1]. "in $targetdynmat.\n";}
          }
          if(-e "$targetdynmat/OUTCAR") {CleanupJobs($targetdynmat, 0, "dynmat");}
        }else{
          #print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
          die "Error in CheckSt: unrecognizable st status: $st->{$prdir}[0] \n";
          #print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
        } 
      } #end of the outmost "if"
    } #end of the for loop
    # check the status of the possible dynmat calculation for the state initial config
    if(lc($st->{"dynmat"}[0]) eq "done") {
      $st->{"status"}[0]="done";
      print "All right, we are ready to make a KMC move\n";
    }elsif(lc($st->{"dynmat"}[0]) eq "killed") {
      print "No more saddle processes are needed. However, the dynmat job for the system has been killed. Deadlock!\n";
    }else{
      print "No more saddle processes are needed. Just wait for the dynmat job for the system to finish off\n";
    }
   }else{
     print "We need to wait for quench/dyn jobs to be finished before moving to next step.\n";
     print "Before that, nerither new searches nor unfinshed searches will be continued.\n";
   }
  }else{ # we need to do more saddle point searches if there is any quota available in the queue
    print "Number of searches left is ".$st->{"NumSearchesLeft"}[0].". ";
    $numallowedjobs=$MaxJobs-$numrunjobs;
    if($numallowedjobs > 0) {
      print "Submitting $numallowedjobs jobs.\n";
      $numprst_new=$numprst+$numallowedjobs;
      for($i=$numprst+1;$i<=$numprst_new;$i++) {
        $prdir=DirName("pr", $i);
        ($numrunjobs,$jobs,$st)=StartOneSP($curdir,$prdir,$SearchesAlgo,$numrunjobs,$MaxJobs,$UseKDB,$DisplaceAlgo,$DisplaceRange,$NN_rcut,$MaxCoordNum,$jobs,$st);
        $numprrun++;
      }
      $st->{"status"}[0]="running";
      $st->{"numprst"}[0]=$numprst_new;
#      print "done.\n";
    }elsif($numallowedjobs < 0) {
      print "There are more jobs running than allowed. Check it.\n";
    }else{
      print "The akmc queue is full.\n";
    }
  }
  return ($numprfin,$numprrun,$numrunjobs,$jobs,$st);
}

# -------------------------------------------------------------------------------
# Analyze the state information and figure out how many processes in each status 
# $NumPrsIn=AnalyzeSt($st);
# -------------------------------------------------------------------------------
sub AnalyzeSt{ 
  my $st=shift;
  my (%NumPrsIn, $i,$numprst,$prdir);
  %NumPrsIn=("search" => 0,"quench" => 0,"dynmat" => 0,"done" => 0,"killed" => 0);
  $numprst=$st->{"numprst"}[0];
  for($i=1;$i<=$numprst;$i++){
    $prdir=DirName("pr",$i);
    if(lc($st->{$prdir}[0]) eq "search"){
      $NumPrsIn{"search"}++;
    }elsif(lc($st->{$prdir}[0]) eq "quench"){
      $NumPrsIn{"quench"}++;
    }elsif(lc($st->{$prdir}[0]) eq "done"){
      $NumPrsIn{"done"}++;
    }elsif(lc($st->{$prdir}[0]) eq "dynmat"){
      $NumPrsIn{"dynmat"}++;
    }elsif(lc($st->{$prdir}[0]) eq "killed"){
      $NumPrsIn{"killed"}++;
    }else{
      die "In AnalyzeSt, unknown status of $prdir\n";
    }
  }
  return \%NumPrsIn;
}

# --------------------------------------------------------------------------------------------------------
# Update state info after working on every folder in st00xx 
# [%rt=WriteEventsTable($curdir,$DynMatFile,$turnoff_dynmat,$Temperature,$RateTableFile,$st,$Prefactor)]
# --------------------------------------------------------------------------------------------------------
sub WriteEventsTable{
  my ($curdir,$DynMatFile,$turnoff_dynmat,$Temperature,$etfile,$st,$Prefactor)=@_;
  my %rt=();
  my ($energy_init,$dynmat_init,$dynmat_sp,$prefactors,$barrier,$finalstate,$rate,$numprst,$i,$j);
  my ($curdyndir,$prdir);
  my $numgoodsp=0;
  if(lc($st->{"status"}[0]) ne "done") {die "In $curdir, the state file shows the state is not done yet.\n";}
  if(!$turnoff_dynmat){
    $curdyndir=$curdir."/dynmat/final/".$DynMatFile;
    if( -e $curdyndir){
      $dynmat_init=$curdyndir;
    }else{
      die "In WriteEventsTable: no initial dynmat :$curdyndir \n";
    }
  }else{
    $Prefactor || die "Dynmat is turned off without an alternative.\n";
  }
  if( -e $etfile) {
    print "In building the event table, the event table file $etfile exists and will be overwritten.\n";
  }
  open (OUT, ">$curdir/$etfile");
  #print OUT  $curdir."  ".$Temperature." Kelvin\n";
  print OUT  substr($curdir,-6)."  ".$Temperature." Kelvin\n"; #st000x is the key to write.
  $energy_init=$st->{$curdir}[0];
  $numprst=$st->{"numprst"}[0];  
  for($i=1;$i<=$numprst;$i++){
    $prdir=DirName("pr",$i);
    if($st->{$prdir}[1] eq "good" && $st->{$prdir}[0] eq "done"){
      $numgoodsp++;
      $rt{$prdir}=[()];
      $barrier=$st->{$prdir}[2]-$energy_init;
      if(!$turnoff_dynmat){
        $dynmat_sp=$curdir."/".$prdir."/dynmat/final/".$DynMatFile;
        $prefactors=Calc_Prefactor($dynmat_init,$dynmat_sp);
      }else{
        $prefactors=$Prefactor;
      }
      $final_state=$curdir."/".$prdir."/".$st->{$prdir}[3];
      # put barrier and prefactor information into st hash
      $rate=Calc_Rate($barrier, $prefactors, $Temperature);
      push @{$rt{$prdir}}, $barrier, $prefactors, $final_state, $rate;
      print OUT  $prdir." "."@{$rt{$prdir}}"."\n"; 
    }else{
      next;
    }    
  }
  close OUT;
  if($numgoodsp == 0) {
    die "In $curdir, no good saddle points were found. The folder or the state file is corrupted.\n";
  }
  if(exists($rt{""})){delete $rt{""};}
  return \%rt;
}

# --------------------------------------------------------------------------------------------------------
# Update state info after working on every folder in st00xx 
# [%rt=ReadEventTable($curdir,$Temperature,$RateTableFile)]
# --------------------------------------------------------------------------------------------------------
sub ReadEventTable{
  my ($curdir,$Temperature,$RateTableFile)=@_;
  my $line="";
  my @line=();
  my %rt=();
  my $temp=0.0;
  my $j=0;
  my $filename=$curdir."/".$RateTableFile;
  my $stname="";  
  my $i=0;
  if(-e $filename){
    open (IN, "<$filename") || die "In ReadEventTable: cannot open $filename\n";
  }else{
    print "Warning: there is no $RateTableFile in the done state of $curdir\n";
    $rt=0;
    return $rt;
  }
  while($line=<IN>){
    $line=~s/^\s+//;
    if($line eq "" || $line eq "\n"){next;} # ignore empty lines
    chomp($line=lc($line));
    @line=split(/\s+/,$line);
    #if(lc($line[-1]) eq "excluded"){next;} # manual control; WriteEventTable will not care.
    $rt{$line[0]}=[()];
    for($j=1;$j<@line;$j++){
      push @{$rt{$line[0]}}, $line[$j];
    }
    $i++;
  }
  if(!$i){$rt=0; return $rt;}
  $stname=lc(substr($curdir,-6)); #st000x is the key
  if(exists($rt{$stname})){
    $temp=$rt{$stname}[0];
    delete $rt{$stname};
  }else{die "No $stname($curdir) information in the Event\n";}
  if($Temperature != $temp){
    for $j (keys %rt){
      $barrier=$rt{$j}[0];
      $prefactor=$rt{$j}[1];
      $rate=Calc_Rate($barrier, $prefactor, $Temperature);
      $rt{$j}[3]=$rate;
    }
  }
  close IN;
  if(exists($rt{""})){delete $rt{""};}
  return \%rt;
}

sub ReadEventTableBolEqu{
  my ($curdir,$Temperature,$RateTableFile)=@_;
  my $line="";
  my @line=();
  my %rt=();
  my $temp=0.0;
  my $j=0;
  my $filename=$curdir."/".$RateTableFile;
  my $stname="";
  my $i=0;
  if(-e $filename){
    open (IN, "<$filename") || die "In ReadEventTableBolEqu: cannot open $filename\n";
  }else{
    print "Warning: there is no $RateTableFile in the done state of $curdir\n";
    $rt=0;
    return $rt;
  }
  while($line=<IN>){
    $line=~s/^\s+//;
    if($line eq "" || $line eq "\n"){next;} # ignore empty lines
    chomp($line=lc($line));
    @line=split(/\s+/,$line);
    if(lc($line[-1]) eq "excluded"){next;} # manual control; WriteEventTable will not care.
    $rt{$line[0]}=[()];
    for($j=1;$j<@line;$j++){
      push @{$rt{$line[0]}}, $line[$j];
    }
    $i++;
  }
  close IN;
  if(!$i){$rt=0; return $rt;}
  if(exists($rt{""})){delete $rt{""};}
  $stname=lc(substr($curdir,-6)); #st000x is the key
  if(exists($rt{$stname})){
    $temp=$rt{$stname}[0];
  }else{die "No $stname($curdir) information in the Event\n";}
  if($Temperature != $temp){# have to rebuild the BolEqu
    $rt=0;
    return $rt;
  }
  return \%rt;
}

# ----------------------------------------------------------------------
# ($selected,$elapse,$enforced,$time_saved)=SelectEvent($rt,$curdir,$prest,$randomnum,$randomnum2,$enforced,
#  $time_saved,$Rdiffmax,$RateTableFile,$Temperature,$RateTableFileBolEqu,$akmc_step,$curstep,$BoltzmanEqu)
# ----------------------------------------------------------------------
sub SelectEvent{
  my ($rt,$curdir,$prest,$rdnum,$rdnum2,$enforced,$t_saved,$Rdiffmax,$RTFile,$Temp,$RTFileBE,$akmc_step,$curstep,$BEqu)=@_;
  my ($rate_total,$ratem_1,$ratem,$selected,$time_elapse,$dominant,$i,$BolEqu,$rtbe);
  $rate_total=$ratem_1=$ratem=0;
  $dominant="";
  $BolEqu=$rtbe=0;
  if($enforced eq ""){
    # check if there is a file named $RTFileBE
    if(-e "$curdir/$RTFileBE"){#manual control is possible
      print "We are trying to read rate information from BolEqu file $RateTableFileBolEqu\n";
      $rtbe=ReadEventTableBolEqu($curdir,$Temp,$RTFileBE);
    }
    # check if curdir is in boltzman equilibrium with other state
    if(!$rtbe){
      $rate_total=Sum_HashArrayDim($rt,3); # the forth element is the rate constant
      $dominant=CheckDominance($rt,3,$rate_total,$BoltzmanEqu);
      if($dominant ne ""){ 
       ($BolEqu,$rtbe)=CheckBoltzman($rt,$curdir,$prest,$BEqu,$dominant,$Rdiffmax,$RTFile,$Temp,$RTFileBE,$akmc_step,$curstep);
      }
    }else{
      $BolEqu=1;
    }
    # handle event selection
    if($BolEqu){
      ($selected,$time_elapse,$enforced,$t_saved)=HandleBoltzmanEquil($rtbe,$rdnum,$rdnum2,$prest,$curdir);
      if($enforced ne ""){
        if($rate_total==0){$rate_total=Sum_HashArrayDim($rt,3);} # the forth element is the rate constant
        $time_elapse=KMC_Time_possion($rate_total,$rdnum);
      }
    }else{ # just a usual selection
      #if($rate_total==0){$rate_total=Sum_HashArrayDim($rt,3);} # the forth element is the rate constant
      $time_elapse=KMC_Time_possion($rate_total,$rdnum);
      print "Selecting a process using kinetic Monte Carlo.\n";
      $selected=KMC_Pickup_Event($rt,3,$rate_total,$rdnum2);
    }
  }else{
    print "The process $enforced is enforced with transition time $t_saved; do it.\n";
    $selected=$enforced;
    $time_elapse=$t_saved;
    $enforced="";
    $t_saved=0;
  }
  return ($selected, $time_elapse,$enforced,$t_saved);
}

# --------------------------------
# $sum=Sum_HashArrayDim($rt,$dim)
# --------------------------------
sub Sum_HashArrayDim{
  my ($rt,$dim)=@_;
  my ($sum,$i);
  $sum=0;
  for $i (keys %$rt){
    $sum+=$rt->{$i}[$dim];
  }
  return $sum;
}

# ---------------------------------------------------------
# $dominant=CheckDominance($rt,$dim,$rate_total,$threshold)
# ---------------------------------------------------------
sub CheckDominance{
  my ($rt,$dim,$rate_total,$threshold)=@_;
  my $dominant="";
  if($rate_total==0){
    print "Warning: in CheckDominance, rate_total=0. Divided by zero is coming. We quit here.\n";
    return $dominant;
  }
  for $i (sort keys %$rt){
    if($rt->{$i}[$dim]/$rate_total >= $threshold){
      $dominant=$i;
      last;
    }
  } 
  return $dominant;
}

# -----------------------------------------------------
# $time_elapse=KMC_Time_possion($rate_total,$randomnum)
# -----------------------------------------------------
sub KMC_Time_possion{
  my ($rate_total,$randomnum)=@_;
  my $time_elapse;
  if($rate_total>0){
    $time_elapse=-1.0*log($randomnum)/$rate_total;
  }else{
    print "Warning in KMC_Time_possion: total rate is zero. Divid by zero is coming! We don't do the time\n";
    $time_elapse=0;
  }
  return $time_elapse;
}

# -----------------------------------------------------------
# $selected=KMC_Pickup_Event($rt,$dim,$rate_total,$randomnum)
# -----------------------------------------------------------
sub KMC_Pickup_Event{
  my ($rt,$dim,$rate_total,$randomnum)=@_;
  my ($selected,$reducedRate,$ratem_1,$ratem,$i);
  $selected="";
  $reducedRate=$rate_total*$randomnum;
  $ratem_1=$ratem=0;
  for $i (sort keys %$rt){
    $ratem=$ratem+$rt->{$i}[$dim];
    if($reducedRate>$ratem_1 && $reducedRate<=$ratem){
      $selected=$i;
      last;
    }
    $ratem_1=$ratem;
  }
  if($selected eq ""){
    print "$total rate: $rate_total; rate(m-1): $ratem_1\n";
    die "Fatal error in KMC_Pickup_Event: Weird, no event is selected.\n";
  }
  return $selected;
}

# ---------------------------------------------------------------------------------------------------------
# ($BolEqu,$rt)=CheckBoltzman($rt,$curdir,$prest,$BoltzmanEqu,$dominant,$Rdiffmax,$RateTableFile,$Temperature,$RateTableFileBolEqu,$akmc_step,$curstep)
# ---------------------------------------------------------------------------------------------------------
sub CheckBoltzman{
  my ($rt,$curdir,$prest,$BoltzmanEqu,$dominant,$Rdiffmax,$RateTableFile,$Temperature,$RateTableFileBolEqu,$akmc_step,$curstep)=@_;
  my ($BolEqu,$same,$rtpre,$rtpre_total,$i,$j,$dominant2,$deltaE,$line,$dummy);
  $BolEqu=$rtpre=0;
  if($prest eq ""){ # there is no previous state,i.e., we are at the very first step
    return ($BolEqu,$rtpre); 
  } 
  $same=CompareTwoPOSCAR($prest."/POSCAR", $rt->{$dominant}[2]."/POSCAR", $Rdiffmax);
  if($same){ # check if the prest also has a dominant process
    $rtpre=ReadEventTable($prest,$Temperature,$RateTableFile);
    if($rtpre==0){return ($BolEqu,$rtpre);}
    print "prest=".$prest."\n";
    $rtpre_total=Sum_HashArrayDim($rtpre,3);
    print "rtpre_total=$rtpre_total\n";
    $dominant2=CheckDominance($rtpre,3,$rtpre_total,$BoltzmanEqu);
    $same=0;
    if($dominant2 ne ""){ 
      print $curdir."/POSCAR"."***".$rtpre->{$dominant2}[2]."/POSCAR"."\n";
      $same=CompareTwoPOSCAR($curdir."/POSCAR", $rtpre->{$dominant2}[2]."/POSCAR", $Rdiffmax);
      if($same){ # bingo
        $deltaE=$akmc_step->{$curstep}[3]-$akmc_step->{$curstep-1}[3];
        $line="";
        for $i (sort keys %$rt){
          if($i eq $dominant){next;}
          $j=@{$rt->{$i}}-2;
          if($deltaE<=0){ # curdir is of the lower energy
            $line.=$curdir."/$i "." @{$rt->{$i}}[0 .. $j] ";
            $line.=Calc_Rate($rt->{$i}[0], $rt->{$i}[1], $Temperature)."\n";
          }else{
            $dummy=$rt->{$i}[0]+$deltaE;
            $line.=$curdir."/$i ".$dummy." @{$rt->{$i}}[1 .. $j] ";
            $line.=Calc_Rate($dummy, $rt->{$i}[1], $Temperature)."\n";
          }
        }
        for $i (sort keys %$rtpre){
          if($i eq $dominant2){next;}
          $j=@{$rtpre->{$i}}-2;
          if($deltaE>0){ # prest is of the lower energy
            $line.=$prest."/$i "." @{$rtpre->{$i}}[0 .. $j] ";
            $line.=Calc_Rate($rtpre->{$i}[0], $rtpre->{$i}[1], $Temperature)."\n";
          }else{
            $dummy=$rtpre->{$i}[0]+$deltaE;
            $line.=$prest."/$i ".$dummy." @{$rtpre->{$i}}[1 .. $j] ";
            $line.=Calc_Rate($dummy, $rtpre->{$i}[1], $Temperature)."\n";
          }
        }
        open(NEWRT, ">$curdir/$RateTableFileBolEqu") || die "In selectevent, cannot open $curdir/$RateTableFileBolEqu\n";
          print NEWRT substr($curdir,-6)." ".$Temperature." ".$dominant."\n".$line;
        close NEWRT;
        open(NEWRT, ">$prest/$RateTableFileBolEqu") || die "In selectevent, cannot open $prest/$RateTableFileBolEqu\n";
          print NEWRT substr($prest,-6)." ".$Temperature." ".$dominant2."\n".$line;
        close NEWRT;
        $rtpre=ReadEventTableBolEqu($curdir,$Temperature,$RateTableFileBolEqu);
        $BolEqu=1;
      }
    } # end of "if($dominant2" block
  }
  return ($BolEqu,$rtpre);
}

# ---------------------------------------------------------------------------------------------------------
# ($selected,$time_elapse,$enforced,$time_saved)=HandleBoltzmanEquil($rt,$randomnum,$randomnum2,$prest,$curdir,$enforced,$time_saved);
# ---------------------------------------------------------------------------------------------------------
sub HandleBoltzmanEquil{
  my ($rt,$randomnum,$randomnum2,$prest,$curdir)=@_;
  my ($selected,$time_elapse,$enforced,$time_saved);
  my ($stdir,$prdir,$dominant);
  $stdir=lc(substr($curdir, -6));
  if(exists($rt->{$stdir})){
    $dominant=$rt->{$stdir}[1];
    delete($rt->{$stdir});
  }else{
    die "Fatal error in HandleBoltzmanEquil, rt hash has to have $stdir keys.\n;"
  }
  $rate_total=Sum_HashArrayDim($rt, 3);
  $time_elapse=KMC_Time_possion($rate_total,$randomnum);
  $selected=KMC_Pickup_Event($rt, 3,$rate_total,$randomnum2);
  $prdir=substr($selected,-6);
  if(!(lc($prdir)=~/^pr\d{4}$/)){die "In HandleBoltzmanEquil, unrecognizable $selected.\n";}
  $stdir=substr($selected,0,-7);
  if($stdir eq $curdir){
    $selected=$prdir;
    $enforced="";
    $time_saved=0;
    print "Selected $selected in $curdir.\n";
  }else{
    $selected=$dominant;
    $enforced=$prdir;
    $time_saved=$time_elapse;
    $time_elapse=0;
    print "Selected $dominant in $curdir.  Remember enforcing $enforced in $prest.\n";
  }
  return ($selected,$time_elapse,$enforced,$time_saved);
}

# ---------------------------------------------------------------------------------------------------------
# WriteAkmcStep($AkmcFile,$akmc_step,$initial)
# ---------------------------------------------------------------------------------------------------------
sub WriteAkmcStep{
  my ($AkmcFile,$akmc_step,$initial)=@_;
  my $AkmcFile=shift;
  my $akmc_step=shift;
  my $initial=shift;
  my $j="";
  my @newline=();
  if($initial) {
    open(OUT, ">$AkmcFile") || die "Error in opening $AkmcFile file for writing KMC step information.\n";
    print OUT "Step stnumber property address stenergy(eV) ftenergy(eV) ChosenPr steptime(s) totaltime(s)"."\n";
  }else{
    open(OUT, ">>$AkmcFile") || die "Error in opening $AkmcFile file for writing KMC step information.\n";
  }
  sub numerically { $a <=> $b}
  @newline= (sort numerically keys %$akmc_step);
  $j=@newline;
  if(!$j) {die "Error in WriteAkmcStep: the akmc_step has no elements.\n";}
  print OUT $newline[$j-1]." "."@{$akmc_step->{$newline[$j-1]}}"."\n";
  close OUT;
}

# ----------------------------------------------------------------------------------------------------------------------------
# Update state by making a move
# ($curstep,$st,$akmc_step)=UpdateSt($curstep,$repeat,$selected,$randomnum,$time,$StFile,$StEnergyFile,$AkmcFile,$st,$akmc_step);
# ----------------------------------------------------------------------------------------------------------------------------
sub UpdateSt{
  my $curstep=shift;
  my $repeat=shift;
  my $selected=shift;
  my $randomnum=shift;
  my $time_elapse=shift;
  my $StFile=shift;
  my $StEnergyFile=shift;
  my $AkmcFile=shift;
  my $st=shift;
  my $akmc_step=shift;
  my $finalst="";
  my ($energy,$force);
  my ($curdir,$stname,$finalstPOSCAR,$POSCAR_j,$total_time,$prestep);
  # save current state information if it is not repeated.
  $curdir=$akmc_step->{$curstep}[2];
  $stname=$akmc_step->{$curstep}[0];
  if($stname ne substr($curdir, -6)) {die "$akmc_step for $curstep is messed up.\n";}
  $energy=$st->{$curdir}[0]; # ReadSt puts whatever in the first line into $st
  if(!($energy=~/^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$/)){die "State energy ($energy) in $curdir is messed up\n";}
  if(!$repeat){
    WriteSt($curdir,$StFile,$st);
    AppendStEnergyFile($curdir,$StEnergyFile,$energy);
  }
  # check out the final state information
  $finalst=$curdir."/".$selected."/".$st->{$selected}[3];
  $finalstPOSCAR=$finalst."/POSCAR";
  #($energy,$force)=GetEF($finalst);
  $energy=$st->{$selected}[4];
  # update POSCAR in the akmc root folder
  print "Process $selected selected from the random number: $randomnum\n";
  system "cp $finalstPOSCAR POSCAR";
  # save the final state energy
  $akmc_step->{$curstep}[4]=$energy;
  $akmc_step->{$curstep}[5]=$selected;
  # accumulate the total time
  $akmc_step->{$curstep}[6]=$time_elapse;
  $akmc_step->{$curstep}[7]=$akmc_step->{$curstep}[7] + $time_elapse;
  $total_time=$akmc_step->{$curstep}[7];
  print "Time for this step is: $time_elapse\n";
  print "Total time elapsed is: $total_time\n";
  # save the done step information
  if($curstep == 1) {
    WriteAkmcStep($AkmcFile, $akmc_step, 1);
  }else {
    WriteAkmcStep($AkmcFile, $akmc_step, 0);
  }
  # put a new line to the akmc_step file and keeps akmc_step from having more than two elements
  $prestep=$curstep-1;
  if(exists($akmc_step->{$prestep})) {delete $akmc_step->{$prestep};}
  if(exists($akmc_step->{""})) {delete $akmc_step->{""};}
  # $curstep=1+(keys %$akmc_step);
  $curstep=$curstep+1;
  # $curstep=DirName("Step",$curstep);
  $akmc_step->{$curstep}=[()];
  push @{$akmc_step->{$curstep}}, "2bcreated", "unknown", "unknown",$energy,"na","na","na",$total_time;
  return ($curstep,$st,$akmc_step);
}
