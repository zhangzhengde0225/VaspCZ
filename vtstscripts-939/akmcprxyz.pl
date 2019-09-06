#!/usr/bin/env perl
#;-*- Perl -*-

use FindBin qw($Bin);
use lib "$Bin";

# Read st.dat file
open ST,"st.dat";
while($line=<ST>){
   $line=~s/^\s+//g;
   @line=split /\s+/,$line;
   $prdir=$line[0],"\n";
   if($prdir =~ m/^pr\d{4}$/){
     push (@prdir, "$prdir");}
 } 

# Temporary directory for making movie
$tmpdir="akmc_pr_tmp";
`rm -rf $tmpdir pr.xyz`;
`mkdir $tmpdir`;
`touch $tmpdir/pr.xyz`;

# Make movie of current configurations
foreach $dir (@prdir){
  print $dir."\n";
  if((-e "$dir/final/CONTCAR") && (-s "$dir/final/CONTCAR")){
    `cp $dir/final/CONTCAR $tmpdir/POSCAR`;
  }elsif((-e "$dir/final/POSCAR") && (-s "$dir/final/POSCAR")){
    `cp $dir/final/POSCAR $tmpdir/POSCAR`;
  }elsif((-e "$dir/CONTCAR") && (-s "$dir/CONTCAR")){
    `cp $dir/CONTCAR $tmpdir/POSCAR`;
  }elsif((-e "$dir/POSCAR") && (-s "$dir/POSCAR")){
    `cp $dir/POSCAR $tmpdir/POSCAR`;
  }
  `cd $tmpdir; $Bin/pos2con.pl POSCAR; $Bin/con2xyz.pl POSCAR.con; rm POSCAR.con; cat POSCAR.xyz | sed 's/Generated with con2xyz/$dir/' >> pr.xyz; rm POSCAR.xyz; rm POSCAR`;

}
`mv $tmpdir/pr.xyz .; rm -rf $tmpdir;`;
#`gzip pr.xyz`

