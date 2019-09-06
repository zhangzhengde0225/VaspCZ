#!/usr/bin/env perl
#;-*- Perl -*-

  die "input a POSCAR or a CONTCAR \n" if @ARGV > 1 ;
  $pos = $ARGV[0] ;  
  open IN , $pos ;
  open OUT , ">$pos.vasp" ;

# Read in the input file
  while (<IN>){$file .= $_ ;}
  close IN ;  
  @file = split /\n/ , $file ;

# Check for new vasp5 format
  $line = $file[5];
  $line=~s/^\s+//;
  @line=split(/\s+/,$line);
#  if ($line=~/^s/i) {
  if ($line[0]=~ /^\d+$/) {
    $index = 0;
  } else {
    $index = 1;
  }

# Get the number of names of atomi types
  if ($index == 0){$line = $file[0] ;}
  else{$line = $file[5] ;}
  chomp($line) ; $line=~s/^\s+//g ; @line=split /\s+/,$line ;
  $types = join " " , @line ;
  $nel = @line ;
  $line = $file[5+$index] ;
  chomp($line) ; $line=~s/^\s+//g ; @line=split /\s+/,$line ;
  $numbers = join "  " , @line ;
# Calculate the total number of atoms
  while($line[$k] != undef){$natoms+=$line[$k++] ;}
 
# Write out the .vasp file
  print OUT "SOME INFORMATION WE DONT CARE ABOUT   ","NCLASS=",$nel,"  ATOM=",$types,"\n" ;
  print OUT "   ",$numbers,"\n" ;
  print OUT "Direct","\n" ;
  print OUT "   ","\n" ;
  for ($i=1; $i<5; $i++){print OUT $file[$i],"\n" ;}
  print OUT "   ","\n" ;
  $sh = 6+$index ;
  for ($i=1;$i<=$natoms; $i++){
    $line = $file[$sh+$i] ; chomp($line) ; $line=~s/^\s+//g ; @line=split /\s+/,$line ;
    printf OUT "%13.8f %11.8f %11.8f %5s\n",@line[0..2],"#$i" ;
  }

  close OUT ;

