#!/usr/bin/env perl
#;-*- Perl -*-

# Calculate the velocity from XDATCAR using forward finite difference

## 1-18-2012, update to new version of 
#### XDATCAR file from vasp5 only

# Open XDATCAR
$zip = $ENV{'VTST_ZIP'};
if($zip eq '') { $zip = 'gzip'; }

if(@ARGV == undef){
    if(-e "XDATCAR") {
        $xdatzipped = 0;
    } elsif(-e "XDATCAR.gz") {
        print " XDATCAR IS ZIPPED, UNZIPPING \n";
        system "gunzip XDATCAR.gz";
        $xdatzipped = 1;
    } elsif(-e "XDATCAR.bz2") {
        print " XDATCAR IS ZIPPED, UNZIPPING \n";
        system "bunzip2 XDATCAR.bz2";
        $xdatzipped = 1;
    } else {
        die " NO XDATCAR IN THIS DIRECTORY \n";
    }
    print " OPEN XDATCAR ... \n";
    open XDAT, "XDATCAR"
    or die " ...  XDATCAR MISSING FROM THIS DIRECTORY \n";
} else {
    print " OPEN ",@ARGV[0]," ... \n";
    open XDAT, @ARGV[0]
    or die " ... ",@ARGV[0]," MISSING FROM THIS DIRECTORY \n";
}

# Write velocity to the output file, "VDATCAR"
open VDAT, ">VDATCAR" ;
print VDAT $line ;

for($i=0; $i<=7; $i++) {
    $line = <XDAT>;
    printf VDAT $line; 

    chomp($line);
    $line =~ s/^\s+//g;
    @line = split /\s+/,$line;
    if($i == 1) { $vectscale = $line[0];}
    if($i == 2) { @cell1[0..2] = @line[0..2]; }
    if($i == 3) { @cell2[0..2] = @line[0..2]; }
    if($i == 4) { @cell3[0..2] = @line[0..2]; }
    if($i == 5) {
        if($line[0] != /\D/) { die "XDATCAR NOT IN VASP5 FORMAT\n"; }
        else {
            @elements = split /\s+/, $line;
            $ntypes = @elements;
        }
    }
    if($i == 6) {
        @neach[0..$ntypes-1] = @line[0..$ntypes-1];
            # Calculate the number of atoms
            while($neach[$k] != undef) {
                $natoms += $neach[$k++];
            }
        }
    }

# Get mass of each type and total mass 
if(-e "OUTCAR") { 
    $zipped = 0;
} elsif(-e "OUTCAR.gz") {
    print " OUTCAR IS ZIPPED, UNZIPPING \n";
    system "gunzip OUTCAR.gz";
    $zipped = 1;
} elsif(-e "OUTCAR.bz2") {
    print " OUTCAR IS ZIPPED, UNZIPPING \n";
    system "bunzip2 OUTCAR.bz2";
    $zipped = 1;
} else {
    die " NO OUTCAR IN THIS DIRECTORY \n";
}
$mline = `grep 'POMASS = ' OUTCAR`;
@mline = split /\s+/, $mline;
@mtype = splice @mline, 3;
$tm = 0;  #collect total mass
$index = -1;
$N_A = 6.0221412927e23;
$kb = 8.617332478e-5;
$evfactor = 1.60217733e-19;

for($i = 0; $i<$ntypes; $i++){
    for($j=0; $j<@neach[$i]; $j++){
        $index++;
        @m[$index] = @mtype[$i]; #store mass of each atom, kg/atom
        $tm += @m[$index]; #store total mass
        $masstmp = $m[$index];
    }
}

# Timestep in seconds (needs to be scaled)
$tline = `grep POTIM OUTCAR`;
chomp($tline);
$tline =~ s/^\s+//g;
@tline = split /\s+/, $tline;
$dt = @tline[2]; #time step in fs

# Read in the initial position from XDATCAR
for($i=0; $i<$natoms; $i++){
    $line = <XDAT>;
    chomp($line);
    $line =~ s/^\s+//g;
    @line = split /\s+/,$line;
    $Rxprev[$i] = @line[0];
    $Ryprev[$i] = @line[1];
    $Rzprev[$i] = @line[2];
#    print($Rxprev[$i]," ",$Ryprev[$i]," ",$Rzprev[$i],"\n");
}

# Read the rest of the file
  for(;;){
    $ke = 0;
    $line = <XDAT>; #get first line and enter loop
    if(!($line=<XDAT>)){ last; }
    for($i=0; $i<$natoms; $i++){
        chomp($line);
        $line =~ s/^\s+//g;
        @line = split /\s+/,$line;
        $Rxnext[$i] = @line[0];
        $Rynext[$i] = @line[1];
        $Rznext[$i] = @line[2];
        if($i < $natoms-1) { $line = <XDAT>; } #get next line, stay in loop
    }

# Calculate the velocity at the intermediate time t+dt
    for($i=0; $i<$natoms; $i++){
# Finite difference of two consecutive positions
        $dx = (@Rxnext[$i] - @Rxprev[$i]);
        $dy = (@Rynext[$i] - @Ryprev[$i]);
        $dz = (@Rznext[$i] - @Rzprev[$i]);
 # Apply periodic boundary conditions
      while($dx<-0.5){$dx+=1.0;} while($dx>0.5){$dx-=1.0;}
      while($dy<-0.5){$dy+=1.0;} while($dy>0.5){$dy-=1.0;}
      while($dz<-0.5){$dz+=1.0;} while($dz>0.5){$dz-=1.0;}
# Transform from direct to cartesian coordinates.
        $dxt = $dx*$cell1[0] + $dy*$cell2[0] + $dz*$cell3[0];
        $dyt = $dx*$cell1[1] + $dy*$cell2[1] + $dz*$cell3[1];
        $dzt = $dx*$cell1[2] + $dy*$cell2[2] + $dz*$cell3[2];
        @v[0] = $dxt*$vectscale/$dt;
        @v[1] = $dyt*$vectscale/$dt;
        @v[2] = $dzt*$vectscale/$dt; #A/fs
        $ke += 0.5*@m[$i]*(@v[0]*@v[0] + @v[1]*@v[1] + @v[2]*@v[2]);
        printf VDAT "%20.15f %20.15f %20.15f\n", @v[0], @v[1], @v[2];
    }
# Calculate total Kinetic Energy and Effective Temp
    $ke = $ke*10e6/$N_A/$evfactor;
    $teff = 2/$kb*$ke/(3*$natoms-3);
    printf VDAT "\nKINETIC ENERGY (eV): %20.15f\n", $ke;
    printf VDAT "TEMP EFF (K): %20.15f\n \n", $teff;

# Copy Rnext to Rprev
    for($i=0; $i<$natoms; $i++){
        @Rxprev[$i] = @Rxnext[$i];
        @Ryprev[$i] = @Rynext[$i];
        @Rzprev[$i] = @Rznext[$i];
    }
print VDAT "\n" ;
} #close (;;) loop
