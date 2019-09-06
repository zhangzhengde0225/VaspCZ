#!/bin/sh

# this script strips md run information from OUTCAR file
# Original from Joel Kress; modified by dano

# usage: mddata.sh OUTCARfile 

  echo "   ****   "
  echo " Striping MD run information from OUTCAR "
  echo "   ****   " 

#temp.press, first column contains the pressure for each timestep
#temp.energy, first column contains the internal energy for each timestep
#temp.temp contains the temperature for each timestep
#for velocity rescaling, temperature should be constant
#for Nose-Hoover, temperature should fluctuate about the target T 
  Bin=`dirname "$0"`

  if [ -e tmp1 ]
  then 
    rm -f tmp1
  fi
  if [ -e tmp2 ]
  then 
    rm -f tmp2
  fi
  if [ -e md.press ]
  then 
    rm -f md.press
  fi
  echo "#      Pressure (kB)   Pullay Stress " > md.press
  grep  'external pressure = ' $1 > tmp1
  sed 's/external pressure =//' tmp1 > tmp2
  sed 's/kB  Pullay stress =//' tmp2 > tmp1
  sed 's/kB//' tmp1 >> md.press
  rm tmp1 tmp2

  if [ -e md.energy ]
  then 
    rm -f md.energy
  fi
  echo "#  Internal Energy    sigma->0  (eV)" > md.energy
  grep  'energy  wi' $1 > tmp1
  sed 's/  energy  without entropy=//' tmp1 > tmp2
  sed 's/ energy(sigma->0) = //' tmp2 >> md.energy
  rm tmp1 tmp2

  if [ -e md.temp ]
  then 
    rm -f md.temp
  fi
  echo "#kinetic Energy   temperature (K)" > md.temp
  grep 'EKIN' $1 > tmp1
  sed 's/  kinetic Energy EKIN   =        //' tmp1 > tmp2
  sed 's/(temperature//' tmp2 > tmp1
  sed 's/K)//' tmp1 >> md.temp
  rm tmp1 tmp2

  if [ -e md.totale ]
  then 
    rm -f md.totale
  fi
  echo "#    Total Energy (eV)" > md.totale
  grep '\% ion-electron   TOTEN  =' $1 > tmp1
  sed 's/\% ion-electron   TOTEN  =//' tmp1 > tmp2
  sed 's/see above//' tmp2 >> md.totale
  rm tmp1 tmp2

  if [ -e md.econs ]
  then 
    rm -f md.econs
  fi
  echo "# Configuration Energy (eV)" > md.econs
  grep '  total energy   ETOTAL = ' $1 > tmp1
  sed 's/ total energy   ETOTAL = //' tmp1 > tmp2
  sed 's/ eV//' tmp2 >> md.econs
  rm tmp1 tmp2

  if [ -e md.sigmaxx ]
  then 
    rm -f md.sigmaxx
  fi
  if [ -e md.sigmayy ]
  then 
    rm -f md.sigmayy
  fi
  if [ -e md.sigmazz ]
  then 
    rm -f md.sigmazz
  fi
  #echo "#  sigmaxx (eV)" > md.sigmaxx
  #echo "#  sigmayy (eV)" > md.sigmayy
  echo "#  sigmazz (eV)" > md.sigmazz
  grep -A 3 'stress matrix ' $1 > tmp1
  grep -B 1 -e -- tmp1 | cut -c 30-40 > tmp2
  sed '/^$/d' tmp2 >> md.sigmazz
  rm tmp1 tmp2

  echo "   ****   "
  echo " Plotting total energy and pressure "
  echo "   ****   " 

  gnuplot "$Bin/mdplot.gnu"  > /dev/null
