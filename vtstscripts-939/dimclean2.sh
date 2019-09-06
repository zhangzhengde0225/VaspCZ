#!/bin/sh

  echo "   ****   "
  echo " Cleaning up after a dimer run "
  echo "   ****   "

  Bin=`dirname "$0"`
  ZIP=${VTST_ZIP-gzip}

  cp INCAR KPOINTS $1
  rm -f vasprun.xml
  if [ ${VTST_STDOUT} ]
  then
    if [ -s ${VTST_STDOUT} ]
    then
      mv ${VTST_STDOUT} $1
    fi
  fi
  if [ ${VTST_STDERR} ]
  then
    rm -f ${VTST_STDERR}
  fi

  mkdir $1/01
  mkdir $1/02

# 1st image
  cd 01
  echo -e " In directory: `pwd`"
  rm -f DOSCAR EIGENVAL IBZKPT OSZICAR PCDAT CHG

  if [ -s WAVECAR ] 
  then
    mv WAVECAR ../$1/01 ; "$ZIP" $1/01/WAVECAR &
  else
    rm -f WAVECAR
  fi
  if [ -s CHGCAR ]
  then
    mv CHGCAR ../$1/01 ; "$ZIP" ../$1/01/CHGCAR &
  else
    rm -f CHGCAR
  fi

  cp out.dat dimer.dat ../$1
  rm -f out.dat dimer.dat

  if [ -s XDATCAR ]
  then
     mv XDATCAR ../$1/01 ; "$ZIP" ../$1/01/XDATCAR &
  fi
  if [ -s OUTCAR ]
  then
     mv OUTCAR ../$1/01 ; "$ZIP" ../$1/01/OUTCAR &
  fi


  "$Bin/pos2con.pl" POSCAR > /dev/null
  "$Bin/con2xyz.pl" POSCAR.con > /dev/null
  mv POSCAR.con POSCAR.xyz ../$1/01 ;
  "$Bin/pos2con.pl" CONTCAR > /dev/null
  "$Bin/con2xyz.pl" CONTCAR.con > /dev/null
  "$Bin/pos2con.pl" CONTCAR.con > /dev/null
  mv CONTCAR.con CONTCAR.xyz ../$1/01 ;
  cp CONTCAR POSCAR ../$1/01

  mv CONTCAR POSCAR

# 2nd image
  cd ../02
  echo -e " In directory: `pwd`"
  rm -f DOSCAR EIGENVAL IBZKPT OSZICAR PCDAT CHG out.dat dimer.dat stdout

  if [ -s WAVECAR ]
  then
    mv WAVECAR ../$1/02 ; "$ZIP" ../$1/02/WAVECAR &
  else
    rm -f WAVECAR
  fi
  if [ -s CHGCAR ]
  then
    mv CHGCAR ../$1/02 ; "$ZIP" ../$1/02/CHGCAR &
  else
    rm -f CHGCAR
  fi

  if [ -s XDATCAR ]
  then
     mv XDATCAR ../$1/02 ; "$ZIP" ../$1/02/XDATCAR &
  fi
  if [ -s OUTCAR ]
  then
     mv OUTCAR ../$1/02 ; "$ZIP" ../$1/02/OUTCAR &
  fi


  "$Bin/pos2con.pl" POSCAR > /dev/null
  "$Bin/con2xyz.pl" POSCAR.con > /dev/null
  mv POSCAR.con POSCAR.xyz ../$1/02 ;
  "$Bin/pos2con.pl" CONTCAR > /dev/null
  "$Bin/con2xyz.pl" CONTCAR.con > /dev/null
  "$Bin/pos2con.pl" CONTCAR.con > /dev/null
  mv CONTCAR.con CONTCAR.xyz ../$1/02 ;
  cp CONTCAR POSCAR ../$1/02

  mv CONTCAR POSCAR

