#!/bin/sh

  echo "   ****   "
  echo " Cleaning up after a dynmat run "
  echo "   ****   "

  Bin=`dirname "$0"`
  ZIP=${VTST_ZIP-gzip}

  rm -f DOSCAR EIGENVAL IBZKPT OSZICAR PCDAT vasprun.xml CHG

  if [ -s WAVECAR ]
  then
    mv WAVECAR $1 ; "$ZIP" $1/WAVECAR &
  else
    rm -f WAVECAR
  fi
  if [ -s CHGCAR ]
  then
    mv CHGCAR $1 ; "$ZIP" $1/CHGCAR &
  else
    rm -f CHGCAR
  fi

  if [ -s XDATCAR ]
  then
     mv XDATCAR $1 ; "$ZIP" $1/XDATCAR &
  fi
  if [ -s OUTCAR ]
  then
     mv OUTCAR $1 ; "$ZIP" $1/OUTCAR &
  fi

  "$Bin/pos2con.pl" POSCAR > /dev/null
  "$Bin/con2xyz.pl" POSCAR.con > /dev/null
  mv POSCAR.con POSCAR.xyz $1 ;
  cp POSCAR INCAR KPOINTS $1
  if [ -s CONTCAR ]
  then
     mv CONTCAR $1 ;
  fi
  if [ -s DISPLACECAR ]
  then
     cp DISPLACECAR $1 ;
  fi
  if [ -s eigs.dat ]
  then
     mv eigs.dat $1 ;
  fi
  if [ -s modes.dat ]
  then
     mv modes.dat $1 ;
  fi
  if [ -s freq.mat ]
  then
     mv freq.mat $1 ;
  fi
  if [ -s freq.dat ]
  then
     mv freq.dat $1 ;
  fi

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
