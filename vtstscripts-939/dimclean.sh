#!/bin/sh

  echo "   ****   "
  echo " Cleaning up after a dimer run "
  echo "   ****   "

  Bin=`dirname "$0"`
  ZIP=${VTST_ZIP-gzip}

  cp INCAR KPOINTS POSCAR CONTCAR $1
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

  echo -e " In directory: `pwd`"
  rm -f DOSCAR EIGENVAL IBZKPT OSZICAR PCDAT CHG

  if [ -s WAVECAR ] 
  then
    mv WAVECAR $1 ; "$ZIP" $1/WAVECAR &
  else
    rm -f WAVECAR
  fi
  if [ -s CHGCAR ]
  then
    mv CHGCAR $1; "$ZIP" $1/CHGCAR &
  else
    rm -f CHGCAR
  fi

  mv DIMCAR $1
 
  if [ -s XDATCAR ]
  then
    "$Bin/xdat2xyz.pl" > /dev/null;
     mv XDATCAR $1 ; "$ZIP" $1/XDATCAR &
  fi

  if [ -s movie.xyz ]
  then
     mv movie.xyz $1
  fi
  if [ -s OUTCAR ]
  then
     mv OUTCAR $1 ; "$ZIP" $1/OUTCAR &
  fi
  if [ -s MODECAR ]
  then
     cp MODECAR $1
  fi
  if [ -s NEWMODECAR ]
  then
     cp NEWMODECAR MODECAR; mv NEWMODECAR $1;
  fi

  "$Bin/pos2con.pl" POSCAR > /dev/null
  "$Bin/con2xyz.pl" POSCAR.con > /dev/null
  mv POSCAR.con POSCAR.xyz $1 ;

  "$Bin/pos2con.pl" CONTCAR > /dev/null
  "$Bin/con2xyz.pl" CONTCAR.con > /dev/null
  mv CONTCAR.con CONTCAR.xyz $1 ;

  if [ -s CENTCAR ]
  then
     cp CENTCAR $1
     mv CENTCAR POSCAR
     rm -f CONTCAR
  else
     mv CONTCAR POSCAR
  fi

