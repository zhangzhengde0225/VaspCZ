#!/bin/sh

  echo "   ****   "
  echo " Cleaning up after a neb run "
  echo "   ****   "

#  DIRS=$(ls -d ??)
  DIRS=$(ls -d [0-9][0-9])
  ZIP=${VTST_ZIP-gzip}

  c=0
  for i in $DIRS
  do
    c=$((c+1))
    mkdir $1/$i
  done

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

  j=0
  for i in $DIRS
  do
    j=$((j+1))
    cd $i
    echo -e " In directory: `pwd`"

    rm -f *.xyz > /dev/null
    rm -f *.vasp > /dev/null
    rm -f *.dat *.eps > /dev/null
    cp POSCAR OUTCAR ../$1/$i
    "$ZIP" ../$1/$i/OUTCAR &

    case $j in
    1)
      ;;
    $c)
      ;;
#    [2-$((c-1))])
    *)
      rm -f DOSCAR EIGENVAL IBZKPT OSZICAR PCDAT CHG OUTCAR
      if [ -f stdout ]
      then 
        rm -f stdout 
      fi
      mv XDATCAR ../$1/$i ; "$ZIP" ../$1/$i/XDATCAR &
      cp POSCAR CONTCAR ../$1/$i
      mv CONTCAR POSCAR

      if [ -s WAVECAR ] 
      then
        mv WAVECAR ../$1/$i ; "$ZIP" ../$1/$i/WAVECAR &
      else
        rm -f WAVECAR
      fi
      if [ -s CHGCAR ] 
      then
         mv CHGCAR ../$1/$i ; "$ZIP" ../$1/$i/CHGCAR &
      else
        rm -f CHGCAR
      fi
      ;;
    esac
    cd ../
  done
