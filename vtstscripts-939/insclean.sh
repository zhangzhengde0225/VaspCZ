#!/bin/sh

# Cleaning up after an Instanton VASP run

  DIRS=$(ls -d [0-9][0-9])

  for i in $DIRS
  do
    cd $i
    pwd
    rm CHG* D* E* IB* OS* PC*
    if [ -e stdout ] 
    then
      rm stdout
    fi
    gzip OUTCAR XDATCAR
    if [ -s WAVECAR ]
    then
      gzip WAVECAR
    else
      rm -f WAVECAR
    fi
    cd ..
  done

  
  cp -r ${DIRS} INCAR KPOINTS $1
  mv ${VTST_STDOUT} vasprun.xml $1
  rm PI*

  for i in $DIRS
  do
    cd $i
    mv CONTCAR POSCAR
    mv NEWMODECAR MODECAR
    rm OUT* XDAT* WAVE* ins*
    cd ..
  done
