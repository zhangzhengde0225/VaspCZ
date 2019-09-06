#!/bin/sh
set -x
set -e

echo "   ****   "
echo " Cleaning up after a normal vasp run "
echo "   ****   " 

VTSTPATH=`dirname "$0"`
ZIP=${VTST_ZIP-gzip}
DEST="$1"

# Create xyz movie from XDATCAR.
"$VTSTPATH/xdat2xyz.pl" > /dev/null
if [ -s movie.xyz ]
then
   mv movie.xyz $DEST
fi

# Create energy and forces plot.
"$VTSTPATH/vef.pl" > /dev/null
if [ -s fe.dat ]
then
   mv fe.dat $DEST
fi
if [ -s vaspout.??? ]
then
   mv vaspout.??? $DEST
fi

# Convert POSCAR and CONTCAR to con and xyz formats.
for file in POSCAR CONTCAR
do
    "$VTSTPATH/pos2con.pl" $file > /dev/null
    "$VTSTPATH/con2xyz.pl" $file.con > /dev/null 
    mv $file.con $file.xyz $DEST
done

# Files to delete.
rm -f EIGENVAL IBZKPT PCDAT vasprun.xml CHG

# Files to copy.
cp POSCAR CONTCAR INCAR KPOINTS $DEST

# Files to move and compress.
for file in WAVECAR CHGCAR AECCAR? PROCAR DOSCAR XDATCAR OUTCAR REPORT \
            ICONST OSZICAR
do
    # If a file exists and has data in it, then keep it.
    if [ -s $file ]; then
        mv $file $DEST

        # Compress it in the background.
        "$ZIP" $DEST/$file &

    # If the file is zero length (or doesn't exist) delete it.
    else
        rm -f $file
    fi
done

# If $VTST_STDOUT is set and the file is of non-zero length, move it to $DEST.
if [ ${VTST_STDOUT} ] 
then
  if [ -s ${VTST_STDOUT} ]
  then
    mv ${VTST_STDOUT} $DEST
  fi
fi

# If $VTST_STDERR is set and exists, delete it.
if [ ${VTST_STDERR} ]
then
  rm -f ${VTST_STDERR}
fi

# Move CONTCAR to POSCAR to prepare for next run.
mv CONTCAR POSCAR
