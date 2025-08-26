#!/bin/bash
interval=1
if [ ! -z $1 ]; then
    echo "Interval argument supplied as $1"
    interval=$1
fi


for i in $(ls -d frames/*);do
    moiety_name=$(echo ${i} | sed 's/frames\///')
    echo "moiety name: ${moiety_name}"
    #dir=$(pwd)
    #echo $dir
    #sed -i '/Na+/d' *.pdbqt
    #sed -i '/Cl-/d' *.pdbqt
    #echo "./scoretraj.exe ./cocomplex_pdbs/${moiety_name}_cocomplex_renamed.pdbqt ${i} ${interval}"
    #cocomplex_pdbs/14_cocomplex_renamed.pdbqt
    ./scoretraj.exe ./cocomplex_pdbs/${moiety_name}_cocomplex_renamed.pdbqt ${i} ${interval}
    mv output.csv ./all_csvs/${moiety_name}.csv
done
