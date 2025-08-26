cd frames_ligand
for i in *;do
    cd ${i}
    echo ${i}
    rm -f *renamed.pdb

    ../../rename_ligand.out ../../cocomplex_pdbs/${i}_cocomplex.pdb
    mv renamed.pdb ../../cocomplex_pdbs/${i}_cocomplex_renamed.pdb
    /home/yao/programs_src/MGLTools-1.5.4/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r ../../cocomplex_pdbs/${i}_cocomplex_renamed.pdb -A None -U None -o ../../cocomplex_pdbs/${i}_cocomplex_renamed.pdbqt

    for j in *.pdb;do
	framenum=$(echo ${j} | sed 's/.pdb//')
        #../../rename_ligand.out ${j} ../../${i}/${j} ../../cocomplex_pdbs/${j}_cocomplex.pdb
        ../../rename_ligand.out ${j}
	mv renamed.pdb ${framenum}_renamed.pdb
	../../rename_ligand.out ../../frames/${i}/${j}
	mv renamed.pdb ../../frames/${i}/${framenum}_renamed.pdb
    done
    cd ..
done
