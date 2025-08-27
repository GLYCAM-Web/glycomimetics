#../Bin/gm.exe -f sample_input_file.txt

cd output

sed '/CONECT/d' 4cst_O1_R_ligand.pdb > 4cst_ligand_noconect.pdb
sed '/CONECT/d' natural_ligand.pdb > natural_ligand_noconect.pdb

net_charge=0
../../Bin/gm2md.exe ../non_gm_input/4cst.mol2 4cst_ligand_noconect.pdb 4cst_O1_R_pdb2glycam.log ../non_gm_input/gaff.dat ../non_gm_input/4cst.frcmod ../non_gm_input/4cst_O1_R_resp_charges.out 4cst_glycam_gaff.frcmod 4cst_glycam_gaff.off ${net_charge}

../../Bin/gm2md.exe ../non_gm_input/natural.mol2 natural_ligand_noconect.pdb natural_pdb2glycam.log ../non_gm_input/gaff.dat ../non_gm_input/natural.frcmod ../non_gm_input/natural_resp_charges.out natural_glycam_gaff.frcmod natural_glycam_gaff.off ${net_charge}
