How to execute this example:
../Bin/gm.exe -f sample_input_file.txt

Detailed explanation:
1. ../Bin/gm.exe is the C++ program that builds putative glycomimetic (GM) compounds, peforms GA searching, predicts initial binding pose and energy, and writes output files. 

2. "sample_input_file.txt" is the input file required for the gm program. File content and explantion:
ComplexPdb:4ljh_chainA.pdbqt
	(1) Path to the protein-carbohydrate PDBQT file

OpenValence:203_A_O1_HO1-moieties-O1_R.pdbqt
	(1) Attach R groups on the following atom: residue 203, chain A, name O1.
	(2) Remove Hydroxyl atom HO1 and all downstream atoms. Then replace with derivative moiety.
	(3) Path to moiety library is "moieties"
	(4) Only select the moiety "O1_R.pdbqt", which technically is a glob pattern. 
		If you want all moieties in a directory, simply say "-pdbqt". This way all pdbqt moiety files will be globbed.

Interval:30
	(1) Deprecated. Will remove in future versions. Leave intact for now

NumThreads:4
	(1) Deprecated. Will remove in future versions. Leave intact for now

OutputPath:output
	(1) Path to the directory where output files will be written

LogFile:sample.log
	(1) Path to log file.

3. Output files:
	(1) Four output files are written for each glycomimetic compound:
		--R_group_receptor.pdb: A PDB file containing the receptor.
		--R_group_ligand.pdb: A PDB file containing the GM ligand.
		--R_group_pdb2glycam.log: Specifies whether an atom has been assigned GLYCAM parameters. If so, type and charge values. 
		--When "R_group" = "natural", files correspond to the endogenous protein-carbohydrate complex from the input PDBQT. 
		--complex_R_group.pdb: A PDB file containing both receptor and GM. For visualization only, not used for MD simlation.

	(2) "complex_best_each_position.pdb": When multiple moieties/open valence positions are analyzed, this PDB contains the best solution.

4. Expected output files are contained in the "expected_output" directory.

