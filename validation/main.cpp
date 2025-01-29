
#include "includes/gmml.hpp"
#include "includes/MolecularModeling/assembly.hpp"
#include "includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp"
#include "includes/ParameterSet/OffFileSpace/offfile.hpp"
#include "includes/ParameterSet/OffFileSpace/offfileresidue.hpp"
#include "includes/ParameterSet/OffFileSpace/offfileprocessingexception.hpp"
#include "includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp"
#include "includes/InputSet/PdbFileSpace/pdbfile.hpp"
#include "includes/InputSet/PdbFileSpace/pdbremarksection.hpp"
#include "includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp"
#include "includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp"
#include "includes/InputSet/PdbqtFileSpace/pdbqtremarkcard.hpp"
#include "includes/utils.hpp"
#include "src/MolecularMetadata/guesses.cc"

#include "../src/vina_bond_by_distance_for_pdb.hpp"
#include "../src/pdb2glycam.hpp"
#include "validation.hpp"

//#include "boost/tokenizer.hpp"
#include <vector>
#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <pthread.h>
#include <iterator>
#include <sstream>
//#include <boost/filesystem.hpp>
#include <functional> //std::greater

typedef std::vector<MolecularModeling::Atom*> AtomVector;

int main(int argc, char* argv[]){
    std::string file_path_str = std::string(argv[1]);
    MolecularModeling::Assembly assemblyA(file_path_str, gmml::InputFileType::PDB); 
    VinaBondByDistanceForPDB(assemblyA, 0);

	char* gemshome = std::getenv("GEMSHOME");
	if (!gemshome){
        std::cout << "GEMSHOME environment variable must be set. Aborting." << std::endl;
        return 0;
    }

    std::string gems_home(gemshome);
	std::string lib1 = gems_home + "/gmml/dat/CurrentParams/leaprc.ff12SB_2014-04-24/amino12.lib";
	std::string lib2 = gems_home + "/gmml/dat/CurrentParams/leaprc.ff12SB_2014-04-24/aminoct12.lib";
	std::string lib3 = gems_home + "/gmml/dat/CurrentParams/leaprc.ff12SB_2014-04-24/aminont12.lib";
    std::vector<std::string> amino_libs = {lib1, lib2, lib3};
    std::string prep = gems_home + "/gmml/dat/prep/GLYCAM_06j-1.prep";

    //Attempt pdb2glycam matching
    std::map<MolecularModeling::Atom*, MolecularModeling::Atom*> actual_template_atom_match;
    AtomVector atoms = assemblyA.GetAllAtomsOfAssembly();
    bool pdb2glycam_successful = pdb2glycam_matching(file_path_str, actual_template_atom_match, atoms, gmml::InputFileType::PDB, amino_libs, prep);

    response r;
	r.pdb2glycam_successful_ = pdb2glycam_successful;
	std::vector<Glycan::Monosaccharide*> monos;
	std::vector<Glycan::Oligosaccharide*> oligos;

	if (pdb2glycam_successful){
    	oligos = assemblyA.ExtractSugars(amino_libs,monos,false,false);
	}

	std::vector<Glycan::Oligosaccharide*> oligos_valid;
    std::vector<Glycan::Monosaccharide*> monos_valid;

    for (unsigned int i = 0; i < oligos.size(); i++){
        Glycan::Oligosaccharide* this_oligo = oligos[i];
        std::string condensed_sequence = this_oligo->IUPAC_name_;
        std::string aglycone = condensed_sequence.substr(condensed_sequence.find_last_of("-")+1);

        //If the aglycone has a protein residue name, skip
        if( std::find( gmml::PROTEINS, ( gmml::PROTEINS + gmml::PROTEINSSIZE ), aglycone ) != ( gmml::PROTEINS + gmml::PROTEINSSIZE ) ){
            std::cout << "Oligo " << " '" << condensed_sequence << "'" << " is N or O linked to protein. Skipping.\n";
            continue;
        }

        oligos_valid.push_back(this_oligo);
        std::vector<Glycan::Monosaccharide*>& this_oligo_monos = oligos[i]->mono_nodes_;
        monos_valid.insert(monos_valid.end(), this_oligo_monos.begin(), this_oligo_monos.end());
    }
	r.oligos_valid_ = oligos_valid;

	ResidueVector residues = assemblyA.GetResidues();
	std::map<MolecularModeling::Residue*, int> rearranged_residues_map;
	std::map<MolecularModeling::Atom*, int> rearranged_atoms_map;
	RearrangeResiduesAndAtoms(residues, monos_valid, rearranged_residues_map, rearranged_atoms_map);

    //Detect available atoms for derivatization
	for (unsigned int i = 0; i < oligos_valid.size(); i++){
        Glycan::Oligosaccharide* oligo = oligos_valid[i];
		//std::cout << "Valid oligosaccharide " << i + 1 << " condensed sequence: " << oligo->IUPAC_name_ << "\n";
		std::vector<Glycan::Monosaccharide*>& this_oligo_monos = oligo->mono_nodes_;		
    	std::vector<available_atom> available_atoms = FindOpenValenceAtoms(this_oligo_monos);  

		for (unsigned int j = 0; j < available_atoms.size(); j++){
			available_atom& aa = available_atoms[j];
			MolecularModeling::Atom* a = aa.atom_;
			MolecularModeling::Residue* r = a->GetResidue();

			int glycam_resum = rearranged_residues_map[r];
			int glycam_atomnum = rearranged_atoms_map[a];
			aa.glycam_resnum_ = std::to_string(glycam_resum);
			aa.glycam_atomnum_ = std::to_string(glycam_atomnum);
		}

		r.available_atoms_.push_back(available_atoms);
    }
	
    std::string output_file_path_str = std::string(argv[2]);
	std::ofstream output_file(output_file_path_str);
	if (output_file.fail()){
		std::cout << "Failed to create " << output_file_path_str << " for writing." << std::endl;
		std::exit(1);
	}

	r.write(output_file);
	output_file.close();
	return 0;
}
