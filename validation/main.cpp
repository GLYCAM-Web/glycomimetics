
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

    //Valid = have sugars and available open valence positions.
    bool is_valid = false, pdb2glycam_available = false;
    std::vector<std::string> comments;

    std::vector<Glycan::Monosaccharide*> monos= std::vector<Glycan::Monosaccharide*>();
    std::vector<Glycan::Oligosaccharide*> oligos = assemblyA.ExtractSugars(amino_libs,monos,false,false);

	std::vector<Glycan::Oligosaccharide*> oligos_unlinked;
    std::vector<Glycan::Monosaccharide*> monos_unlinked;

    for (unsigned int i = 0; i < oligos.size(); i++){
        Glycan::Oligosaccharide* this_oligo = oligos[i];
        std::string condensed_sequence = this_oligo->IUPAC_name_;
        std::string aglycone = condensed_sequence.substr(condensed_sequence.find_last_of("-")+1);

        //If the aglycone has a protein residue name, skip
        if( std::find( gmml::PROTEINS, ( gmml::PROTEINS + gmml::PROTEINSSIZE ), aglycone ) != ( gmml::PROTEINS + gmml::PROTEINSSIZE ) ){
            std::cout << " Oligo " << " '" << condensed_sequence << "'" << " is N or O linked to protein. Skipping.\n";
            continue;
        }

        oligos_unlinked.push_back(this_oligo);
        std::vector<Glycan::Monosaccharide*>& this_oligo_monos = oligos[i]->mono_nodes_;
        monos_unlinked.insert(monos_unlinked.end(), this_oligo_monos.begin(), this_oligo_monos.end());
    }

    if (oligos_unlinked.empty()){
		is_valid = false;
        comments.push_back("Not elegible for pdb2glycam because no sugars were detected");
    }

    if (!pdb2glycam_successful){
        comments.push_back("Pdb2glycam matching failed. Cannot use this feature");
		pdb2glycam_available = false;
    }
    else{
        pdb2glycam_available = true;
    }

    //Detect available atoms for derivatization
    std::vector<available_atom> available_atoms = FindOpenValenceAtoms(monos_unlinked);  

	ResidueVector residues = assemblyA.GetResidues();
	std::map<MolecularModeling::Residue*, int> rearranged_residues_map;
	std::map<MolecularModeling::Atom*, int> rearranged_atoms_map;
	RearrangeResiduesAndAtoms(residues, monos_unlinked, rearranged_residues_map, rearranged_atoms_map);

	//std::string pdb_atomnum_, glycam_resname_, glycam_resnum_, glycam_atomnum_;
	for (unsigned int i = 0; i < available_atoms.size(); i++){
		available_atom& aa = available_atoms[i];
		MolecularModeling::Atom* a = aa.atom_;
		MolecularModeling::Residue* r = a->GetResidue();

		int glycam_resum = rearranged_residues_map[r];
		int glycam_atomnum = rearranged_atoms_map[a];
		aa.glycam_resnum_ = std::to_string(glycam_resum);
		aa.glycam_atomnum_ = std::to_string(glycam_atomnum);
	}
	
    std::string output_file_path_str = std::string(argv[2]);
	std::ofstream output_file(output_file_path_str);
	if (output_file.fail()){
		std::cout << "Failed to create " << output_file_path_str << " for writing." << std::endl;
		std::exit(1);
	}

    int oligo_index = 0;
	std::cout << "Oligos size: " << oligos_unlinked.size() << " \n";
	for (unsigned int i = 0; i < oligos_unlinked.size(); i++){
        Glycan::Oligosaccharide* oligo = oligos_unlinked[i];
		output_file << "Oligosaccharide " << oligo_index + 1 << " condensed sequence: " << oligo->IUPAC_name_ << "\n";
		oligo_index++;
    }

	std::cout << "Num avail atoms: " << available_atoms.size() << std::endl;
    for (unsigned int i = 0; i < available_atoms.size(); i++){
        available_atom& atom = available_atoms[i];
		//std::string residue_index_str_, atom_name_, atom_to_replace_;
		std::cout << "Open for derivatization: " << atom.residue_index_ << "-" << atom.atom_name_ << "-" << atom.atom_to_replace_ << std::endl;
		atom.print_attribute(output_file);
    }
	output_file.close();

    if (available_atoms.empty()){
        comments.push_back("No available positions for modification detected. For now must be ring hydroxyl/amino groups");
		is_valid = false;
    }
    else{
        is_valid = true;
    }
 
    response this_response(is_valid, pdb2glycam_available, available_atoms, comments);
	return 0;
}
