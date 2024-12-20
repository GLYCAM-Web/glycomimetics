
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

//#include "boost/tokenizer.hpp"
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <pthread.h>
#include <iterator>
#include <sstream>
//#include <boost/filesystem.hpp>
#include <functional> //std::greater

struct available_atom{
    available_atom(std::string residue_index_str, std::string resname, std::string chain_id, std::string atom_name, std::string atom_to_replace){
    	this->residue_index_str_ = residue_index_str;
		this->resname_ = resname;
		this->chain_id_ = chain_id;
		this->atom_name_ = atom_name;
		this->atom_to_replace_ = atom_to_replace;
    }
	void print_attribute(std::ofstream& output){
		//std::cout << "Available open valence atom: " << std::endl;
		//std::cout << "Residue index: " << this->residue_index_str_ << std::endl;
		//std::cout << "Atom to modify: " << this->atom_name_ << std::endl;
		//std::cout << "Downstream atom to replace: " << this->atom_to_replace_ << std::endl;
		//std::cout << "New option: " << this->residue_index_str_ << "-" << this->atom_name_ << "-" << this->atom_to_replace_ << std::endl;
		output << this->residue_index_str_ << "-" << this->resname_ << "-" << this->chain_id_ << "-" << this->atom_name_ << "-" << this->atom_to_replace_ << std::endl;
	}
    std::string residue_index_str_, resname_, chain_id_, atom_name_, atom_to_replace_;
};

struct response{
    response(bool valid, bool pdb2glycam_available, std::vector<available_atom>& available_atoms, std::vector<std::string>& comments){
        this->is_valid_ = valid;
		this->pdb2glycam_available_ = pdb2glycam_available;
		this->available_atoms_ = available_atoms; 
		this->comments_ = comments;
    }
    bool is_valid_ = false;
    bool pdb2glycam_available_ = false;
    std::vector<available_atom> available_atoms_;
    std::vector<std::string> comments_;
};

bool CheckAttachmentQualification(AtomVector& cycle_atoms, AtomVector& sa_arm, AtomVector& visited_atoms, MolecularModeling::Atom* current_atom){
	if (std::find(cycle_atoms.begin(), cycle_atoms.end(), current_atom) != cycle_atoms.end()) return false;
	if (std::find(sa_arm.begin(), sa_arm.end(), current_atom) == sa_arm.end()) return false;
	visited_atoms.push_back(current_atom);
	AtomVector neighbors = current_atom->GetNode()->GetNodeNeighbors();

	bool all_childs_true = true;
	for (unsigned int i = 0; i < neighbors.size(); i++){
		MolecularModeling::Atom* n = neighbors[i];
		if (std::find(visited_atoms.begin(), visited_atoms.end(), n) != visited_atoms.end()) continue;
		if (std::find(cycle_atoms.begin(), cycle_atoms.end(), n) != cycle_atoms.end()) continue;
		if (std::find(sa_arm.begin(), sa_arm.end(), n) == sa_arm.end()) continue;

		bool child_true = CheckAttachmentQualification(cycle_atoms, sa_arm, visited_atoms, n);

		if (!child_true) all_childs_true = false;
	}
	std::cout << "Current atom: " << current_atom->GetId() << " return " << all_childs_true << std::endl;

	return all_childs_true;
}

std::vector<std::pair<MolecularModeling::Atom*, AtomVector> > FindAttachReplacePairs(AtomVector& cycle_atoms, std::vector<AtomVector>& side_atoms){
	std::vector<std::pair<MolecularModeling::Atom*, AtomVector> > attach_replace_pairs;
	for (unsigned int i = 0; i < side_atoms.size(); i++){
		AtomVector& sa_arm = side_atoms[i];

		for (unsigned int j = 0; j < sa_arm.size(); j++){
			MolecularModeling::Atom* sa = sa_arm[j];
			if (sa == NULL) continue;
			std::cout << "\nThis sa: " << sa->GetId() << std::endl;
			AtomVector neighbors = sa->GetNode()->GetNodeNeighbors();

			AtomVector possible_replacements;
			for (unsigned int k = 0; k < neighbors.size(); k++){
				MolecularModeling::Atom* n = neighbors[k];
				AtomVector visited_atoms(1, sa);

				bool qualify = CheckAttachmentQualification(cycle_atoms, sa_arm, visited_atoms, n);
				if (qualify){
					possible_replacements.push_back(n);		
				}
			}

			int num_replacements = possible_replacements.size();
			if (num_replacements == 0){
				std::cout << "Cannot find replacement atom for " << sa->GetId() << "\n";
			}
			else if (num_replacements >= 2){
				std::cout << "Multiple replacement atoms for " << sa->GetId() << "\n";
				for (unsigned int x = 0; x < possible_replacements.size(); x++){
					std::cout << "It is: " << possible_replacements[x]->GetId() << std::endl;
				}

			}
			else{
				std::cout << "Replacement for " << sa->GetId() << " is " << possible_replacements[0]->GetId() << "\n";
				attach_replace_pairs.emplace_back(std::make_pair(sa, possible_replacements));
			}
		}

	}

	return attach_replace_pairs;
}

std::vector<available_atom> FindOpenValenceAtoms(std::vector<Glycan::Monosaccharide*>& monos){
	std::vector<available_atom> available_atoms;

    for (unsigned int i = 0; i < monos.size(); i++){
		std::cout << "find1" << std::endl;
        Glycan::Monosaccharide* mono = monos[i];

		AtomVector& cycle_atoms = mono->cycle_atoms_;
		MolecularModeling::Residue* this_residue = cycle_atoms[0]->GetResidue();

		std::cout << "find2" << std::endl;
		mono->InitiateDetectionOfCompleteSideGroupAtoms();
		std::cout << "find3" << std::endl;
		std::vector<AtomVector>& side_atoms = mono->side_atoms_;
		std::vector<std::pair<MolecularModeling::Atom*, AtomVector> > attach_replace_pairs = FindAttachReplacePairs(cycle_atoms, side_atoms);
		std::cout << "find4" << std::endl;

		std::string condensed_name = mono->sugar_name_.monosaccharide_short_name_;
		//std::cout << "Mono condensed name: " << condensed_name << " , " << std::endl;

		std::string residue_id = this_residue->GetId();
        std::vector<std::string> underscore_split_token = gmml::Split(residue_id, "_");
        std::string pdb_residue_index = underscore_split_token[2];
		std::string chain_id = this_residue->GetChainID();
		std::string pdb_resname = this_residue->GetName(); 
		AtomVector this_residue_atoms = this_residue->GetAtoms();

		//Now: writing logic to select possible open valence atoms.
		/*for (unsigned int j = 0; j < cycle_atoms.size(); j++){
	    	MolecularModeling::Atom* cycle_atom = cycle_atoms[j];
	    	AtomVector cycle_neighbors = cycle_atom->GetNode()->GetNodeNeighbors();

	    	for (unsigned int k = 0; k < cycle_neighbors.size(); k++){
	        	MolecularModeling::Atom* neighbor = cycle_neighbors[k];
				std::string neighbor_element = neighbor->GetElementSymbol();

				if (std::find(cycle_atoms.begin(), cycle_atoms.end(), neighbor) != cycle_atoms.end()) continue;
				if (std::find(this_residue_atoms.begin(), this_residue_atoms.end(), neighbor) == this_residue_atoms.end()) continue; 
				if (neighbor_element != "N" && neighbor_element != "O") continue;
				
		    	available_atoms.emplace_back(available_atom(residue_index, resname, chain_id, cycle_atom->GetName(), neighbor->GetName()));
	    	}
		}*/
    }
	std::exit(1);
	return available_atoms;
}

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
	std::cout << "cp-1\n";
    bool pdb2glycam_successful = pdb2glycam_matching(file_path_str, actual_template_atom_match, atoms, gmml::InputFileType::PDB, amino_libs, prep);
	std::cout << "cp-00" << std::endl;

    //Valid = have sugars and available open valence positions.
    bool is_valid = false, pdb2glycam_available = false;
    std::vector<std::string> comments;

    std::vector<Glycan::Monosaccharide*> monos= std::vector<Glycan::Monosaccharide*>();
    std::vector<Glycan::Oligosaccharide*> oligos = assemblyA.ExtractSugars(amino_libs,monos,false,false);

    if (oligos.empty()){
		is_valid = false;
        comments.push_back("Not elegible for pdb2glycam because no sugars were detected");
    }

	std::cout << "cp0" << std::endl;

    if (!pdb2glycam_successful){
        comments.push_back("Pdb2glycam matching failed. Cannot use this feature");
		pdb2glycam_available = false;
    }
    else{
        pdb2glycam_available = true;
    }

	std::cout << "cp1" << std::endl;
    //Detect available atoms for derivatization
    std::vector<available_atom> available_atoms = FindOpenValenceAtoms(monos);  
	std::cout << "cp2" << std::endl;
	
    //Obtain oligosaccharide condensed sequence
    //std::vector<Glycan::Oligosaccharide*> oligos = assemblyA.ExtractSugars(amino_libs,monos,false,false);
    std::map<Glycan::Oligosaccharide*, std::vector<std::string> > oligo_id_map;
    std::map<Glycan::Oligosaccharide*, std::vector<MolecularModeling::Residue*> > oligo_residue_map;

    for (std::vector<Glycan::Oligosaccharide*>::iterator oligo_it = oligos.begin(); oligo_it != oligos.end(); oligo_it++){
            std::vector<std::string> empty_vector = std::vector<std::string>();
            oligo_id_map[*oligo_it] = empty_vector;
            std::vector<MolecularModeling::Residue*> empty_residue_vector = std::vector<MolecularModeling::Residue*>();
            oligo_residue_map[*oligo_it] = empty_residue_vector;
    }


    std::string output_file_path_str = std::string(argv[2]);
	std::ofstream output_file(output_file_path_str);
	if (output_file.fail()){
		std::cout << "Failed to create " << output_file_path_str << " for writing." << std::endl;
		std::exit(1);
	}

    int oligo_index = 0;
    for (std::map<Glycan::Oligosaccharide*, std::vector<MolecularModeling::Residue*> >::iterator mapit = oligo_residue_map.begin(); mapit != oligo_residue_map.end(); mapit++){
        Glycan::Oligosaccharide* oligo = mapit->first;
	output_file << "Oligosaccharide " << oligo_index + 1 << " condensed sequence: " << oligo->IUPAC_name_ << "\n";
	oligo_index++;
    }

	std::cout << "Num avail atoms: " << available_atoms.size() << std::endl;
    for (unsigned int i = 0; i < available_atoms.size(); i++){
        available_atom& atom = available_atoms[i];
		//std::string residue_index_str_, atom_name_, atom_to_replace_;
		std::cout << "Open for derivatization: " << atom.residue_index_str_ << "-" << atom.atom_name_ << "-" << atom.atom_to_replace_ << std::endl;
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
}
