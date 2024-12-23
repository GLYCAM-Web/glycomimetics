#ifndef VALIDATION_HPP
#define VALIDATION_HPP

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
//#include "src/MolecularMetadata/guesses.cc"

//#include "../src/vina_bond_by_distance_for_pdb.hpp"
//#include "../src/pdb2glycam.hpp"

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

struct available_atom{
    available_atom(MolecularModeling::Atom* atom, std::string pdb_resname, std::string chain_id, std::string residue_index, std::string atom_name, std::string atom_to_replace){
		this->atom_ = atom;
		this->pdb_resname_ = pdb_resname;
		this->chain_id_ = chain_id;
    	this->residue_index_ = residue_index;
		this->atom_name_ = atom_name;
		this->atom_to_replace_ = atom_to_replace;
    }
	void print_attribute(std::ofstream& output){
		output << this->pdb_resname_ << "-" << this->chain_id_ << "-" << this->residue_index_ << "-" << this->atom_name_ << "-" << this->pdb_atomnum_ << "-" << this->glycam_resname_ << "-" << glycam_resnum_ << "-" << this->glycam_atomnum_ << "\n";
	}
	MolecularModeling::Atom* atom_ = NULL;
    std::string pdb_resname_, chain_id_, residue_index_, atom_name_, atom_to_replace_;
	std::string pdb_atomnum_, glycam_resname_, glycam_resnum_, glycam_atomnum_;
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

	return all_childs_true;
}

std::vector<std::pair<MolecularModeling::Atom*, AtomVector> > FindAttachReplacePairs(AtomVector& cycle_atoms, std::vector<AtomVector>& side_atoms){
	std::vector<std::pair<MolecularModeling::Atom*, AtomVector> > attach_replace_pairs;
	for (unsigned int i = 0; i < side_atoms.size(); i++){
		AtomVector& sa_arm = side_atoms[i];

		for (unsigned int j = 0; j < sa_arm.size(); j++){
			MolecularModeling::Atom* sa = sa_arm[j];
			if (sa == NULL) continue;

			AtomVector neighbors = sa->GetNode()->GetNodeNeighbors();
			std::string sa_element = sa->GetElementSymbol();

			if (sa_element != "O" && sa_element != "N") continue;
			if (sa_element == "O"){
				MolecularModeling::Atom* h = NULL;
				for (unsigned int k = 0; k < neighbors.size(); k++){
					MolecularModeling::Atom* n = neighbors[k];
					if (std::find(cycle_atoms.begin(), cycle_atoms.end(), n) != cycle_atoms.end()) continue;
        			if (std::find(sa_arm.begin(), sa_arm.end(), n) == sa_arm.end()) continue;

					std::string n_element = n->GetElementSymbol();
					if (n_element == "H"){
						h = n;
						break;
					}
				}

				if (h != NULL){
					AtomVector replacements(1, h);
					attach_replace_pairs.emplace_back(std::make_pair(sa, replacements));
				}
			} 
			else if (sa_element == "N"){
				MolecularModeling::Atom* heavy = NULL;
				int num_qualify = 0;

				for (unsigned int k = 0; k < neighbors.size(); k++){
                    MolecularModeling::Atom* n = neighbors[k];
					if (std::find(cycle_atoms.begin(), cycle_atoms.end(), n) != cycle_atoms.end()) continue;
        			if (std::find(sa_arm.begin(), sa_arm.end(), n) == sa_arm.end()) continue;

                    std::string n_element = n->GetElementSymbol();
                    if (n_element == "H") continue;
					AtomVector visited_atoms(1, sa);
					bool qualify = CheckAttachmentQualification(cycle_atoms, sa_arm, visited_atoms, n);
					if (qualify){ 
						num_qualify++;
						heavy = n;
					}
                }

				if (num_qualify == 0){
					std::cout << "Cannot find replacement atom for " << sa->GetId() << "\n";
				}
				else if (num_qualify >= 2){
                	std::cout << "Multiple replacement atoms for " << sa->GetId() << "\n";
            	}
				else{
					AtomVector replacements(1, heavy);
                	//std::cout << "Replacement for " << sa->GetId() << " is " << heavy->GetId() << "\n";
					attach_replace_pairs.emplace_back(std::make_pair(sa, replacements));
            	}
			}

		}

	}

	return attach_replace_pairs;
}

std::vector<available_atom> FindOpenValenceAtoms(std::vector<Glycan::Monosaccharide*>& monos){
	std::vector<available_atom> available_atoms;

    for (unsigned int i = 0; i < monos.size(); i++){
        Glycan::Monosaccharide* mono = monos[i];

		AtomVector& cycle_atoms = mono->cycle_atoms_;
		MolecularModeling::Residue* this_residue = cycle_atoms[0]->GetResidue();

		mono->InitiateDetectionOfCompleteSideGroupAtoms();
		std::vector<AtomVector>& side_atoms = mono->side_atoms_;
		std::vector<std::pair<MolecularModeling::Atom*, AtomVector> > attach_replace_pairs = FindAttachReplacePairs(cycle_atoms, side_atoms);

		std::string condensed_name = mono->sugar_name_.monosaccharide_short_name_;

		std::string residue_id = this_residue->GetId();
        std::vector<std::string> underscore_split_token = gmml::Split(residue_id, "_");

		std::string pdb_resname = this_residue->GetName(); 
		std::string chain_id = this_residue->GetChainID();
        std::string pdb_residue_index = underscore_split_token[2];
		
		for (unsigned int j = 0; j < attach_replace_pairs.size(); j++){
			std::pair<MolecularModeling::Atom*, AtomVector>& this_pair = attach_replace_pairs[j];
			MolecularModeling::Atom* open_atom = this_pair.first;
			AtomVector& replace_atoms = this_pair.second; 
			std::string atom_name = open_atom->GetName();
			std::string atom_to_replace = replace_atoms[0]->GetName();

			std::string atom_id = open_atom->GetId();
			std::vector<std::string> underscore_split_tokens = gmml::Split(atom_id, "_");
			int pdb_atomnum = std::stoi(underscore_split_tokens[1]) + 1; //Atom index now starts at zero?

			available_atoms.emplace_back(available_atom(open_atom, pdb_resname, chain_id, pdb_residue_index, atom_name, atom_to_replace));
			available_atoms.back().pdb_atomnum_ = std::to_string(pdb_atomnum);
			available_atoms.back().glycam_resname_ = condensed_name;
		}
    }
	//std::exit(1);
	return available_atoms;
}

void RearrangeResiduesAndAtoms(ResidueVector& residues, std::vector<Glycan::Monosaccharide*>& monos, std::map<MolecularModeling::Residue*, int>& rearranged_residues_map, std::map<MolecularModeling::Atom*, int>& rearranged_atoms_map){
	ResidueVector rearranged_residues, proteins, ions, waters, sugars;	

	for (unsigned int i = 0; i < monos.size(); i++){
		MolecularModeling::Residue* r = monos[i]->cycle_atoms_[0]->GetResidue();
		if (std::find(sugars.begin(), sugars.end(), r) == sugars.end()){
			sugars.push_back(r);
		}
	}

	std::vector<std::string> common_ions = {"cu", "zn", "cd", "mo", "mg", "k", "ca", "fe", "fe2", "co", "ni"};

	for (unsigned int i = 0; i < residues.size(); i++){
		MolecularModeling::Residue* r = residues[i];
		std::string resname_lower = r->GetName();
		for (char &c : resname_lower){
        	c = std::tolower(c);
    	}

		if (r->CheckIfProtein()){
			proteins.push_back(r);
		}
		else if (std::find(common_ions.begin(), common_ions.end(), resname_lower) != common_ions.end()){
			ions.push_back(r);
		}	
		else if (resname_lower == "hoh" || resname_lower == "wat"){
			//std::cout << "Water residue: " << r->GetId() << std::endl;
			waters.push_back(r);
		}
		else if (std::find(sugars.begin(), sugars.end(), r) == sugars.end()){
			std::cout << "Residue " << r->GetId() << " is not protein, ions, waters, or sugars. Ignored.\n";
        }
	}
	
	rearranged_residues.insert(rearranged_residues.end(), proteins.begin(), proteins.end());
	rearranged_residues.insert(rearranged_residues.end(), ions.begin(), ions.end());
	rearranged_residues.insert(rearranged_residues.end(), waters.begin(), waters.end());
	rearranged_residues.insert(rearranged_residues.end(), sugars.begin(), sugars.end());

	int atom_index = 0;
	for (unsigned int i = 0; i < rearranged_residues.size(); i++){
		MolecularModeling::Residue* r = rearranged_residues[i];
		rearranged_residues_map[r] = i+1;

		AtomVector atoms = r->GetAtoms();
		for (unsigned int j = 0; j < atoms.size(); j++){
			MolecularModeling::Atom* a = atoms[j];
			atom_index++;
			rearranged_atoms_map[a] = atom_index;
		}		
	}

	return;
}
#endif //VALIDATION_HPP
