#ifndef GLYCAN_HPP 
#define GLYCAN_HPP

#include "includes/Glycan/oligosaccharide.hpp"
#include "includes/Glycan/monosaccharide.hpp"
#include "includes/MolecularModeling/atom.hpp"
#include "utility.hpp"

void DetectAnomericCarbonAndOxygen(Glycan::Monosaccharide* mono, MolecularModeling::Atom*& ano_carbon, MolecularModeling::Atom*& ano_oxygen){
	AtomVector& ca = mono->cycle_atoms_;
	//AtomVector& sa = mono->side_atoms_;

	for (unsigned int i = 0; i < ca.size(); i++){
		MolecularModeling::Atom* a = ca[i];
		std::string el = a->GetElementSymbol();
		if (el != "C") continue;

		MolecularModeling::Atom* ano_o = NULL;
		int NO_count = 0, cycle_NO_count = 0, side_NO_count = 0;
		AtomVector n = a->GetNode()->GetNodeNeighbors();

		for (unsigned int j = 0; j < n.size(); j++){
			MolecularModeling::Atom* na = n[j];
			std::string nel = na->GetElementSymbol();

			if (nel != "N" && nel != "O") continue;
			NO_count++;

			if (std::find(ca.begin(), ca.end(), na) != ca.end()){
 				cycle_NO_count++;
			}
			else{
			 	side_NO_count++;
				ano_o = na;
			}
		}

		if (NO_count == 2 && cycle_NO_count == 1 && side_NO_count == 1){
			ano_carbon = a;
			ano_oxygen = ano_o;
		}
	}
	return;
}

AtomVector DetectPhiAnomericTorsion(Glycan::Monosaccharide* mono, MolecularModeling::Atom* anomeric_carbon, MolecularModeling::Atom* anomeric_oxygen, MolecularModeling::Atom* moiety_head_atom){
	AtomVector anomeric_phi_torsion (4, NULL);
	AtomVector& ca = mono->cycle_atoms_;

	AtomVector cn = anomeric_carbon->GetNode()->GetNodeNeighbors();
	for (unsigned int i = 0; i < cn.size(); i++){
		MolecularModeling::Atom* a = cn[i];
		if (std::find(ca.begin(), ca.end(), a) == ca.end()) continue;
		std::string el = a->GetElementSymbol();
		if (el == "O" || el == "N") anomeric_phi_torsion[0] = a;
	}

	anomeric_phi_torsion[1] = anomeric_carbon;
	anomeric_phi_torsion[2] = anomeric_oxygen;
	anomeric_phi_torsion[3] = moiety_head_atom;

	return anomeric_phi_torsion;
}

double ScoreAnomericPhiTorsion(Glycan::Monosaccharide* mono, AtomVector& phi, int thread_index){
	double angle_degrees = GetDihedral(phi[0], phi[1], phi[2], phi[3], thread_index);
	return 0.00; //TODO: Get the math expression form, based on alpha/beta, D/L, etc.  
}
#endif //GLYCAN_HPP
