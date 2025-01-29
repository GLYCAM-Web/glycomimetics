#ifndef GLYCAN_HPP 
#define GLYCAN_HPP

#include "includes/Glycan/oligosaccharide.hpp"
#include "includes/Glycan/monosaccharide.hpp"
#include "includes/MolecularModeling/atom.hpp"
#include "utility.hpp"
#include <math.h>

const double chi_coeff = 1.0;
const double chi_cutoff = 0.0;

void DetectAnomericCarbonAndOxygen(Glycan::Monosaccharide* mono, MolecularModeling::Atom*& ring_oxygen, MolecularModeling::Atom*& ano_carbon, MolecularModeling::Atom*& ano_oxygen, AtomVector& atoms_replaced){
	AtomVector& ca = mono->cycle_atoms_;
	//AtomVector& sa = mono->side_atoms_;

	for (unsigned int i = 0; i < ca.size(); i++){
		MolecularModeling::Atom* a = ca[i];
		std::string el = a->GetElementSymbol();
		if (el != "C") continue;

		MolecularModeling::Atom* ring_o = NULL, *ano_o = NULL;
		int NOS_count = 0, cycle_NOS_count = 0, side_NOS_count = 0;
		AtomVector n = a->GetNode()->GetNodeNeighbors();

		for (unsigned int j = 0; j < n.size(); j++){
			MolecularModeling::Atom* na = n[j];
			if (std::find(atoms_replaced.begin(), atoms_replaced.end(), na) != atoms_replaced.end()) continue;
			std::string nel = na->GetElementSymbol();

			if (nel != "N" && nel != "O" && nel != "S") continue;
			NOS_count++;

			if (std::find(ca.begin(), ca.end(), na) != ca.end()){
				ring_o = na;
 				cycle_NOS_count++;
			}
			else{
			 	side_NOS_count++;
				ano_o = na;
			}
		}
		if (NOS_count == 2 && cycle_NOS_count == 1 && side_NOS_count == 1){
			ring_oxygen = ring_o;
			ano_carbon = a;
			ano_oxygen = ano_o;
		}
	}
	return;
}

//Borrowed from Anita Vina-Carb source code
double phi_alpha_energy(double phi_angle)
{
	double LH = 2.97696467271672, Lc = -199.494365163839, LW = 677.808323900125, RH = 102.253303636096, Rc = 170.599580473404, RW = 1696.78443699429, aH = 10.7448005875571, ac = -105.313553566706, aW = 4724.58364072706, bH = 3.67344580413578, bc = 6.20116671874232, bW = 1347.72056251564, cH = 2.06094652655659, cc = 91.6553021324274, cW = 1500.02002601097, Off = 1.00501e-30, dH = 6.19388683252667, dc = -22.9786969888816, dW = 2122.27783139301, eH = -2.11153017593601, ec = 83.6019123356148, eW = 1254.13371108961, fH = -98.0013005657107, fc = 170.012289132741, fW = 1598.73272567307, Leftx, Rightx, ax, bx, cx, dx, ex, fx, x, Totx;
	x=phi_angle;
	Leftx = LH * exp(-pow((x-(Lc)),2.0)/LW);
	Rightx = RH * exp(-pow((x-(Rc)),2.0)/RW);
	ax = aH * exp(-pow((x-(ac)),2.0)/aW);
	bx = bH * exp(-pow((x-(bc)),2.0)/bW);
	cx = cH * exp(-pow((x-(cc)),2.0)/cW);
	dx = dH * exp(-pow((x-(dc)),2.0)/dW);
	ex = eH * exp(-pow((x-(ec)),2.0)/eW);
	fx = fH * exp(-pow((x-(fc)),2.0)/fW);
	Totx = Rightx + Leftx + ax + bx + cx + dx + ex + fx;
	return Totx;
}

double phi_beta_energy(double phi_angle)
{
	double Lc = -330.769995527134, aH = 5.93533323829663, ac = -152.080139620062, aW = 6049.77220005964, bH = 22.467372096061, bc = -23.5159916173247, bW = 606.89715970453, cH = 10.0360057033439, cc = 120.962836525241, cW = 4037.89330459304, LH = 450.540038600828, LW = 4449.7622241787, RH = 23.7118506901333, Rc = 304.624980492529, RW = 8375.1929028027, /*Off = -2.27829251796721,*/ Off = -2.1283, dH = -18.1406478565247, dc = -24.2677756921736, dW = 543.050986049266, eH = 5.88226333077368, ec = 19.6321032903376, eW = 897.92664572344, Leftx, Rightx, ax, bx, cx, dx, ex, x, Totx;
	x=phi_angle;
	Leftx = LH * exp(-pow((x-(Lc)),2.0)/LW);
	Rightx = RH * exp(-pow((x-(Rc)),2.0)/RW);
	ax = aH * exp(-pow((x-(ac)),2.0)/aW);
	bx = bH * exp(-pow((x-(bc)),2.0)/bW);
	cx = cH * exp(-pow((x-(cc)),2.0)/cW);
	dx = dH * exp(-pow((x-(dc)),2.0)/dW);
	ex = eH * exp(-pow((x-(ec)),2.0)/eW);
	Totx = Rightx + Leftx + ax + bx + cx + dx + ex + Off;
	return Totx;
}

double ScoreAnomericPhiTorsion(Glycan::Monosaccharide* mono, AtomVector& phi, int coord_index){
	//return 0.00;
	if (phi.empty()) return 0.00;
	double angle_degrees = GetDihedral(phi[0], phi[1], phi[2], phi[3], coord_index);

	Glycan::SugarName& sn = mono->sugar_name_;
	std::string& iso = sn.isomer_;
	std::string& conf = sn.configuration_;

	if (iso == "L") angle_degrees *= -1.00;
	double chi_energy = (conf == "a") ? phi_alpha_energy(angle_degrees) : phi_beta_energy(angle_degrees);

	//chi_cutoff applies to each individual torsional energy.
	if (chi_energy <= chi_cutoff) chi_energy = 0.00;
	//chi_coeff_applies to the sum of all individual energy values. If there's only one it doesn't matter.
	//But in the future if there are multiple, move the chi_coeff scaling to a place to apply to the total chi energy. 
	chi_energy *= chi_coeff;
	//std::cout << angle_degrees << "-" << chi_energy << std::endl;
	return chi_energy;
}
#endif //GLYCAN_HPP
