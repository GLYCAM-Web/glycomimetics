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

#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
typedef std::vector<MolecularModeling::Atom*> AtomVector;

void GenerateReceptorLigandRange(ResidueVector& complex_residues, std::pair<int, int>& rec_begin_end, std::pair<int, int>& lig_begin_end){
	rec_begin_end.first = 1;
	int num_protein_res = 0;

	std::vector<int> lig_indices;
	//std::vector<std::string> others = {"HOH", "WAT", "CA", "CL", "NA"};
	std::vector<std::string> others = {"CA"};
	std::vector<std::string> filter = {"HOH", "WAT"};

	for (unsigned int i = 0; i < complex_residues.size(); i++){
        MolecularModeling::Residue* r = complex_residues[i];
        int index = i + 1;

        std::string name = r->GetName();
		if (r->CheckIfProtein() || std::find(others.begin(), others.end(), name) != others.end()){
			num_protein_res++;
		}
		else if(std::find(filter.begin(), filter.end(), name) == filter.end()){
			lig_indices.push_back(index);
		}
    }
	rec_begin_end.second = num_protein_res;
	lig_begin_end.first = lig_indices[0];
	lig_begin_end.second = lig_indices.back();
	
}

std::pair<std::string, std::string> GenerateReceptorLigandMask(ResidueVector& complex_residues){
	std::stringstream rec;
	rec << ":";
	std::stringstream lig;
	lig << ":";
	std::vector<std::string> others = {"HOH", "WAT", "CA", "CL", "NA"};

	std::vector<int> rec_indices;
	std::vector<int> lig_indices;

	for (unsigned int i = 0; i < complex_residues.size(); i++){
		MolecularModeling::Residue* r = complex_residues[i];
		int index = i + 1;

		std::string name = r->GetName();
		if (r->CheckIfProtein() || std::find(others.begin(), others.end(), name) != others.end()){
			rec_indices.push_back(index);
		}
		else{
			lig_indices.push_back(index);
		}
	}

	std::vector<std::pair<int, int>> r_consecutive;
	int previous = rec_indices[0];
	for (unsigned int i = 0; i < rec_indices.size(); i++){
		int current = rec_indices[i]; 	
		if (i == 0){
			r_consecutive.push_back(std::make_pair(current, -1));
		}
		else if (current != previous + 1 || i == rec_indices.size() -1){
			std::pair<int, int>& this_segment = r_consecutive.back();
			this_segment.second = previous;
			r_consecutive.push_back(std::make_pair(current, -1));
		}
		previous = current;
	}

	std::vector<std::pair<int, int>> l_consecutive;
	previous = lig_indices[0];
    for (unsigned int i = 0; i < lig_indices.size(); i++){
        int current = lig_indices[i];
		std::cout << "Previous & current: " << previous << " & "<< current << std::endl;
        if (i == 0){
            l_consecutive.push_back(std::make_pair(current, -1));
        }
        else if (current != previous + 1 || i == lig_indices.size() -1){
            std::pair<int, int>& this_segment = l_consecutive.back();
            this_segment.second = previous;
            l_consecutive.push_back(std::make_pair(current, -1));
        }
        previous = current;
    }

	for (unsigned int i = 0; i < r_consecutive.size(); i++){
		std::pair<int, int>& this_segment = r_consecutive[i];
		int begin = this_segment.first;
		int end = this_segment.second;

		if (begin == end || end == -1){
			rec << begin << ",";
		}
		else{
			rec << begin << "-" << end << ",";
		}
	}

	for (unsigned int i = 0; i < l_consecutive.size(); i++){
        std::pair<int, int>& this_segment = l_consecutive[i];
        int begin = this_segment.first;
        int end = this_segment.second;

        if (begin == end || end == -1){
            lig << begin << ",";
        }
        else{
            lig << begin << "-" << end << ",";
        }
    }

	std::string rec_mask = rec.str();
	std::string lig_mask = lig.str();
	//Remove trailing comma ","
	rec_mask.erase(rec_mask.length()-1, 1);
	lig_mask.erase(lig_mask.length()-1, 1);

	return std::make_pair(rec_mask, lig_mask);
}

void WriteCalphaRestraintSection(std::ofstream& input_file, int num_receptor_residues){
	//return; //Temporary
	input_file << "Restrain all protein C-alpha atoms\n";
	input_file << "10.0\n";
	input_file << "FIND\n";
	input_file << "CA * * *\n";
	input_file << "SEARCH\n";
	input_file << "RES 1 ";
	input_file << num_receptor_residues;
	input_file << "\n";
	input_file << "END\n";
	return;
}

void WriteWaterRestraintSection(std::ofstream& input_file, std::vector<int>& water_oxygen_indices){
	//return; //Temporary
    input_file << "Restrain binding site water oxygens\n100.0" << std::endl;

    for (unsigned int i = 0; i < water_oxygen_indices.size(); i++){
        int& this_water_oxygen_index = water_oxygen_indices[i];
        input_file << "ATOM " << this_water_oxygen_index << " " << this_water_oxygen_index << std::endl;
    }

    input_file << "END" << std::endl;

    return;
}

void WriteStrandBreakRestraintSection(std::ofstream& input_file, std::vector<std::pair<int, int> >& strand_breaks){
    input_file << "Restrain strand breaks\n500.0" << std::endl;

    for (unsigned int i = 0; i < strand_breaks.size(); i++){
        std::pair<int, int>& this_break = strand_breaks[i];
        int& this_c_index = this_break.first;
        int& next_n_index = this_break.second;
        input_file << "ATOM " << this_c_index << " " << this_c_index << std::endl;
        input_file << "ATOM " << next_n_index << " " << next_n_index << std::endl;
    }

    input_file << "END" << std::endl;

    return;
}

void WriteTruncatedEndRestraintSection(std::ofstream& heating_input_file, std::pair<int, int>& head_n_tail_c_indices){
	//return; //Temporary
    int& head_n_index = head_n_tail_c_indices.first;
    int& tail_c_index = head_n_tail_c_indices.second;
    heating_input_file << "Restrain truncated end\n500.0" << std::endl;
    heating_input_file << "ATOM " << head_n_index << " " << head_n_index << std::endl;
    heating_input_file << "ATOM " << tail_c_index << " " << tail_c_index << std::endl;
    heating_input_file << "END" << std::endl;
}

std::pair<int, int> FindHeadNTailC(MolecularModeling::Assembly& receptor_pdb){
    std::pair<int, int> head_n_tail_c_indices;
    ResidueVector receptor_residues = receptor_pdb.GetResidues();
    ResidueVector protein_residues;
    for (unsigned int i = 0; i < receptor_residues.size(); i++){
        MolecularModeling::Residue* residue = receptor_residues[i];
        if (residue->CheckIfProtein()){
            protein_residues.push_back(residue);
        }
    }

    AtomVector all_receptor_atoms = receptor_pdb.GetAllAtomsOfAssembly();
    AtomVector first_protein_residue_atoms = protein_residues[0]->GetAtoms();
    AtomVector last_protein_residue_atoms = protein_residues.back()->GetAtoms();

    for (unsigned int i = 0; i < first_protein_residue_atoms.size(); i++){
        MolecularModeling::Atom* atom = first_protein_residue_atoms[i];
        if (atom->GetName() == "N"){
            AtomVector::iterator it = std::find(all_receptor_atoms.begin(), all_receptor_atoms.end(), atom);
            int head_n_index = std::distance(all_receptor_atoms.begin(), it) + 1;
            head_n_tail_c_indices.first = head_n_index;
            break;
        }
    }

    for (unsigned int i = 0; i < last_protein_residue_atoms.size(); i++){
        MolecularModeling::Atom* atom = last_protein_residue_atoms[i];
        if (atom->GetName() == "C"){
            AtomVector::iterator it = std::find(all_receptor_atoms.begin(), all_receptor_atoms.end(), atom);
            int tail_c_index = std::distance(all_receptor_atoms.begin(), it) + 1;
            head_n_tail_c_indices.second = tail_c_index;
            break;
        }
    }

    return head_n_tail_c_indices;
}

void WriteStepTwoRelaxationInputFile(std::string step2_input_path, std::pair<int, int>& head_n_tail_c_indices, std::vector<std::pair<int, int>>& strand_breaks, std::vector<int>& water_oxygens, int num_receptor_residues){
    std::ofstream step2_input_file(step2_input_path);
    if (step2_input_file.fail()){
        std::cout << "Failed to open " << step2_input_path << " for writing." << std::endl;
        std::exit(1);
    }

    std::string control_section = R"(Dan Roe Step Two Relaxation
&cntrl
  ig = -1,
  imin = 0,
  irest = 0,
  ntx = 1,
  ntb = 1,
  cut = 8.0,
  ntc = 2,
  ntf = 2,
  temp0 = 300.0,
  ntt = 3,
  gamma_ln = 5,
  nstlim = 15000,
  dt = 0.001,
  ntwx = 5000,
  ntpr = 5000,
  ntwe = 5000,
  ioutfm = 1,
  iwrap = 0,
  ntr = 1,
  restraintmask = '!:WAT= & !@H=',
  restraint_wt = 5.0,
 &end)";
    step2_input_file << control_section << std::endl;

	WriteCalphaRestraintSection(step2_input_file, num_receptor_residues);
    WriteTruncatedEndRestraintSection(step2_input_file, head_n_tail_c_indices);

    if (!strand_breaks.empty()){
        WriteStrandBreakRestraintSection(step2_input_file, strand_breaks);
    }
    if (!water_oxygens.empty()){
        WriteWaterRestraintSection(step2_input_file, water_oxygens);
    }

    step2_input_file << "END" << std::endl;

    return;

}

void WriteStepSixRelaxationInputFile(std::string step6_input_path, std::pair<int, int>& head_n_tail_c_indices, std::vector<std::pair<int, int>>& strand_breaks, std::vector<int>& water_oxygens, int num_receptor_residues){
    std::ofstream step6_input_file(step6_input_path);
    if (step6_input_file.fail()){
        std::cout << "Failed to open " << step6_input_path << " for writing." << std::endl;
        std::exit(1);
    }

    std::string control_section = R"(Dan Roe Step Six Relaxation
 &cntrl
  ig = -1,
  imin = 0,
  irest = 0,
  ntx = 1,
  ntb = 2,
  ntp = 1,
  cut = 8.0,
  ntc = 2,
  ntf = 2,
  temp0 = 300.0,
  ntt = 3,
  gamma_ln = 5,
  nstlim = 2500000,
  dt = 0.001,
  ntwx = 1000,
  ntpr = 1000,
  ntwe = 1000,
  ioutfm = 1,
  iwrap = 0,
  ntr = 1,
  restraintmask = '!:WAT= & !@H=',
  restraint_wt = 1.0,
 &end)";
    step6_input_file << control_section << std::endl;

	WriteCalphaRestraintSection(step6_input_file, num_receptor_residues);
    WriteTruncatedEndRestraintSection(step6_input_file, head_n_tail_c_indices);

    if (!strand_breaks.empty()){
        WriteStrandBreakRestraintSection(step6_input_file, strand_breaks);
    }
    if (!water_oxygens.empty()){
        WriteWaterRestraintSection(step6_input_file, water_oxygens);
    }

    step6_input_file << "END" << std::endl;

    return;

}

void WriteStepSevenRelaxationInputFile(std::string step7_input_path, std::pair<int, int>& head_n_tail_c_indices, std::vector<std::pair<int, int>>& strand_breaks, std::vector<int>& water_oxygens, int num_receptor_residues){
    std::ofstream step7_input_file(step7_input_path);
    if (step7_input_file.fail()){
        std::cout << "Failed to open " << step7_input_path << " for writing." << std::endl;
        std::exit(1);
    }

    std::string control_section = R"(Dan Roe Step Seven Relaxation
 &cntrl
  imin = 0,
  irest = 0,
  ntx = 5,
  ntb = 2,
  ntp = 1,
  cut = 8.0,
  ntc = 2,
  ntf = 2,
  temp0 = 300.0,
  ntt = 3,
  gamma_ln = 5,
  nstlim = 2500000,
  dt = 0.001,
  ntwx = 1000,
  ntpr = 1000,
  ntwe = 1000,
  ioutfm = 1,
  iwrap = 0,
  ntr = 1,
  restraintmask = '!:WAT= & !@H=',
  restraint_wt = 0.5,
 &end)";
    step7_input_file << control_section << std::endl;

	WriteCalphaRestraintSection(step7_input_file, num_receptor_residues);
    WriteTruncatedEndRestraintSection(step7_input_file, head_n_tail_c_indices);

    if (!strand_breaks.empty()){
        WriteStrandBreakRestraintSection(step7_input_file, strand_breaks);
    }
    if (!water_oxygens.empty()){
        WriteWaterRestraintSection(step7_input_file, water_oxygens);
    }

    step7_input_file << "END" << std::endl;

    return;

}

void WriteStepEightRelaxationInputFile(std::string step8_input_path, std::pair<int, int>& head_n_tail_c_indices, std::vector<std::pair<int, int>>& strand_breaks, std::vector<int>& water_oxygens, int num_receptor_residues){
    std::ofstream step8_input_file(step8_input_path);
    if (step8_input_file.fail()){
        std::cout << "Failed to open " << step8_input_path << " for writing." << std::endl;
        std::exit(1);
    }

    std::string control_section = R"(Dan Roe Step Eight Relaxation
 &cntrl
  imin = 0,
  irest = 0,
  ntx = 5,
  ntb = 2,
  ntp = 1,
  cut = 8.0,
  ntc = 2,
  ntf = 2,
  temp0 = 300.0,
  ntt = 3,
  gamma_ln = 5,
  nstlim = 2500000,
  dt = 0.001,
  ntwx = 1000,
  ntpr = 1000,
  ntwe = 1000,
  ioutfm = 1,
  iwrap = 0,
  ntr = 1,
  restraintmask = '@CA,C,O,N | @%Cg,Os,Cy,Oy',
  restraint_wt = 0.5,
 &end)";
    step8_input_file << control_section << std::endl;

	WriteCalphaRestraintSection(step8_input_file, num_receptor_residues);
    WriteTruncatedEndRestraintSection(step8_input_file, head_n_tail_c_indices);

    if (!strand_breaks.empty()){
        WriteStrandBreakRestraintSection(step8_input_file, strand_breaks);
    }
    if (!water_oxygens.empty()){
        WriteWaterRestraintSection(step8_input_file, water_oxygens);
    }

    step8_input_file << "END" << std::endl;

    return;

}

void WriteStepNineRelaxationInputFile(std::string step9_input_path, std::pair<int, int>& head_n_tail_c_indices, std::vector<std::pair<int, int>>& strand_breaks, std::vector<int>& water_oxygens, int num_receptor_residues){
    std::ofstream step9_input_file(step9_input_path);
    if (step9_input_file.fail()){
        std::cout << "Failed to open " << step9_input_path << " for writing." << std::endl;
        std::exit(1);
    }

    std::string control_section = R"(Dan Roe Step Nine Relaxation
 &cntrl
  imin = 0,
  irest = 0,
  ntx = 5,
  ntb = 2,
  ntp = 1,
  cut = 8.0,
  ntc = 2,
  ntf = 2,
  temp0 = 300.0,
  ntt = 3,
  gamma_ln = 5,
  nstlim = 2500000,
  dt = 0.002,
  ntwx = 1000,
  ntpr = 1000,
  ntwe = 1000,
  ioutfm = 1,
  iwrap = 0,
  ntr = 1,
 &end)";
    step9_input_file << control_section << std::endl;

	WriteCalphaRestraintSection(step9_input_file, num_receptor_residues);
    WriteTruncatedEndRestraintSection(step9_input_file, head_n_tail_c_indices);

    if (!strand_breaks.empty()){
        WriteStrandBreakRestraintSection(step9_input_file, strand_breaks);
    }
    if (!water_oxygens.empty()){
        WriteWaterRestraintSection(step9_input_file, water_oxygens);
    }

    step9_input_file << "END" << std::endl;

    return;

}

void WriteStepTenProductionInputFile(std::string step10_input_path, int num_complex_atoms, std::pair<int, int>& head_n_tail_c_indices, std::vector<std::pair<int, int>>& strand_breaks, std::vector<int>& water_oxygens, int num_receptor_residues){
    std::ofstream step10_input_file(step10_input_path);
    if (step10_input_file.fail()){
        std::cout << "Failed to open " << step10_input_path << " for writing." << std::endl;
        std::exit(1);
    }

    std::string control_section = R"(Dan Roe Step Ten Relaxation
 &cntrl
  imin = 0,
  irest = 0,
  ntx = 5,
  ntb = 2,
  ntp = 1,
  cut = 8.0,
  ntc = 2,
  ntf = 2,
  temp0 = 300.0,
  ntt = 3,
  gamma_ln = 5,
  nstlim = 50000000,
  dt = 0.002,
  ntwx = 5000,
  ntpr = 50000,
  ntwe = 50000,
  ntwr = 500000,
  ntwprt = 0,
  ioutfm = 1,
  iwrap = 0,
  ntr = 1,)";

    step10_input_file << control_section << std::endl;
	step10_input_file << "ntwprt = " << num_complex_atoms << std::endl;
	step10_input_file << "&end" << std::endl;

	WriteCalphaRestraintSection(step10_input_file, num_receptor_residues);
    WriteTruncatedEndRestraintSection(step10_input_file, head_n_tail_c_indices);

    if (!strand_breaks.empty()){
        WriteStrandBreakRestraintSection(step10_input_file, strand_breaks);
    }
    if (!water_oxygens.empty()){
        WriteWaterRestraintSection(step10_input_file, water_oxygens);
    }

    step10_input_file << "END" << std::endl;

    return;

}

//std::pair<int, int> rec_begin_end;
//std::pair<int, int> lig_begin_end;

void WriteGBSAComplexInputFile(std::string gbsa_complex_input_file_path, int num_complex_residues, int num_receptor_residues, std::pair<int, int>& rec_begin_end, std::pair<int, int>& lig_begin_end){
    std::ofstream gbsa_complex_input_file(gbsa_complex_input_file_path);
    if (gbsa_complex_input_file.fail()){
        std::cout << "Failed to open " << gbsa_complex_input_file_path << " for writing." << std::endl;
        std::exit(1);
    }

    std::string constant_section = R"(File generated by C++
&cntrl
 dec_verbose=0, ntb=0,
 surften=0.0072, extdiel=80.0, ncyc=0,
 cut=999.0, gbsa=2, imin=5, idecomp=1,
 igb=2, intdiel=1
/)";

    gbsa_complex_input_file << constant_section << std::endl;

    gbsa_complex_input_file << "Residues considered as REC" << std::endl;
    //gbsa_complex_input_file << "RRES 1 " << num_receptor_residues << std::endl;
    gbsa_complex_input_file << "RRES " << rec_begin_end.first << " "<< rec_begin_end.second << std::endl;
    gbsa_complex_input_file << "END" << std::endl;

    gbsa_complex_input_file << "Residues considered as LIG" << std::endl;
    //gbsa_complex_input_file << "LRES " << num_receptor_residues + 1 << " " << num_complex_residues << std::endl;
    gbsa_complex_input_file << "LRES " << lig_begin_end.first << " " << lig_begin_end.second << std::endl;
    gbsa_complex_input_file << "END" << std::endl;

	int real_num_protein = rec_begin_end.second;
	int real_num_ligand = lig_begin_end.second - lig_begin_end.first + 1;
	int real_num_complex = real_num_protein + real_num_ligand;

    gbsa_complex_input_file << "Residues to print" << std::endl;
    gbsa_complex_input_file << "RES 1 " << real_num_complex << std::endl;
    gbsa_complex_input_file << "END" << std::endl;

    gbsa_complex_input_file << "END" << std::endl;

    gbsa_complex_input_file.close();
}

void WriteGBSALigandInputFile(std::string gbsa_ligand_input_file_path, int num_complex_residues, int num_receptor_residues, std::pair<int, int>& lig_begin_end){
    std::ofstream gbsa_ligand_input_file(gbsa_ligand_input_file_path);
    if (gbsa_ligand_input_file.fail()){
        std::cout << "Failed to open " << gbsa_ligand_input_file_path << " for writing." << std::endl;
        std::exit(1);
    }

    std::string constant_section = R"(File generated by C++
&cntrl
 nsnb=99999, dec_verbose=0, ntb=0,
 surften=0.0072, extdiel=80.0, ncyc=0,
 cut=999.0, gbsa=2, imin=5, idecomp=1,
 igb=2, intdiel=1,
/)";

    gbsa_ligand_input_file << constant_section << std::endl;

    gbsa_ligand_input_file << "Residues considered as LIG" << std::endl;
    //gbsa_ligand_input_file << "LRES 1 " << num_complex_residues - num_receptor_residues << std::endl;

	int end_index = lig_begin_end.second - lig_begin_end.first + 1;
	gbsa_ligand_input_file << "LRES 1 " << end_index << std::endl;
    gbsa_ligand_input_file << "END" << std::endl;

    gbsa_ligand_input_file << "Residues to print" << std::endl;
    //gbsa_ligand_input_file << "RES 1 " << num_complex_residues - num_receptor_residues << std::endl;
	gbsa_ligand_input_file << "RES 1 " << end_index << std::endl;
    gbsa_ligand_input_file << "END" << std::endl;

    gbsa_ligand_input_file << "END" << std::endl;

    gbsa_ligand_input_file.close();
}

void WriteGBSAReceptorInputFile(std::string gbsa_receptor_input_file_path, int num_receptor_residues, std::pair<int, int>& rec_begin_end){
    std::ofstream gbsa_receptor_input_file(gbsa_receptor_input_file_path);
    if (gbsa_receptor_input_file.fail()){
        std::cout << "Failed to open " << gbsa_receptor_input_file_path << " for writing." << std::endl;
        std::exit(1);
    }

    std::string constant_section = R"(File generated by C++
&cntrl
 nsnb=99999, dec_verbose=0, ntb=0,
 surften=0.0072, extdiel=80.0, ncyc=0,
 cut=999.0, gbsa=2, imin=5, idecomp=1,
 igb=2, intdiel=1,
/)";

    gbsa_receptor_input_file << constant_section << std::endl;

    gbsa_receptor_input_file << "Residues considered as REC" << std::endl;
    //gbsa_receptor_input_file << "RRES 1 " << num_receptor_residues << std::endl;
    gbsa_receptor_input_file << "RRES " << rec_begin_end.first << " "<< rec_begin_end.second << std::endl;
    gbsa_receptor_input_file << "END" << std::endl;

    gbsa_receptor_input_file << "Residues to print" << std::endl;
    //gbsa_receptor_input_file << "RES 1 " << num_receptor_residues << std::endl;
    gbsa_receptor_input_file << "RES " << rec_begin_end.first << " "<< rec_begin_end.second << std::endl;
    gbsa_receptor_input_file << "END" << std::endl;

    gbsa_receptor_input_file << "END" << std::endl;

    gbsa_receptor_input_file.close();
}

void WriteGBSAQuasiHarmonicInputFile(std::string gbsa_entropy_input_file_path, int num_complex_residues, int num_receptor_residues, std::pair<int, int>& rec_begin_end, std::pair<int, int>& lig_begin_end){
	std::ofstream gbsa_entropy_input_file(gbsa_entropy_input_file_path);
	if (gbsa_entropy_input_file.fail()){
        std::cout << "Failed to open " << gbsa_entropy_input_file_path << " for writing." << std::endl;
        std::exit(1);
    }

	gbsa_entropy_input_file << "#File generated by C++\n";
	gbsa_entropy_input_file << "trajin _MMPBSA_complex.mdcrd\n";
	gbsa_entropy_input_file << "reference _MMPBSA_avgcomplex.pdb\n";
	/*gbsa_entropy_input_file << "rms mass reference :1-" << num_complex_residues << "\n";
	gbsa_entropy_input_file << "matrix mwcovar name comp.matrix :1-" << num_complex_residues << "\n";
	gbsa_entropy_input_file << "analyze matrix comp.matrix out _MMPBSA_complex_entropy.out thermo reduce\n";
	gbsa_entropy_input_file << "rms mass reference :1-" << num_receptor_residues << "\n";
	gbsa_entropy_input_file << "matrix mwcovar name rec.matrix :1-" << num_receptor_residues << "\n";
	gbsa_entropy_input_file << "analyze matrix rec.matrix out _MMPBSA_receptor_entropy.out thermo reduce\n";
	gbsa_entropy_input_file << "rms mass reference :" << num_receptor_residues + 1 << "-" << num_complex_residues << "\n";
	gbsa_entropy_input_file << "matrix mwcovar name lig.matrix :" << num_receptor_residues + 1 << "-" << num_complex_residues << "\n";
	gbsa_entropy_input_file << "analyze matrix lig.matrix out _MMPBSA_ligand_entropy.out thermo reduce\n";*/

	int real_num_protein = rec_begin_end.second;
	int real_num_ligand = lig_begin_end.second - lig_begin_end.first + 1;
	int real_num_complex = real_num_protein + real_num_ligand;

	gbsa_entropy_input_file << "rms mass reference :1-" << real_num_complex << "\n";
    gbsa_entropy_input_file << "matrix mwcovar name comp.matrix :1-" << real_num_complex << "\n";
    gbsa_entropy_input_file << "analyze matrix comp.matrix out _MMPBSA_complex_entropy.out thermo reduce\n";
    gbsa_entropy_input_file << "rms mass reference :" << rec_begin_end.first << "-" << rec_begin_end.second << "\n";
    gbsa_entropy_input_file << "matrix mwcovar name rec.matrix :" << rec_begin_end.first << "-" << rec_begin_end.second << "\n";
    gbsa_entropy_input_file << "analyze matrix rec.matrix out _MMPBSA_receptor_entropy.out thermo reduce\n";
    gbsa_entropy_input_file << "rms mass reference :" << lig_begin_end.first << "-" << lig_begin_end.second << "\n";
    gbsa_entropy_input_file << "matrix mwcovar name lig.matrix :" << lig_begin_end.first << "-" << lig_begin_end.second << "\n";
    gbsa_entropy_input_file << "analyze matrix lig.matrix out _MMPBSA_ligand_entropy.out thermo reduce\n";

	gbsa_entropy_input_file.close();
}

void WriteGBSAGeneralInputFile(std::string gbsa_general_input_file_path){
    std::ofstream gbsa_general_input_file(gbsa_general_input_file_path);
    if (gbsa_general_input_file.fail()){
        std::cout << "Failed to open " << gbsa_general_input_file_path << " for writing." << std::endl;        
        std::exit(1);
    }
    std::string constant_part = R"(File generated by C++
&general
startframe=1, endframe=10000, interval=1,
verbose=2,entropy=0,
/
&gb
igb=2
/
&decomp
idecomp=1
dec_verbose=0
/)"; 

    gbsa_general_input_file << constant_part << std::endl;
    gbsa_general_input_file.close();
}

std::vector<int> DetectBindingSiteWaterOxygens(AtomVector& complex_atoms){
    std::vector<int> water_oxygen_indices;
    
    for (unsigned int i = 0; i < complex_atoms.size(); i++){
        MolecularModeling::Atom* this_atom = complex_atoms[i];
        std::string element = this_atom->GetElementSymbol();
        std::string residue_name = this_atom->GetResidue()->GetName();

        if ((residue_name == "HOH" || residue_name == "WAT") && element == "O"){
            water_oxygen_indices.push_back(i + 1);
        }
    }

    return water_oxygen_indices; 
}

std::vector<std::pair<int, int> > DetectStrandBreaks(ResidueVector& receptor_residues, AtomVector& receptor_atoms){
    std::vector<std::pair<int, int> > strand_breaks;
    int last_residue_index = receptor_residues.size() -1;

    for(unsigned int i = 0; i < last_residue_index; i++){
        MolecularModeling::Residue* this_residue = receptor_residues[i];
        MolecularModeling::Residue* next_residue = receptor_residues[i+1];
        AtomVector this_residue_atoms = this_residue->GetAtoms();
        AtomVector next_residue_atoms = next_residue->GetAtoms();

        for (unsigned int j = 0; j < this_residue_atoms.size(); j++){
            MolecularModeling::Atom* atom = this_residue_atoms[j];

            if (atom->GetName() == "C"){
                for (unsigned int k = 0; k < next_residue_atoms.size(); k++){
                    MolecularModeling::Atom* next_atom = next_residue_atoms[k];

                    if (next_atom->GetName() == "N"){
                        double distance = atom->GetCoordinate()->Distance(next_atom->GetCoordinate());

                        if (distance > 2.0){
                            AtomVector::iterator it1 = std::find(receptor_atoms.begin(), receptor_atoms.end(), atom);
                            int atom_index1 = std::distance(receptor_atoms.begin(), it1) + 1;

                            AtomVector::iterator it2 = std::find(receptor_atoms.begin(), receptor_atoms.end(), next_atom);
                            int atom_index2 = std::distance(receptor_atoms.begin(), it2) + 1;
                            strand_breaks.emplace_back(std::pair<int, int>(atom_index1, atom_index2));
                        }
                    }
                }    
            }
        }
    }
    return strand_breaks;
}

void WriteStepOneMinimizationInputFile(std::string step1_input_path, std::pair<int, int>& head_n_tail_c_indices, std::vector<std::pair<int, int>>& strand_breaks, std::vector<int>& water_oxygens){
    std::ofstream step1_input_file(step1_input_path);
    if (step1_input_file.fail()){
        std::cout << "Failed to open " << step1_input_path << " for writing." << std::endl;
        std::exit(1);
    }

    std::string constant_section = R"(File generated by C++: Dan Roe Step One
&cntrl
  imin = 1,
  maxcyc = 10000,
  ncyc = 5000,
  cut = 8.0,
  ntwr = 1000,
  ntpr = 1000,
  ntr = 1,
  restraintmask = '!:WAT= & !@H=',
  restraint_wt = 5.0,
 &end)";

    step1_input_file << constant_section << std::endl;

	WriteTruncatedEndRestraintSection(step1_input_file, head_n_tail_c_indices);

    if (!strand_breaks.empty()){
        WriteStrandBreakRestraintSection(step1_input_file, strand_breaks);
    }
    if (!water_oxygens.empty()){
        WriteWaterRestraintSection(step1_input_file, water_oxygens);
    }
	step1_input_file << "END" << std::endl;
    step1_input_file.close();

    return;
}

void WriteStepThreeMinimizationInputFile(std::string step3_input_path, std::pair<int, int>& head_n_tail_c_indices, std::vector<std::pair<int, int>>& strand_breaks, std::vector<int>& water_oxygens){
    std::ofstream step3_input_file(step3_input_path);
    if (step3_input_file.fail()){
        std::cout << "Failed to open " << step3_input_path << " for writing." << std::endl;
        std::exit(1);
    }

    std::string constant_section = R"(File generated by C++: Dan Roe Step Three
&cntrl
  imin = 1,
  maxcyc = 10000,
  ncyc = 5000,
  cut = 8.0,
  ntwr = 5000,
  ntpr = 5000,
  ntr = 1,
  restraintmask = '!:WAT= & !@H=',
  restraint_wt = 2.0,
 &end)";

    step3_input_file << constant_section << std::endl;

	WriteTruncatedEndRestraintSection(step3_input_file, head_n_tail_c_indices);

    if (!strand_breaks.empty()){
        WriteStrandBreakRestraintSection(step3_input_file, strand_breaks);
    }
    if (!water_oxygens.empty()){
        WriteWaterRestraintSection(step3_input_file, water_oxygens);
    }
	step3_input_file << "END" << std::endl;
    step3_input_file.close();

    return;
}

void WriteStepFourMinimizationInputFile(std::string step4_input_path, std::pair<int, int>& head_n_tail_c_indices, std::vector<std::pair<int, int>>& strand_breaks, std::vector<int>& water_oxygens){
    std::ofstream step4_input_file(step4_input_path);
    if (step4_input_file.fail()){
        std::cout << "Failed to open " << step4_input_path << " for writing." << std::endl;
        std::exit(1);
    }

    std::string constant_section = R"(File generated by C++: Dan Roe Step Four
 &cntrl
  imin = 1,
  maxcyc = 10000,
  ncyc = 5000,
  cut = 8.0,
  ntwr = 5000,
  ntpr = 5000,
  ntr = 1,
  restraintmask = '!:WAT= & !@H=',
  restraint_wt = 0.1,
 &end)";

    step4_input_file << constant_section << std::endl;

    WriteTruncatedEndRestraintSection(step4_input_file, head_n_tail_c_indices);

    if (!strand_breaks.empty()){
        WriteStrandBreakRestraintSection(step4_input_file, strand_breaks);
    }
    if (!water_oxygens.empty()){
        WriteWaterRestraintSection(step4_input_file, water_oxygens);
    }
	step4_input_file << "END" << std::endl;
    step4_input_file.close();

    return;
}

void WriteStepFiveMinimizationInputFile(std::string step5_input_path, std::pair<int, int>& head_n_tail_c_indices, std::vector<std::pair<int, int>>& strand_breaks, std::vector<int>& water_oxygens){
    std::ofstream step5_input_file(step5_input_path);
    if (step5_input_file.fail()){
        std::cout << "Failed to open " << step5_input_path << " for writing." << std::endl;
        std::exit(1);
    }

    std::string constant_section = R"(File generated by C++: Dan Roe Step Five
&cntrl
  imin = 1,
  maxcyc = 10000,
  ncyc = 5000,
  cut = 8.0,
  ntwr = 5000,
  ntpr = 5000,
  ntr = 1,
 &end)";

    step5_input_file << constant_section << std::endl;

    WriteTruncatedEndRestraintSection(step5_input_file, head_n_tail_c_indices);

    if (!strand_breaks.empty()){
        WriteStrandBreakRestraintSection(step5_input_file, strand_breaks);
    }
    if (!water_oxygens.empty()){
        WriteWaterRestraintSection(step5_input_file, water_oxygens);
    }
	step5_input_file << "END" << std::endl;
    step5_input_file.close();

    return;
}

//argv:1, tleap_receptor_pdb; 2, tleap_cocomplex_pdb ; 3, step 1 input file; 4, step 2 input file ; 5, step 3 input file ; 6, step 4 input file ; 7, step 5 input file ; 8. step 6 input file ; 9 step 7 input file
//argv 10: step 8 input file; argv 11: step 9 input file; argv 12: step 10 production input file; argv:13, GBSA general input file; 14, GBSA complex input file; 15, GBSA ligand input file; 16, GBSA receptor input file

int main(int argc, char* argv[])
{
    std::string receptor_pdb_path(argv[1]);
    std::string complex_pdb_path(argv[2]);

    MolecularModeling::Assembly receptor_pdb(receptor_pdb_path, gmml::InputFileType::PDB);
    MolecularModeling::Assembly complex_pdb(complex_pdb_path, gmml::InputFileType::PDB);

    AtomVector receptor_atoms = receptor_pdb.GetAllAtomsOfAssembly();
    int num_receptor_atoms = receptor_atoms.size();
    AtomVector complex_atoms = complex_pdb.GetAllAtomsOfAssembly();
    int num_complex_atoms = complex_atoms.size();

    ResidueVector complex_residues = complex_pdb.GetResidues();
    ResidueVector receptor_residues = receptor_pdb.GetResidues();
    int num_complex_residues = complex_residues.size(); 
    int num_receptor_residues = receptor_residues.size();

    std::vector<std::pair<int, int> > strand_breaks = DetectStrandBreaks(receptor_residues, receptor_atoms);
    std::vector<int> binding_site_water_oxygen_indices = DetectBindingSiteWaterOxygens(complex_atoms);
    std::pair<int, int> head_n_tail_c_indices = FindHeadNTailC(receptor_pdb);

    std::string step1_input_path(argv[3]);
    WriteStepOneMinimizationInputFile(step1_input_path, head_n_tail_c_indices, strand_breaks, binding_site_water_oxygen_indices);

    std::string step2_input_path(argv[4]);
    WriteStepTwoRelaxationInputFile(step2_input_path, head_n_tail_c_indices, strand_breaks, binding_site_water_oxygen_indices, num_receptor_residues);

    std::string step3_input_path(argv[5]);
    WriteStepThreeMinimizationInputFile(step3_input_path, head_n_tail_c_indices, strand_breaks, binding_site_water_oxygen_indices);

    std::string step4_input_path(argv[6]);
    WriteStepFourMinimizationInputFile(step4_input_path, head_n_tail_c_indices, strand_breaks, binding_site_water_oxygen_indices);

    std::string step5_input_path(argv[7]);
    WriteStepFiveMinimizationInputFile(step5_input_path, head_n_tail_c_indices, strand_breaks, binding_site_water_oxygen_indices);

    std::string step6_input_path(argv[8]);
    WriteStepSixRelaxationInputFile(step6_input_path, head_n_tail_c_indices, strand_breaks, binding_site_water_oxygen_indices, num_receptor_residues);

    std::string step7_input_path(argv[9]);
    WriteStepSevenRelaxationInputFile(step7_input_path, head_n_tail_c_indices, strand_breaks, binding_site_water_oxygen_indices, num_receptor_residues);

    std::string step8_input_path(argv[10]);
    WriteStepEightRelaxationInputFile(step8_input_path, head_n_tail_c_indices, strand_breaks, binding_site_water_oxygen_indices, num_receptor_residues);

    std::string step9_input_path(argv[11]);
    WriteStepNineRelaxationInputFile(step9_input_path, head_n_tail_c_indices, strand_breaks, binding_site_water_oxygen_indices, num_receptor_residues);

    std::string step10_input_path(argv[12]);
	WriteStepTenProductionInputFile(step10_input_path, num_complex_atoms, head_n_tail_c_indices, strand_breaks, binding_site_water_oxygen_indices, num_receptor_residues);

    std::string gbsa_general_input_file_path(argv[13]);
    WriteGBSAGeneralInputFile(gbsa_general_input_file_path);

	//std::pair<std::string, std::string> r_l_mask = GenerateReceptorLigandMask(complex_residues);
	//std::cout << "Receptor mask: " << r_l_mask.first << std::endl;
	//std::cout << "Ligand mask: " << r_l_mask.second << std::endl;

	std::pair<int, int> rec_begin_end;
	std::pair<int, int> lig_begin_end;
	GenerateReceptorLigandRange(complex_residues, rec_begin_end, lig_begin_end);

	std::cout << "Rec range: " << rec_begin_end.first << "-" << rec_begin_end.second << std::endl;
	std::cout << "Lig range: " << lig_begin_end.first << "-" << lig_begin_end.second << std::endl;

    std::string gbsa_complex_input_file_path(argv[14]);
    WriteGBSAComplexInputFile(gbsa_complex_input_file_path, num_complex_residues, num_receptor_residues, rec_begin_end, lig_begin_end);

    std::string gbsa_ligand_input_file_path(argv[15]);
    WriteGBSALigandInputFile(gbsa_ligand_input_file_path, num_complex_residues, num_receptor_residues, lig_begin_end);

    std::string gbsa_receptor_input_file_path(argv[16]);
    WriteGBSAReceptorInputFile(gbsa_receptor_input_file_path, num_receptor_residues, rec_begin_end);
	
	std::string gbsa_entropy_input_file_path(argv[17]);
	WriteGBSAQuasiHarmonicInputFile(gbsa_entropy_input_file_path, num_complex_residues, num_receptor_residues, rec_begin_end, lig_begin_end);

}
