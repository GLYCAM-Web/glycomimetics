#include<string>
#include<fstream>
#include<cstdlib>
#include<pthread.h>
#include<iterator>
#include<cmath>

#include "../src/vina_atom_data.hpp"
#include "../src/utility.hpp" //glob
//#include "/cm/shared/apps/gems/gmml/includes/gmml.hpp"
//#include "/cm/shared/apps/gems/gmml/includes/MolecularModeling/assembly.hpp"
//#include "/cm/shared/apps/gems/gmml/includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp"
//#include "/cm/shared/apps/gems/gmml/includes/utils.hpp"

#include "includes/gmml.hpp"
#include "includes/MolecularModeling/assembly.hpp"
#include "includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp"
#include "includes/utils.hpp"


namespace PtheadTrajScoring{
    struct thread_args{
        std::vector<std::string>* file_path_vector_ptr = NULL;
        int start_index = 0;
        int end_index = 1;
        int thread_index = 0;
	std::map<int, std::vector<double> >* frame_score_map_ptr = NULL;
	int interval = 1;
	std::map<int, std::vector<int> >* bonding_map_ptr = NULL;
	std::vector<std::string>* atom_types_vector_ptr = NULL;
	std::string pdb_dir_path;

        thread_args(std::string& pdb_dir, std::vector<std::string>* vector_ptr, int index_start, int index_end, int index_thread, std::map<int, std::vector<double> >* score_map_ptr, int reading_interval, std::map<int, std::vector<int> >* bonding_map, std::vector<std::string>* atom_types){
	    pdb_dir_path = pdb_dir;
            file_path_vector_ptr = vector_ptr;
            start_index = index_start;
            end_index = index_end;
            thread_index = index_thread;
	    frame_score_map_ptr = score_map_ptr;
	    interval = reading_interval;
	    bonding_map_ptr = bonding_map;
	    atom_types_vector_ptr = atom_types;
        }
    };
}

void ScoreTraj(std::string pdb_dir_path, std::vector<std::string>* file_paths_vector_ptr, int start_index, int end_index, std::map<int, std::vector<double> >* frame_score_map_ptr, int interval, std::map<int, std::vector<int> >* bonding_map_ptr, std::vector<std::string>* atom_types_vector_ptr){
    for (unsigned int i = start_index; i <= end_index; i++){
	std::string pdb_full_path = pdb_dir_path; 
	pdb_full_path+="/";
	std::string file_path = file_paths_vector_ptr->at(i);
        pdb_full_path += file_path;

        unsigned int extension_index = file_path.find(".pdb");
        std::string frame_number_str = file_path.substr(0, extension_index);
        int frame_number = std::stoi(frame_number_str);
        std::map<int, std::vector<int> > bonding_map = *bonding_map_ptr;
        std::vector<std::string> atom_types = *atom_types_vector_ptr;
	
        if ( (frame_number -1) % interval == 0){
            MolecularModeling::Assembly cocomplex = MolecularModeling::Assembly(pdb_full_path, gmml::InputFileType::PDB);
            //VinaBondByDistance(cocomplex, vina_atom_data);
	    //Get bonding from map
	    AtomVector all_atoms = cocomplex.GetAllAtomsOfAssembly();
	    for (unsigned int j = 0; j < all_atoms.size(); j++){
                MolecularModeling::Atom* atom = all_atoms[j];
                MolecularModeling::AtomNode* new_node = new MolecularModeling::AtomNode();
                new_node->SetAtom(atom);
                atom->SetNode(new_node);

		atom->MolecularDynamicAtom::SetAtomType(atom_types[j]);
            }

	    for (std::map<int, std::vector<int> >::iterator mapit = bonding_map.begin(); mapit != bonding_map.end(); mapit++){
	        const int& atom_index = mapit->first;
	        std::vector<int>& neighbor_indices = mapit->second;

	        MolecularModeling::AtomNode* atom_node = all_atoms[atom_index]->GetNode();
	        for (unsigned int j = 0; j < neighbor_indices.size(); j++){
	            int& neighbor_index = neighbor_indices[j];
	            atom_node->AddNodeNeighbor(all_atoms[neighbor_index]);
	        }

	    }
	    
            std::vector<MolecularModeling::Residue*> residues = cocomplex.GetResidues();
	    std::vector<std::string> receptor_residue_names = {"CA", "WAT"};

            AtomVector receptor_atoms;
            AtomVector ligand_atoms;
	    std::vector<MolecularModeling::Residue*> receptor_residues;
	    std::vector<MolecularModeling::Residue*> ligand_residues;
            for (unsigned int j = 0; j < residues.size(); j++){
                AtomVector residue_atoms = residues[j]->GetAtoms();
	        std::string resname = residues[j]->GetName();

                if (residues[j]->CheckIfProtein() || std::find(receptor_residue_names.begin(), receptor_residue_names.end(), resname) != receptor_residue_names.end()){
                    receptor_atoms.insert(receptor_atoms.end(), residue_atoms.begin(), residue_atoms.end());
	            receptor_residues.push_back(residues[j]);
                }
                else{
                    ligand_atoms.insert(ligand_atoms.end(), residue_atoms.begin(), residue_atoms.end());
	            ligand_residues.push_back(residues[j]);
                }
            }
            //std::cout << "Num recetor atoms: " << receptor_atoms.size() << std::endl;

	    //Total, gau1, gau2, repulsion, phobic, hbond, catpi
	    VinaScorePrerequisites prerequisites(ligand_atoms, receptor_atoms);
	    std::vector<double> vina_scores = VinaScoreInPlace(prerequisites, 0);
            (*frame_score_map_ptr)[frame_number] = vina_scores;
            std::cout << frame_number << "\t" << vina_scores[0] << std::endl;
	}
    }
}

void* ScoringTrajByMultipleThreadsStartRoutine(void* arg_ptr){
    PtheadTrajScoring::thread_args* argstruct = (PtheadTrajScoring::thread_args*) arg_ptr;
    ScoreTraj(argstruct->pdb_dir_path, argstruct->file_path_vector_ptr, argstruct->start_index, argstruct->end_index, argstruct->frame_score_map_ptr, argstruct->interval, argstruct->bonding_map_ptr, argstruct->atom_types_vector_ptr);
    return NULL;
}

int main(int argc, char* argv[]){ //argv[1] is sample pdbqt path, argv[2] is directory path, argv[3] is interval

    std::string sample_pdbqt_path = std::string(argv[1]);
    std::string dir_path = std::string(argv[2]);
    int interval = std::stoi(std::string(argv[3]));

    std::string extension = "pdb";
    std::vector<std::string> file_paths = glob(dir_path, extension);

    //Obtain Autodock atom type and build bonding for the sample pdbqt file. 
    std::map<int, std::vector<int> > bonding_map;
    std::vector<std::string> atom_types;
    MolecularModeling::Assembly sample_pdbqt = MolecularModeling::Assembly(sample_pdbqt_path, gmml::InputFileType::PDBQT);
    VinaBondByDistance(sample_pdbqt, vina_atom_data);

    AtomVector atoms = sample_pdbqt.GetAllAtomsOfAssembly();
    AtomVector::iterator begin_pos = atoms.begin();

    for (unsigned int j = 0; j < atoms.size(); j++){
        MolecularModeling::Atom* atom = atoms[j];
	atom_types.push_back(atom->MolecularDynamicAtom::GetAtomType());
	AtomVector neighbors = atom->GetNode()->GetNodeNeighbors();

	for (unsigned int k = 0; k < neighbors.size(); k++){
	    MolecularModeling::Atom* neighbor = neighbors[k];
	    AtomVector::iterator neighbor_it = std::find(atoms.begin(), atoms.end(), neighbor);
	    int neighbor_index = std::distance(begin_pos, neighbor_it);
	    bonding_map[j].push_back(neighbor_index);
	}
    }

    std::cout << "In the beginning map size; " << bonding_map.size() << std::endl;
    std::cout << "Number of paths globbed: " << file_paths.size() << std::endl;

    std::string csv_filename = "output.csv";
    std::ofstream output_csv(csv_filename);
    if (output_csv.fail()){
        std::cout << "Failed to open file for writing" << std::endl;
	return 0;
    }

    std::cout << "Opened csv " << csv_filename << std::endl;
    std::vector<std::string> amino_libs;
    amino_libs.push_back("/home/yao/GLYCAM_Dev_Env/V_2/Web_Programs/gems/gmml/dat/CurrentParams/leaprc.ff12SB_2014-04-24/amino12.lib");
    amino_libs.push_back("/home/yao/GLYCAM_Dev_Env/V_2/Web_Programs/gems/gmml/dat/CurrentParams/leaprc.ff12SB_2014-04-24/aminoct12.lib");
    amino_libs.push_back("/home/yao/GLYCAM_Dev_Env/V_2/Web_Programs/gems/gmml/dat/CurrentParams/leaprc.ff12SB_2014-04-24/aminont12.lib");

    int num_threads = 14;
    std::vector<PtheadTrajScoring::thread_args> argstruct;
    std::map<int, std::vector<double> > frame_score_map;

    int start = 0, end = 0, structs_per_thread = file_paths.size()/num_threads;
    int last_thread_index = num_threads -1; 

    for (unsigned int i = 0; i < num_threads; i++){
	end = (i < last_thread_index) ? (start + structs_per_thread -1) : (file_paths.size()-1);
	PtheadTrajScoring::thread_args new_struct(dir_path, &file_paths, start, end, i, &frame_score_map, interval, &bonding_map, &atom_types);
	argstruct.push_back(new_struct);
	start = end + 1;
    }

    pthread_t tid[num_threads];
    //pthread_mutex_t lock;
    for (unsigned int i = 0; i < num_threads; i++){
        //Create a child thread
        int success_status = pthread_create(&tid[i], NULL, &ScoringTrajByMultipleThreadsStartRoutine, &argstruct[i]);
    }

    for (unsigned int i = 0; i < num_threads; i++){
        pthread_join(tid[i], NULL);
    }
    //pthread_mutex_destroy(&lock);
    
    std::map<int, std::vector<double> > frame_running_avg_map;
    std::map<int, std::vector<double> > frame_running_stdev_map;
    std::vector<double> running_sums(frame_score_map.begin()->second.size(), 0);
    std::vector<double> running_avgs(frame_score_map.begin()->second.size(), 0);
    std::vector<double> running_variance(frame_score_map.begin()->second.size(), 0);
    std::vector<double> running_stdevs(frame_score_map.begin()->second.size(), 0);

    for (std::map<int, std::vector<double> >::iterator mapit = frame_score_map.begin(); mapit != frame_score_map.end(); mapit++){
	int num_entries = std::distance(frame_score_map.begin(), mapit) + 1;
	for (unsigned int i = 0; i < mapit->second.size(); i++){
	    running_sums[i] += mapit->second[i];
	    if (i == 0){  //If it's the 1st term, which is total affinity, subtract reference affinity
	        running_avgs[i] = (double) (running_sums[i] / num_entries);
	    }
	    else{  //Otherwise, do not substarct reference affinity
	        running_avgs[i] = (double) (running_sums[i] / num_entries);
	    }
     
            if (i == 0){
	        running_variance[i] += std::pow(mapit->second[i] - running_avgs[i], 2);
	    }
	    else{
	        running_variance[i] += std::pow(mapit->second[i] - running_avgs[i], 2);
	    }

	    if (num_entries > 1){
	        running_stdevs[i] = std::pow(running_variance[i] / (num_entries -1), 0.5);
	    }
	}
	frame_running_avg_map[mapit->first] = running_avgs;
	frame_running_stdev_map[mapit->first] = running_stdevs;
    }

    output_csv << std::left << std::setw(7) << "Frame" << std::left << std::setw(7) << "RAvg" << std::left << std::setw(7) << "Stdev"; 
    output_csv << std::left << std::setw(7) << "Gau1" << std::left << std::setw(7) << "Stdev";
    output_csv << std::left << std::setw(7) << "Gau2" << std::left << std::setw(7) << "Stdev";
    output_csv << std::left << std::setw(7) << "Rep" << std::left << std::setw(7) << "Stdev";
    output_csv << std::left << std::setw(7) << "Pho" << std::left << std::setw(7) << "Stdev"; 
    output_csv << std::left << std::setw(7) << "Hb" << std::left << std::setw(7) << "Stdev";
    output_csv << std::left << std::setw(7) << "Catpi" << std::left << std::setw(7) << "Stdev";
    output_csv << std::left << std::setw(7) << "Pipi" << std::left << std::setw(7) << "Stdev";
    output_csv << std::left << std::setw(7) << "CHpi" << std::left << std::setw(7) << "Stdev";
    output_csv << std::left << std::setw(7) << "Tot" << std::left << std::setw(7) << std::endl;

    for (std::map<int, std::vector<double> >::iterator mapit = frame_running_avg_map.begin(); mapit != frame_running_avg_map.end(); mapit++){
	output_csv << std::left << std::setw(7) << mapit->first;
	for (unsigned int i = 0; i < mapit->second.size(); i++){
            output_csv << std::left << std::setw(7) << std::setprecision(2) << std::fixed << mapit->second[i] 
		       << std::left << std::setw(7) << std::setprecision(2) << std::fixed << frame_running_stdev_map[mapit->first][i];
	}
	output_csv << std::left << std::setw(7) << std::setprecision(2) << std::fixed << mapit->second[0] << "\n";
    }
    //output_csv << "\n";
    output_csv.close();
    return 0;
}

