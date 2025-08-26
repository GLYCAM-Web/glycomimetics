#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <sstream>

std::vector<std::string> tokenize_by_delim(std::string str, std::string delim){
	std::vector<std::string> tokens;
	std::string token;
	size_t pos = 0;

	while ((pos = str.find(delim)) != std::string::npos) {
		token = str.substr(0, pos);
    	tokens.emplace_back(token);
    	str.erase(0, pos + delim.length());
	}

	tokens.emplace_back(str);
	return tokens;
}

void read_section(std::ifstream& per_frame_breakdown, std::string& line, std::string& delim, std::vector<double>& totals_per_frame){
	while(std::getline(per_frame_breakdown, line)){
		//std::cout << "This line: " << line << std::endl;
        if (line.find(delim) == std::string::npos) break; //Two blank lines at end of file, skip.

       	std::vector<std::string> tokens = tokenize_by_delim(line, delim);
       	double total = std::stod(tokens.back());
       	totals_per_frame.emplace_back(total);

       	/*std::cout << "This line: " << line << std::endl;
       	std::cout << "Tokens back " << tokens.back() << std::endl;
       	for (unsigned int i = 0; i < tokens.size(); i++){
           	std::cout << "Tokens: " << tokens[i] << std::endl;
       	}*/
	}
	return;
}

//double ravg, rstdev, r_exp_avg, r_ie, r_c2;
void process_section(std::vector<double>& totals_per_frame, double& avg, double& stdev, double& exp_avg, double& ie, double& c2){
	double num_frames = (double) totals_per_frame.size();
    avg = 0;
    for (unsigned int i = 0; i < totals_per_frame.size(); i++){
        avg += totals_per_frame[i];
    }
    avg /= num_frames;

    stdev = 0;
    for (unsigned int i = 0; i < totals_per_frame.size(); i++){
        stdev += (totals_per_frame[i] - avg) * (totals_per_frame[i] - avg);
    }
    stdev = std::sqrt(stdev / num_frames);

    double R = 0.001987, T = 298.15; //R=0.001987kcal/mol, T=298.15K
    exp_avg = 0;
    for (unsigned int i = 0; i < totals_per_frame.size(); i++){
        exp_avg += std::exp((totals_per_frame[i] - avg) / (R*T) );
    }
    exp_avg /= num_frames;

    ie = R * T * std::log(exp_avg);
    c2 = stdev * stdev / (2 * R * T);
}

int main(int argc, char* argv[]){
	std::string compound_name(argv[1]);
	std::string gbsa_file_path(argv[2]);
	std::ifstream per_frame_breakdown(gbsa_file_path);

	if (per_frame_breakdown.fail()){
		std::cout << "Interaction entropy: Could not open " << gbsa_file_path << " for reading." << std::endl;
		std::exit(1);
	}

	std::string line;
	std::string delim = ",";
	bool contains_receptor = false, contains_ligand = false, contains_delta = false;
	std::vector<double> ip_energies;
	std::vector<double> il_energies;
	std::vector<double> rl_energies;

	while(std::getline(per_frame_breakdown, line)){
		if (line.find("Receptor Energy Terms") != std::string::npos){
			contains_receptor =  true;
			std::getline(per_frame_breakdown, line);
			read_section(per_frame_breakdown, line, delim, ip_energies);		
		}
		else if (line.find("Ligand Energy Terms") != std::string::npos){
			contains_ligand = true;
			std::getline(per_frame_breakdown, line);
			read_section(per_frame_breakdown, line, delim, il_energies);		
		}
		if (line.find("DELTA Energy Terms") != std::string::npos){
			contains_delta = true;
			std::getline(per_frame_breakdown, line);
			read_section(per_frame_breakdown, line, delim, rl_energies);		
		}
	}
	per_frame_breakdown.close();

	if (!contains_receptor){
		std::cout << "Could not find Receptor Energy Terms section from the input file. Check your input file" << std::endl;
		std::exit(1);
	}
	if (!contains_ligand){
		std::cout << "Could not find Ligand Energy Terms section from the input file. Check your input file" << std::endl;
		std::exit(1);
	}
	if (!contains_delta){
		std::cout << "Could not find DELTA Energy Terms section from the input file. Check your input file" << std::endl;
		std::exit(1);
	}

	double ravg, rstdev, r_exp_avg, r_ie, r_c2;
	double lavg, lstdev, l_exp_avg, l_ie, l_c2;
	double rlavg, rlstdev, rl_exp_avg, rl_ie, rl_c2;

	process_section(ip_energies, ravg, rstdev, r_exp_avg, r_ie, r_c2);
	process_section(il_energies, lavg, lstdev, l_exp_avg, l_ie, l_c2);
	process_section(rl_energies, rlavg, rlstdev, rl_exp_avg, rl_ie, rl_c2);

	double stdev_cutoff = 3.58269; //15 kj/mol = 3.58269 kcal/mol. Literature says if stdev > 15 KJ/mol it's impossible to converge the entropy. 
	std::string stable = (rlstdev >= stdev_cutoff) ? "false" : "true";
	
	std::cout << compound_name << " " << rlavg << " " << rlstdev << " " << r_ie << " " << l_ie << " " << rl_ie
							   << " " << rlavg + rl_ie
							   << " " << r_c2 << " " << l_c2 << " " << rl_c2
							   << " " << rlavg + rl_c2
							   << std::endl;
	return 0;
}
