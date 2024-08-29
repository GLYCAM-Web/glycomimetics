#include <string>
#include "includes/CentralDataStructure/Readers/Pdb/pdbFile.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/Selections/residueSelections.hpp"

//Borrowed from 026.test.editPDB.cpp in GMML, originally written by Oliver. 
int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " inputFile.pdb\n";
        std::cout << "Example: " << argv[0] << " tests/inputs/4mbz.pdb\n";
        std::exit(EXIT_FAILURE);
    }
    // requirement: a chain ID for every single ATOM entry, and all ligand atoms should be put in a single residue.
    pdb::PdbFile pdbFile(argv[1]); // PdbFile is an "Ensemble" (made up of "Assemblies"), but if you want to just set
                                   // every molecule to have any chain ID you can do:
    for (auto& residue : pdbFile.getResidues())
    {
        static_cast<pdb::PdbResidue*>(residue)->setChainId("Y");
        // std::cout << "Set chain of " << residue->getStringId() << "\n";
    }
    // ResidueTypes are guessed upon input. Using that guess to find the ligand, can improve this if you need:
    std::vector<cds::Residue*> ligandResidues = cdsSelections::selectResiduesByType(pdbFile.getResidues(), cds::ResidueType::Undefined);
    if (ligandResidues.empty())
    {
        std::cout << "Error: no ligand residues found in input file\n";
        std::exit(EXIT_FAILURE);
    }

    cds::Residue* firstLigandResidue = ligandResidues.front(); // Everything has the same residue number as the first one.
    for (auto& ligandResidue : ligandResidues) // Each MODEL in PdbFile is converted into an "Assembly"
    {
        // std::cout << "Renumbering and renaming " << ligandResidue->getStringId() << "\n";
        ligandResidue->setNumber(firstLigandResidue->getNumber());
        ligandResidue->setName(firstLigandResidue->getName());
    }

	std::string output_pdb_path(argv[1]);
	int pdb_pos = output_pdb_path.find(".pdb");
	if (pdb_pos != std::string::npos){
		output_pdb_path.erase(pdb_pos, 4);
	}
	output_pdb_path += "_reformatted.pdb";
    pdbFile.Write(output_pdb_path);
    return 0;
}
