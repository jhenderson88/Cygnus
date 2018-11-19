#include "DataReader.h"

DataReader::DataReader(Nucleus* nucl, const char *datafilename){
	exptData.clear();
	fNucleus = *nucl;
	ReadDataFile(datafilename);
}

void	DataReader::ReadDataFile(const char* datafilename){

	std::ifstream infile(datafilename);

	int nExpt = 0;

	int initial_state;
	int final_state;
	double counts;
	double uncertainty;

	ExperimentData tmpExpt;

	std::string line;
	while(std::getline(infile,line)){
		std::size_t found = line.find("!"); // Comment
		if(found == std::string::npos){
			std::istringstream ss(line);
			found = line.find("EXPT");
			if(found != std::string::npos && nExpt !=0){
				exptData.push_back(tmpExpt);	
				tmpExpt.ClearData();
				nExpt++;
				continue;
			}
			else if(found != std::string::npos){
				nExpt++;
				continue;
			}
			found = line.find("END");
			if(found != std::string::npos){
				exptData.push_back(tmpExpt);
				tmpExpt.ClearData();
				break;
			}
			ss >> initial_state >> final_state >> counts >> uncertainty;
			if(initial_state < fNucleus.GetNstates() && final_state < fNucleus.GetNstates() && initial_state >= 0 && final_state >= 0){
				tmpExpt.AddData(initial_state,final_state,counts,uncertainty);
			}
			else{
				if(initial_state < fNucleus.GetNstates())
					std::cout << "State " << initial_state << " out of range" << std::endl;
				if(final_state < fNucleus.GetNstates())
					std::cout << "State " << final_state << " out of range" << std::endl;
				fNucleus.PrintNucleus();
			}
		}
	}	

}
