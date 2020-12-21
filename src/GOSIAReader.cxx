#include "GOSIAReader.h"

GOSIAReader::GOSIAReader(Nucleus* nucl, const char *datafilename){
	gosiaData.clear();
	fNucleus = *nucl;
	ReadGOSIAFile(datafilename);
}

void	GOSIAReader::ReadGOSIAFile(const char* datafilename){

	std::ifstream infile(datafilename);

	int nExpt = 0;

	int initial_state;
	int final_state;
	double iJ;
	double fJ;
	double counts;
	double uncertainty;	// Not required for GOSIA

	ExperimentData tmpExpt;

	bool	flag = false;
	int	counter = 0;
	int	exptCounter = 0;

	std::string line;
	while(std::getline(infile,line)){
		if(flag){
			if(counter > 0){	// First GOSIA line after NORMALIZED YIELD is whitespace
				std::istringstream ss(line);
				ss	>> initial_state;
				ss	>> final_state;
				ss	>> iJ;
				ss	>> fJ;
				ss	>> counts;
				ss	>> uncertainty; 
				initial_state--;	// GOSIA to Cygnus numbering
				final_state--;		// GOSIA to Cygnus numbering
	
				if(initial_state < fNucleus.GetNstates() && final_state < fNucleus.GetNstates() && initial_state >= 0 && final_state >= 0){
					tmpExpt.AddData(initial_state,final_state,counts,uncertainty);
				}
				if(initial_state == 1 && final_state == 0){
					flag = false;
					nExpt++;
					gosiaData.push_back(tmpExpt);
					continue;
				}
			}
			counter++;
		} 
		std::size_t found = line.find("NORMALIZED YIELD"); // Comment
		if(found != std::string::npos){
			flag = true;
			counter = 0;
			tmpExpt.ClearData();
		}
	}	
}

void	GOSIAReader::Print() const{

	for(size_t i=0;i<gosiaData.size();i++)
		gosiaData.at(i).Print();

}
