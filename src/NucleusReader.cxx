#include "NucleusReader.h"

NucleusReader::NucleusReader() : fNucleus(NULL)
{
}

NucleusReader::NucleusReader(const char *filename) : fNucleus(NULL)
{
	ReadNucleusFile(filename);	
}

NucleusReader::NucleusReader(const NucleusReader &n) : fNucleus(n.fNucleus)
{
}

NucleusReader& NucleusReader::operator = (const NucleusReader &n)
{

	fNucleus	= n.fNucleus;

	return *this;

}

void NucleusReader::ReadNucleusFile(const char *filename){
	
	std::ifstream infile(filename);
	
	int counter = 0;
	int tmpA, tmpZ, tmpS, tmpL;
	int tmpI, tmpP; 
	double tmpE, tmpJ;

	int tmpI1, tmpI2, tmp_Lambda;
	std::string tmp_Lambda_s;
	double tmp_ME;

	std::string line;
	while(std::getline(infile,line)){
		std::size_t found = line.find("!");
		if(found == std::string::npos){
			std::istringstream ss(line);

			if(counter == 0){
				ss >> tmpA >> tmpZ >> tmpS >> tmpL;
				fNucleus = new Nucleus(tmpA,tmpZ,tmpS,tmpL);
			}
			else if(counter <= tmpS){
				ss >> tmpI >> tmpE >> tmpJ >> tmpP;
				fNucleus->SetState(tmpI,tmpE,tmpJ,tmpP);	
			}
			else{
				std::size_t foundE = line.find("E");
				std::size_t foundM = line.find("M");
				if(foundE != std::string::npos || foundM != std::string::npos){ // Found a transition specified E or M
					ss >> tmpI1 >> tmpI2 >> tmp_ME >> tmp_Lambda_s;
					if(foundE != std::string::npos && foundM !=std::string::npos){
						std::cout << "Both E and M found, not possible. Skipping line!" << std::endl;
						continue;
					}
					tmp_Lambda = tmp_Lambda_s[1] - '0'; // Since we know the second digit is the multipolarity
					if(foundM != std::string::npos)
						tmp_Lambda += 6;	
				}
				else{
					ss >> tmpI1 >> tmpI2 >> tmp_ME >> tmp_Lambda;	
				}
				fNucleus->SetMatrixElement(tmp_Lambda-1,tmpI1,tmpI2,tmp_ME);
			}

			counter++;
		}

	}

	//fNucleus->PrintNucleus();

}
