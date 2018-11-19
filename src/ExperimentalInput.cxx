#include "ExperimentalInput.h"

ExptData::ExptData(const ExptData& e){

 	StateIndex_I	= e.StateIndex_I;
  	StateIndex_F	= e.StateIndex_F;
	Counts		= e.Counts;
	UpUncertainty 	= e.UpUncertainty;
	DnUncertainty	= e.DnUncertainty;
	

}
ExptData& ExptData::operator = (const ExptData& e){

 	StateIndex_I	= e.StateIndex_I;
  	StateIndex_F	= e.StateIndex_F;
	Counts		= e.Counts;
	UpUncertainty 	= e.UpUncertainty;
	DnUncertainty	= e.DnUncertainty;

	return *this;

}
ExperimentData::ExperimentData(const ExperimentData& e){

	Data.resize(e.Data.size());
	for(size_t d=0;d<e.Data.size();d++)
		Data.at(d) = e.Data.at(d);		
	thetaCM		= e.thetaCM;
	

}
ExperimentData& ExperimentData::operator = (const ExperimentData& e){

	Data.resize(e.Data.size());
	for(size_t d=0;d<e.Data.size();d++)
		Data.at(d) = e.Data.at(d);		
	thetaCM		= e.thetaCM;

	return *this;

}

void ExperimentData::Print() const {

	if(thetaCM < 0)
		std::cout << "Theta undefined" << std::endl;
	else
		std::cout << "Theta [CM]: " << thetaCM << " degrees" << std::endl;

	std::cout 	<< std::setw(9) << std::left << "Initial"
			<< std::setw(7) << std::left << "Final"
			<< std::setw(10) << std::left << "Counts"
			<< std::setw(10) << std::left << "UpUncert."
			<< std::setw(10) << std::left << "DnUncert."
			<< std::endl;
	for(unsigned int i=0;i<Data.size();i++)
		Data.at(i).Print();

}

void ExptData::Print() const {
	std::cout 	<< std::setw(9) << std::left << StateIndex_I 
			<< std::setw(7) << std::left << StateIndex_F 
			<< std::setw(10) << std::left << Counts 
			<< std::setw(10) << UpUncertainty 
			<< std::setw(10) << DnUncertainty 
			<< std::endl; 
}
