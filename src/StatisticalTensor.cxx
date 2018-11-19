#include "StatisticalTensor.h"

void StateTensor::AddElement(double rho, double k, double kappa){

	fTensor.push_back(rho);
	fK.push_back(k);
	fKappa.push_back(kappa);		

}

size_t StateTensor::IndexFromKkappa(double k, double kappa){

	for(size_t i = 0; i < fK.size(); i++)
		if(fK.at(i) == k && fKappa.at(i) == kappa)
			return i;
	
	return 0; 

}
