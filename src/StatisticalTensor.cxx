#include "StatisticalTensor.h"

//*************************** State tensor ************************************//

//	Copy operator
StateTensor::StateTensor(const StateTensor& s){

	fStateIndex	= s.fStateIndex;
	fTensor		= s.fTensor;
	fK		= s.fK;	
	fKappa		= s.fKappa;
	fMaxK		= s.fMaxK;	

}
//	Assignment operator:
StateTensor& StateTensor::operator = (const StateTensor &s){

	fStateIndex	= s.fStateIndex;
	fTensor		= s.fTensor;
	fK		= s.fK;	
	fKappa		= s.fKappa;
	fMaxK		= s.fMaxK;

	return		*this;

}

void StateTensor::AddElement(double rho, double k, double kappa){

	fTensor.push_back(rho);
	fK.push_back(k);
	fKappa.push_back(kappa);		

	if(k > fMaxK)
		fMaxK = k;

}

size_t StateTensor::IndexFromKkappa(double k, double kappa){

	for(size_t i = 0; i < fK.size(); i++)
		if(fK.at(i) == k && fKappa.at(i) == kappa)
			return i;
	
	return 0; 

}

//************************** Statistical Tensor *********************************// 

//	Copy operator
StatisticalTensor::StatisticalTensor(const StatisticalTensor& s){

	fStateTensors.clear();
	fStateTensors.resize(s.fStateTensors.size());
	for(size_t i=0;i<s.fStateTensors.size();i++)
		fStateTensors[i] = s.fStateTensors[i];
	
}
//	Assignment operator:
StatisticalTensor& StatisticalTensor::operator = (const StatisticalTensor &s){

	fStateTensors.clear();
	fStateTensors.resize(s.fStateTensors.size());
	for(size_t i=0;i<s.fStateTensors.size();i++)
		fStateTensors[i] = s.fStateTensors[i];

	return		*this;

}
