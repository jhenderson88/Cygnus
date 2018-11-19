#include "StoppingPower.h"

StoppingPower::StoppingPower(){
	fitted = false;
}

StoppingPower::StoppingPower(const StoppingPower& s){

	fitted 		= s.fitted;
	stoppingP	= s.stoppingP;
	stoppingE	= s.stoppingE;
	fit		= s.fit;

}
StoppingPower& StoppingPower::operator = (const StoppingPower& s){

	fitted 		= s.fitted;
	stoppingP	= s.stoppingP;
	stoppingE	= s.stoppingE;
	fit		= s.fit;

	return *this;

}

void StoppingPower::FitStoppingPowers(){

	TGraph gTmp;
	for(size_t e = 0; e < stoppingE.size(); e++)
		gTmp.SetPoint((int)e,stoppingE.at(e),stoppingP.at(e));

	TF1 *tmpFit = new TF1("stopping_fit","pol3",stoppingE.at(0),stoppingE.at(stoppingE.size()-1));
	gTmp.Fit(tmpFit,"RQ0");

	fit = *tmpFit;

	fitted = true;

}
