#include "ScalingFitFCN.h"


double ScalingFitFCN::operator()(const double* par){

	double chisq = 0;

	double	scaling	= par[0];

	for(size_t i=0;i<exptData.size();i++){
		chisq	+= pow((scaling * calcData.at(i) - exptData.at(i))/exptUnc.at(i),2);
	}

	return	chisq;

}
