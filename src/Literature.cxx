#include "Literature.h"

LitLifetime::LitLifetime(const LitLifetime& lt){

	StateIndex	= lt.StateIndex;
	Lifetime	= lt.Lifetime;
	DnUncertainty 	= lt.DnUncertainty;
	UpUncertainty	= lt.UpUncertainty;

}
LitLifetime& LitLifetime::operator = (const LitLifetime& lt){

	StateIndex	= lt.StateIndex;
	Lifetime	= lt.Lifetime;
	DnUncertainty 	= lt.DnUncertainty;
	UpUncertainty	= lt.UpUncertainty;

	return *this;

}

LitBranchingRatio::LitBranchingRatio(const LitBranchingRatio& lb){

	StateIndex_I1	= lb.StateIndex_I1;
	StateIndex_F1	= lb.StateIndex_F1;
	StateIndex_F2	= lb.StateIndex_F2;
	BranchingRatio	= lb.BranchingRatio;
	DnUncertainty	= lb.DnUncertainty;
	UpUncertainty	= lb.UpUncertainty;

}
LitBranchingRatio& LitBranchingRatio::operator = (const LitBranchingRatio& lb){

	StateIndex_I1	= lb.StateIndex_I1;
	StateIndex_F1	= lb.StateIndex_F1;
	StateIndex_F2	= lb.StateIndex_F2;
	BranchingRatio	= lb.BranchingRatio;
	DnUncertainty	= lb.DnUncertainty;
	UpUncertainty	= lb.UpUncertainty;

	return *this;	

}

LitMixingRatio::LitMixingRatio(const LitMixingRatio& lm){

	StateIndex_I	= lm.StateIndex_I;
	StateIndex_F	= lm.StateIndex_F;
	MixingRatio	= lm.MixingRatio;
	DnUncertainty	= lm.DnUncertainty;
	UpUncertainty	= lm.UpUncertainty;

}
LitMixingRatio& LitMixingRatio::operator = (const LitMixingRatio& lm){

	StateIndex_I	= lm.StateIndex_I;
	StateIndex_F	= lm.StateIndex_F;
	MixingRatio	= lm.MixingRatio;
	DnUncertainty	= lm.DnUncertainty;
	UpUncertainty	= lm.UpUncertainty;


	return *this;

}
