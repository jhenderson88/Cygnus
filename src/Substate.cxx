#include "Substate.h"

Substate::Substate(double M, int ss, int s){

	fM		= M;
	fSubstateIndex	= ss;
	fStateIndex	= s;
	fMirrorIndex	= ss;

}


Substate::Substate(const Substate& ss){

	fM		= ss.fM;
	fSubstateIndex	= ss.fSubstateIndex;
	fStateIndex	= ss.fStateIndex;
	fMirrorIndex	= ss.fMirrorIndex;
	fConnections.clear();
	fConnections.resize(ss.fConnections.size());
	for(size_t c = 0; c < fConnections.size(); c++)
		fConnections.at(c) = ss.fConnections.at(c);

}

Substate& Substate::operator = (const Substate& ss){

	fM		= ss.fM;
	fSubstateIndex	= ss.fSubstateIndex;
	fStateIndex	= ss.fStateIndex;
	fMirrorIndex	= ss.fMirrorIndex;
	fConnections.clear();
	fConnections.resize(ss.fConnections.size());
	for(size_t c = 0; c < fConnections.size(); c++)
		fConnections.at(c) = ss.fConnections.at(c);

	return *this;

}
