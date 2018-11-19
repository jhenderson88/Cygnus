#include "State.h"

State::State(double E, double J){

	fStateEnergy 	= E;
	fJ		= J;

}


State::State(const State& s){

	fStateEnergy	= s.fStateEnergy;
	fJ		= s.fJ;
	fEta		= s.fEta;
	fPsi		= s.fPsi;

}

State& State::operator = (const State& s){

	fStateEnergy	= s.fStateEnergy;
	fJ		= s.fJ;
	fEta		= s.fEta;
	fPsi		= s.fPsi;

	return *this;

}
