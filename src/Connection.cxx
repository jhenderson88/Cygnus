#include "Connection.h"

Connection::Connection(const Connection& c){
	
	fXi.resize(c.fXi.size());
	for(size_t i = 0; i < fXi.size(); i++)	
		fXi[i] = c.fXi[i];
	fPsi.resize(c.fPsi.size());
	for(size_t i = 0; i < fPsi.size(); i++)	
		fPsi[i]	= c.fPsi[i];
	fZeta.resize(c.fZeta.size());
	for(size_t i = 0; i < fZeta.size(); i++)	
		fZeta[i] = c.fZeta[i];
	fLambda.resize(c.fLambda.size());
	for(size_t i = 0; i < fLambda.size(); i++)
		fLambda[i] = c.fLambda[i];
	fConnectedState		= c.fConnectedState;

}
Connection& Connection::operator = (const Connection& c){

	fXi.resize(c.fXi.size());
	for(size_t i = 0; i < fXi.size(); i++)	
		fXi[i] = c.fXi[i];
	fPsi.resize(c.fPsi.size());
	for(size_t i = 0; i < fPsi.size(); i++)	
		fPsi[i]	= c.fPsi[i];
	fZeta.resize(c.fZeta.size());
	for(size_t i = 0; i < fZeta.size(); i++)	
		fZeta[i] = c.fZeta[i];
	fLambda.resize(c.fLambda.size());
	for(size_t i = 0; i < fLambda.size(); i++)
		fLambda[i] = c.fLambda[i];
	fConnectedState		= c.fConnectedState;

	return *this;

}

void Connection::AddLambda(int lambda, double xi, double psi, double zeta){
	
	fXi[lambda]		= xi;
	fPsi[lambda]		= psi;
	fZeta[lambda]	 	= zeta;
	fLambda[lambda]		= true;

}

void Connection::Clear(){

	fXi.clear();
	fPsi.clear();
	fZeta.clear();
	fLambda.clear();

}
