#include "Nucleus.h"

// Empty constructor
Nucleus::Nucleus()
{
}

Nucleus::Nucleus(int Z, int A, int nS, int nL){
	SetZ(Z);
	SetA(A);
	SetMaxLambda(nL);
	SetNstates(nS);
}
Nucleus::Nucleus(const Nucleus& n) : MatrixElements(n.MatrixElements.size()), MatrixElementsUL(n.MatrixElementsUL.size()), MatrixElementsLL(n.MatrixElementsLL.size()) {

	nStates 	= n.nStates;			// Number of states in the nucleus
	maxLambda	= n.maxLambda;			// Number of multipolarities to be considered

	nucleusZ	= n.nucleusZ;			// Proton number of the nucleus
	nucleusA 	= n.nucleusA;			// Nucleon number of the nucleus

	LevelEnergies 	= n.LevelEnergies;
	LevelJ  	= n.LevelJ;
	LevelP		= n.LevelP;


	for(std::size_t i=0;i<n.MatrixElements.size();i++){
		MatrixElements.at(i).ResizeTo(n.MatrixElements.at(i).GetNcols(),n.MatrixElements.at(i).GetNrows());
		MatrixElements.at(i) = n.MatrixElements.at(i);
	}

	for(std::size_t i=0;i<n.MatrixElementsUL.size();i++){
		MatrixElementsUL.at(i).ResizeTo(n.MatrixElementsUL.at(i).GetNcols(),n.MatrixElementsUL.at(i).GetNrows());
		MatrixElementsUL.at(i) = n.MatrixElementsUL.at(i);
	}

	for(std::size_t i=0;i<n.MatrixElementsLL.size();i++){
		MatrixElementsLL.at(i).ResizeTo(n.MatrixElementsLL.at(i).GetNcols(),n.MatrixElementsLL.at(i).GetNrows());
		MatrixElementsLL.at(i) = n.MatrixElementsLL.at(i);
	}

}
Nucleus& Nucleus::operator = (const Nucleus& n){

	nStates 	= n.nStates;			// Number of states in the nucleus
	maxLambda	= n.maxLambda;			// Number of multipolarities to be considered

	nucleusZ	= n.nucleusZ;			// Proton number of the nucleus
	nucleusA 	= n.nucleusA;			// Nucleon number of the nucleus

	LevelEnergies 	= n.LevelEnergies;
	LevelJ  	= n.LevelJ;
	LevelP		= n.LevelP;


	MatrixElements.resize((size_t)n.MatrixElements.size());
	for(size_t i=0;i<n.MatrixElements.size();i++){
		MatrixElements.at(i).ResizeTo(n.MatrixElements.at(i).GetNrows(),n.MatrixElements.at(i).GetNcols());
		MatrixElements.at(i) = n.MatrixElements.at(i);
	}

	MatrixElementsUL.resize((size_t)n.MatrixElementsUL.size());
	for(size_t i=0;i<n.MatrixElementsUL.size();i++){
		MatrixElementsUL.at(i).ResizeTo(n.MatrixElementsUL.at(i).GetNrows(),n.MatrixElementsUL.at(i).GetNcols());
		MatrixElementsUL.at(i) = n.MatrixElementsUL.at(i);
	}

	MatrixElementsLL.resize((size_t)n.MatrixElementsLL.size());
	for(size_t i=0;i<n.MatrixElementsLL.size();i++){
		MatrixElementsLL.at(i).ResizeTo(n.MatrixElementsLL.at(i).GetNrows(),n.MatrixElementsLL.at(i).GetNcols());
		MatrixElementsLL.at(i) = n.MatrixElementsLL.at(i);
	}

	return *this;

}


void Nucleus::SetMaxLambda(int mL){
	maxLambda = mL;
	MatrixElements.resize(mL);
	MatrixElementsUL.resize(mL);
	MatrixElementsLL.resize(mL);	
}

void Nucleus::SetNstates(int nS){
	nStates = nS;

	LevelEnergies.resize(nS);
	LevelJ.resize(nS);
	LevelP.resize(nS);

	for(int i=0;i<maxLambda;i++){
		MatrixElements[i].ResizeTo(nS,nS);
		MatrixElementsUL[i].ResizeTo(nS,nS);
		MatrixElementsLL[i].ResizeTo(nS,nS);
	}
	
}

void Nucleus::SetState(int s, double E, double J, int P){
	SetStateE(s,E);
	SetStateJ(s,J);
	SetStateP(s,P);
}

void Nucleus::SetStateE(int s, double E){
	if(s < nStates && s >= 0)
		LevelEnergies[s] = E;	
	else
		std::cout << "Energy: State " << s << " out of range, " << nStates << " declared" << std::endl;
}
void Nucleus::SetStateP(int s, int P){
	if(s < nStates && s >= 0)
		LevelP[s] = P;	
	else
		std::cout << "Parity: State " << s << " out of range, " << nStates << " declared" << std::endl;
}
void Nucleus::SetStateJ(int s, double J){
	if(s < nStates && s >= 0)
		LevelJ[s] = J;	
	else
		std::cout << "Spin: State " << s << " out of range, " << nStates << " declared" << std::endl;
}

void Nucleus::SetMatrixElement(int l, int s1, int s2, double me){

	std::string mult[7] = {"E1","E2","E3","E4","E5","E6","M1"};
	int	lambda = l+1; // This lambda is only used for checking spin/parity conservation
	if(l >= 6) // Magnetic - opposite parity conventions to electric
		lambda -= 5;
	int	dP_Lambda = TMath::Power(-1,lambda);
	int	dP = GetLevelP().at(s1)/GetLevelP().at(s2);
	if(dP_Lambda != dP){
		std::cout << "Parity conservation broken! Check state parities and transition multipolarities." << std::endl;
		std::cout 	<< std::setw(25) << std::left << "Transition mult:" 
				<< std::setw(25) << std::left << "P change (transition):" 
				<< std::setw(25) << std::left << "P change (levels):"
				<< std::endl;
		std::cout 	<< std::setw(25) << std::left << mult[l] 
				<< std::setw(25) << std::left << dP_Lambda 
				<< std::setw(25) << std::left << dP 
				<< std::endl;
		return;
	}
	if(l >= 6)
		lambda--;
	int	dJ = TMath::Abs(GetLevelJ().at(s1) - GetLevelJ().at(s2));
	if(dJ > lambda){
		std::cout << "Angular momentum conservation broken! Check state spins and transition multipolarities." << std::endl;
		std::cout 	<< std::setw(25) << std::left << "Transition mult:" 
				<< std::setw(30) << std::left << "Max J change (transition):" 
				<< std::setw(25) << std::left << "J change (levels):"
				<< std::endl;
		std::cout 	<< std::setw(20) << std::left << mult[l] 
				<< std::setw(30) << std::left << lambda 
				<< std::setw(20) << std::left << dJ 
				<< std::endl;
		return;
	}

	if(l < maxLambda && l >= 0){
		if(s1 < nStates && s1 >= 0 && s2 < nStates && s2 >=0){
			MatrixElements.at(l)[s1][s2]=me;
			MatrixElements.at(l)[s2][s1]=me;
		}
		else{
			if(s1 >= nStates || s1 < 0)
				std::cout << "Matrix element: State " <<  s1 << " out of range, " << nStates << " declared" << std::endl;
			if(s2 >= nStates || s2 < 0)
				std::cout << "Matrix element: State " <<  s2 << " out of range, " << nStates << " declared" << std::endl;
		}
	}
	else{
		std::cout << "Matrix element: Multipolarity " << l << " out of range, " << maxLambda << " declared" << std::endl;
		MultipolarityReminder();
	}

}

void Nucleus::MultipolarityReminder() const{
	std::cout << "E1 = 1\nE2 = 2\nE3 = 3\nE4 = 4\nE5 = 5\nE6 = 5\nM1 = 7\nM2 = 8" << std::endl;
}

void Nucleus::PrintNucleus() const{

	std::string mult[7] = {"E1","E2","E3","E4","E5","E6","M1"};

	std::cout << "Z: " << GetZ() << ", A: " << GetA() << std::endl;
	std::cout << "Number of states: " << GetNstates() << std::endl;
	std::cout << "Max multipolarity: " << GetMaxLambda() << std::endl;

	std::cout << "States [E, J, Pi]:" << std::endl;
	for(unsigned int s = 0; s < GetLevelEnergies().size();s++)
		std::cout << GetLevelEnergies().at(s) << "\t" << GetLevelJ().at(s) << "\t" << GetLevelP().at(s) << std::endl;


	std::cout << "Matrix elements: " << std::endl;
	for(unsigned int l = 0; l < GetMatrixElements().size(); l++){
		double max = 0;
		for(int i=0;i<GetMatrixElements().at(l).GetNcols();i++){
			for(int j=0;j<GetMatrixElements().at(l).GetNrows();j++){
				if(GetMatrixElements().at(l)[j][i] > max) max = GetMatrixElements().at(l)[j][i];
			}
		}

		if(max > 0){
			std::cout << "Lambda: " << mult[l] << std::endl;
			GetMatrixElements().at(l).Print();
		}
		
	}

}

void Nucleus::PrintState(int s) const{
	std::cout << GetLevelEnergies().at(s) << "\t" << GetLevelJ().at(s) << "\t" << GetLevelP().at(s) << std::endl;
}	
