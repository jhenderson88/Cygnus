#include "TransitionRates.h"

TransitionRates::TransitionRates() : fNucleus(NULL)
{
}

TransitionRates::TransitionRates(Nucleus *nucl) : fNucleus(nucl)
{

	StateJ = nucl->GetLevelJ();
	StateE = nucl->GetLevelEnergies();

	SetMatrixElements();
	
}

TransitionRates::TransitionRates(const TransitionRates& t){

	fNucleus 			= t.fNucleus;

	MatrixElements.resize(t.MatrixElements.size());
	for(size_t i=0;i<MatrixElements.size();i++){
		MatrixElements[i].ResizeTo(t.MatrixElements[i].GetNcols(),t.MatrixElements[i].GetNrows());
		MatrixElements[i] 		= t.MatrixElements[i];
	}
	TransitionStrengths.resize(t.TransitionStrengths.size());
	for(size_t i=0;i<TransitionStrengths.size();i++){
		TransitionStrengths[i].ResizeTo(t.TransitionStrengths[i].GetNcols(),t.TransitionStrengths[i].GetNrows());
		TransitionStrengths[i] 		= t.TransitionStrengths[i];
	}
	TransitionStrengths_Abs.resize(t.TransitionStrengths_Abs.size());
	for(size_t i=0;i<TransitionStrengths_Abs.size();i++){
		TransitionStrengths_Abs[i].ResizeTo(t.TransitionStrengths_Abs[i].GetNcols(),t.TransitionStrengths_Abs[i].GetNrows());
		TransitionStrengths_Abs[i] 	= t.TransitionStrengths_Abs[i];
	}
	TransitionAmplitudes.resize(t.TransitionAmplitudes.size());
	for(size_t i=0;i<TransitionAmplitudes.size();i++){
		TransitionAmplitudes[i].ResizeTo(t.TransitionAmplitudes[i].GetNcols(),t.TransitionAmplitudes[i].GetNrows());
		TransitionAmplitudes[i] 	= t.TransitionAmplitudes[i];
	}
	SummedTransitionStrengths.ResizeTo(t.SummedTransitionStrengths.GetNcols(),t.SummedTransitionStrengths.GetNrows());
	SummedTransitionStrengths 	= t.SummedTransitionStrengths;
	Lifetimes.ResizeTo(t.Lifetimes.GetNcols(),t.Lifetimes.GetNrows());
	Lifetimes 			= t.Lifetimes;
	BranchingRatios.ResizeTo(t.BranchingRatios.GetNcols(),t.BranchingRatios.GetNrows());
	BranchingRatios 		= t.BranchingRatios;
	MixingRatios.ResizeTo(t.MixingRatios.GetNcols(),t.MixingRatios.GetNrows());
	MixingRatios 			= t.MixingRatios;

	StateJ 				= t.StateJ; 
	StateE 				= t.StateE;

	StateLifetimes.ResizeTo(t.StateLifetimes.GetNrows());
	StateLifetimes 			= t.StateLifetimes;
	StateDecayProb.ResizeTo(t.StateDecayProb.GetNrows());
	StateDecayProb 			= t.StateDecayProb;

	nDecays				= t.nDecays;

}

TransitionRates& TransitionRates::operator = (const TransitionRates& t){

	fNucleus 			= t.fNucleus;

	MatrixElements.resize(t.MatrixElements.size());
	for(size_t i=0;i<MatrixElements.size();i++){
		MatrixElements[i].ResizeTo(t.MatrixElements[i].GetNcols(),t.MatrixElements[i].GetNrows());
		MatrixElements[i] 		= t.MatrixElements[i];
	}
	TransitionStrengths.resize(t.TransitionStrengths.size());
	for(size_t i=0;i<TransitionStrengths.size();i++){
		TransitionStrengths[i].ResizeTo(t.TransitionStrengths[i].GetNcols(),t.TransitionStrengths[i].GetNrows());
		TransitionStrengths[i] 		= t.TransitionStrengths[i];
	}
	TransitionStrengths_Abs.resize(t.TransitionStrengths_Abs.size());
	for(size_t i=0;i<TransitionStrengths_Abs.size();i++){
		TransitionStrengths_Abs[i].ResizeTo(t.TransitionStrengths_Abs[i].GetNcols(),t.TransitionStrengths_Abs[i].GetNrows());
		TransitionStrengths_Abs[i] 	= t.TransitionStrengths_Abs[i];
	}
	TransitionAmplitudes.resize(t.TransitionAmplitudes.size());
	for(size_t i=0;i<TransitionAmplitudes.size();i++){
		TransitionAmplitudes[i].ResizeTo(t.TransitionAmplitudes[i].GetNcols(),t.TransitionAmplitudes[i].GetNrows());
		TransitionAmplitudes[i] 	= t.TransitionAmplitudes[i];
	}
	SummedTransitionStrengths.ResizeTo(t.SummedTransitionStrengths.GetNcols(),t.SummedTransitionStrengths.GetNrows());
	SummedTransitionStrengths 	= t.SummedTransitionStrengths;
	Lifetimes.ResizeTo(t.Lifetimes.GetNcols(),t.Lifetimes.GetNrows());
	Lifetimes 			= t.Lifetimes;
	BranchingRatios.ResizeTo(t.BranchingRatios.GetNcols(),t.BranchingRatios.GetNrows());
	BranchingRatios 		= t.BranchingRatios;
	MixingRatios.ResizeTo(t.MixingRatios.GetNcols(),t.MixingRatios.GetNrows());
	MixingRatios 			= t.MixingRatios;

	StateJ 				= t.StateJ; 
	StateE 				= t.StateE;

	StateLifetimes.ResizeTo(t.StateLifetimes.GetNrows());
	StateLifetimes 			= t.StateLifetimes;
	StateDecayProb.ResizeTo(t.StateDecayProb.GetNrows());
	StateDecayProb 			= t.StateDecayProb;

	nDecays				= t.nDecays;

	return *this;

}

// SYNTAX NOTE:
// COLUMN 	== INITIAL STATE
// ROW		== FINAL STATE
// TMATRIX[ROW][COLUMN]

void TransitionRates::SetMatrixElements(){
	
	double nbarns[7] = {1,2,3,4,5,6,0};
	double multfactor[7] = {0.629e-15,816e-12,1760e-6,5882,2.89e10,1.95e17,56.8e-15};	// Factors for B --> Lifetime
	int power[7] = {3,5,7,9,11,13,3};

	MatrixElements.resize(fNucleus->GetMatrixElements().size());
	for(unsigned int i=0;i<MatrixElements.size();i++){
		MatrixElements.at(i).ResizeTo(fNucleus->GetMatrixElements().at(i).GetNcols(),fNucleus->GetMatrixElements().at(i).GetNrows());
		MatrixElements.at(i) = fNucleus->GetMatrixElements().at(i);
	}

	TransitionStrengths.resize(fNucleus->GetMatrixElements().size());
	TransitionStrengths_Abs.resize(fNucleus->GetMatrixElements().size());
	TransitionAmplitudes.resize(fNucleus->GetMatrixElements().size());
	for(unsigned int i=0;i<MatrixElements.size();i++){
		TransitionStrengths.at(i).ResizeTo(fNucleus->GetMatrixElements().at(i).GetNcols(),fNucleus->GetMatrixElements().at(i).GetNrows());
		TransitionStrengths_Abs.at(i).ResizeTo(fNucleus->GetMatrixElements().at(i).GetNcols(),fNucleus->GetMatrixElements().at(i).GetNrows());
		TransitionAmplitudes.at(i).ResizeTo(fNucleus->GetMatrixElements().at(i).GetNcols(),fNucleus->GetMatrixElements().at(i).GetNrows());
		for(int x=0;x<TransitionStrengths.at(i).GetNcols();x++){
			for(int y=0;y<TransitionStrengths.at(i).GetNrows();y++){
				if(x>y){
					TransitionStrengths.at(i)[y][x] = TMath::Power(MatrixElements.at(i)[y][x],2) / (2 * StateJ.at(x) + 1) * TMath::Power(100,nbarns[i]); 
					TransitionStrengths_Abs.at(i)[y][x] = (TMath::Power(TMath::Abs(StateE.at(x)-StateE.at(y)),power[i]) * TransitionStrengths.at(i)[y][x]) / multfactor[i];
					if(i==0)
						TransitionAmplitudes.at(i)[y][x] = 398.77393 * TMath::Power(TMath::Abs(StateE.at(x)-StateE.at(y)),3./2.) / (2 * StateJ.at(x) + 1);
					else if(i==2)
						TransitionAmplitudes.at(i)[y][x] = 3.5002636 * TMath::Power(TMath::Abs(StateE.at(x)-StateE.at(y)),5./2.) / (2 * StateJ.at(x) + 1);
					else if(i==3)
						TransitionAmplitudes.at(i)[y][x] = 0.0238913 * TMath::Power(TMath::Abs(StateE.at(x)-StateE.at(y)),7./2.) / (2 * StateJ.at(x) + 1);
					else if(i==6)
						TransitionAmplitudes.at(i)[y][x] = 4.1923861 * TMath::Power(TMath::Abs(StateE.at(x)-StateE.at(y)),3./2.) / (2 * StateJ.at(x) + 1);
					else
						TransitionAmplitudes.at(i)[y][x] = 0.;
				}
			}
		}
	}

	Lifetimes.ResizeTo(TransitionStrengths.at(0).GetNcols(),TransitionStrengths.at(0).GetNrows());
	SummedTransitionStrengths.ResizeTo(TransitionStrengths.at(0).GetNcols(),TransitionStrengths.at(0).GetNrows());
	for(int x=0;x<TransitionStrengths.at(0).GetNcols();x++){
		for(int y=0;y<TransitionStrengths.at(0).GetNrows();y++){
			SummedTransitionStrengths[y][x] = 0;
			for(unsigned int i=0;i<TransitionStrengths.size();i++){
				SummedTransitionStrengths[y][x] += (TMath::Power(TMath::Abs(StateE.at(x)-StateE.at(y)),power[i]) * TransitionStrengths.at(i)[y][x]) / multfactor[i];
			}
			if(SummedTransitionStrengths[y][x]>0)
				Lifetimes[y][x] = (1 / SummedTransitionStrengths[y][x]);
			else
				Lifetimes[y][x] = 0;
		}
	}

	for(int x=0;x<SummedTransitionStrengths.GetNrows();x++){
		for(int y=x+1;y<SummedTransitionStrengths.GetNcols();y++){
			if(SummedTransitionStrengths[y][x] > 0)
				nDecays++;
		}
	}
	
	StateLifetimes.ResizeTo(SummedTransitionStrengths.GetNcols());	
	StateDecayProb.ResizeTo(SummedTransitionStrengths.GetNcols());
	for(int x=0;x<SummedTransitionStrengths.GetNcols();x++){
		double tmp = 0;
		for(int y=0;y<SummedTransitionStrengths.GetNrows();y++)
			tmp+=SummedTransitionStrengths[y][x] * 1e-12;
		StateDecayProb[x] = tmp;
		if(tmp>0)
			StateLifetimes[x] = 1/tmp;
		else
			StateLifetimes[x] = 0;
	}
		

	BranchingRatios.ResizeTo(fNucleus->GetNstates(),fNucleus->GetNstates());
	for(int x=0;x<BranchingRatios.GetNcols();x++){
		double sum = SumColumn(SummedTransitionStrengths,x);
		if(sum > 0){
			for(int y=0;y<BranchingRatios.GetNrows();y++)
				BranchingRatios[y][x] = SummedTransitionStrengths[y][x] / sum;
		}
		else
			for(int y=0;y<BranchingRatios.GetNrows();y++)
				BranchingRatios[y][x] = 0;

	}
	

	MixingRatios.ResizeTo(fNucleus->GetNstates(),fNucleus->GetNstates());
	for(int x=0;x<MixingRatios.GetNcols();x++){
		for(int y=0;y<MixingRatios.GetNrows();y++){
			if(TransitionStrengths_Abs.at(6)[y][x] != 0 && TransitionStrengths_Abs.at(1)[y][x] !=0)
				MixingRatios[y][x] = TransitionStrengths_Abs.at(1)[y][x] / TransitionStrengths_Abs.at(6)[y][x];
		}
	}

}

void TransitionRates::Print() const{

	std::string mult[7] = {"E1","E2","E3","E4","E5","E6","M1"};
	std::string MEUnits[7] = {"eb^(1/2)","eb","eb^(3/2)","eb^(2)","eb^(5/2)","eb^3","uN"};	
	std::string TSUnits[7] = {"e2fm2","e2fm4","e2fm6","e2fm8","e2fm10","e2fm12","uN2"};	
	
	std::cout << "\nMatrix elements" << std::endl;
	for(unsigned int i=0;i<MatrixElements.size();i++){
		if(MiscFunctions::GetMaxAbsMatrix(MatrixElements.at(i)) == 0)
			continue;
		std::cout << mult[i] << "\t" << MEUnits[i] << std::endl;
		MiscFunctions::PrintMatrixNucleus(MatrixElements.at(i),*fNucleus);
	}

	std::cout << "\nReduced transition strengths" << std::endl;
	for(unsigned int i=0;i<MatrixElements.size();i++){
		if(MiscFunctions::GetMaxAbsMatrix(MatrixElements.at(i)) == 0)
			continue;
		std::cout << mult[i] << "\t" << TSUnits[i] << std::endl;
		MiscFunctions::PrintMatrixNucleus(TransitionStrengths.at(i),*fNucleus);
	}

	std::cout << "\nReduced decay rate" << std::endl;
	for(unsigned int i=0;i<MatrixElements.size();i++){
		if(MiscFunctions::GetMaxAbsMatrix(MatrixElements.at(i)) == 0)
			continue;
		std::cout << mult[i] << "\t/ps" << std::endl;
		MiscFunctions::PrintMatrixNucleus(TransitionStrengths_Abs.at(i) * 1e-12,*fNucleus);
	}

	std::cout << "\nTotal transition strengths (/ps)" << std::endl;
	MiscFunctions::PrintMatrixNucleus(SummedTransitionStrengths * 1e-12,*fNucleus);

	std::cout << "\nEffective lifetimes (ps)" << std::endl;
	MiscFunctions::PrintMatrixNucleus(Lifetimes * 1e12,*fNucleus);

	std::cout << "\nState transition strengths (/ps)" << std::endl;
	MiscFunctions::PrintVectorNucleus(StateDecayProb,*fNucleus,"Transition strength");

	std::cout << "\nState lifetimes (ps)" << std::endl;
	MiscFunctions::PrintVectorNucleus(StateLifetimes,*fNucleus,"Lifetime");

	std::cout << "\nBranching ratios (normalized):" << std::endl;
	MiscFunctions::PrintMatrixNucleus(BranchingRatios,*fNucleus);

	std::cout << "\nMixing ratios:" << std::endl;
	MiscFunctions::PrintMatrixNucleus(MixingRatios,*fNucleus);

}

double TransitionRates::SumColumn(TMatrixD mat, int col) const{
	double sum = 0;
	for(int y=0;y<mat.GetNrows();y++)
		sum += mat[y][col];
	return sum;
}
double TransitionRates::MaxAbsMatrix(TMatrixD mat) const{
	double max=0;
	for(int x=0;x<mat.GetNcols();x++)
		for(int y=0;y<mat.GetNrows();y++)
			if(TMath::Abs(mat[y][x])>0)
				max = mat[y][x];
	return max;
}
