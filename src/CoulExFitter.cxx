#include "CoulExFitter.h"


CoulExFitter::CoulExFitter()
{

	ClearAll();

	first		= true;

	maxIter		= 500;
	maxCalls	= 500;
	fitTolerance	= 0.001;
	nThreads	= 1;

	verbose		= false;

}

void CoulExFitter::DoFit(const char* method, const char *algorithm){

	CoulExMinFCN theFCN(exptData);

	theFCN.SetMatrixElements(matrixElements);
	theFCN.SetScalingParameters(scalingParameters);

	theFCN.SetLitLifetimes(litLifetimes);
	theFCN.SetLitBranching(litBranchingRatios);
	theFCN.SetLitMixing(litMixingRatios);	

	theFCN.SetPointCalcs(pointCalcs);

	theFCN.SetNucleus(&fNucleus);

	theFCN.SetCorrectionFactors(correctionFactors);

	theFCN.SetIter(maxIter);
	theFCN.SetCalls(maxCalls);

	theFCN.SetNthreads(nThreads);

	theFCN.SetVerbose(verbose);

	theFCN.SetupCalculation();

	//parameters.push_back(0.15);
	//par_LL.push_back(0.05);
	//par_UL.push_back(0.3);
	//first = false;

	for(size_t m=0;m<matrixElements.size();m++)
		matrixElements.at(m).Print();

	parameters.clear();
	par_LL.clear();
	par_UL.clear();
	for(unsigned int i=0;i<matrixElements.size();i++){
		parameters.push_back(matrixElements.at(i).GetMatrixElement());
		par_LL.push_back(matrixElements.at(i).GetMatrixElementLowerLimit());
		par_UL.push_back(matrixElements.at(i).GetMatrixElementUpperLimit());
	}
	for(unsigned int i=0;i<scalingParameters.size();i++){
		parameters.push_back(scalingParameters.at(i).GetScalingParameter());
		par_LL.push_back(scalingParameters.at(i).GetScalingLowerLimit());
		par_UL.push_back(scalingParameters.at(i).GetScalingUpperLimit());
	}	

	std::cout 	<< std::setw(12) << std::left << "Parameters:" 
			<< std::endl;
	for(unsigned int i=0;i<matrixElements.size();i++){
		std::cout	<< std::setw(11) << std::left << "Matrix El." 
				<< std::setw(4) << std::left << i+1;
	}
	for(unsigned int i=0;i<scalingParameters.size();i++){
		std::cout	<< std::setw(11) << std::left << "Scaling " 
				<< std::setw(4) << std::left << i+1;
	}
	std::cout	<< std::endl;
	for(unsigned int i=0;i<parameters.size();i++){
		std::cout 	<< std::setw(11) << std::left << parameters.at(i) 
				<< std::setw(4) << std::left << "";
	}
	std::cout << std::endl;

	theFCN.SetNpar(parameters.size());

	ROOT::Math::Minimizer *min =
			ROOT::Math::Factory::CreateMinimizer(method, algorithm);
	//		ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");
	//		ROOT::Math::Factory::CreateMinimizer("Minuit2","Simplex");
	//		ROOT::Math::Factory::CreateMinimizer("Minuit2","Combined");
	//		ROOT::Math::Factory::CreateMinimizer("Minuit2","Scan");
	//		ROOT::Math::Factory::CreateMinimizer("Minuit2","Fumili");
	//		ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "ConjugateFR");
	//		ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "ConjugatePR");
	//		ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "BFGS");
	//		ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "BFGS2");
	//		ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "SteepestDescent");
	ROOT::Math::Functor f(theFCN,parameters.size());

	std::cout << "Iterations: " << maxIter << std::endl;
	std::cout << "Calls: " << maxCalls << std::endl;

	min->SetMaxFunctionCalls(maxCalls);
	min->SetMaxIterations(maxIter);
	min->SetTolerance(fitTolerance);
	min->SetFunction(f);
	//min->SetPrintLevel(printLevel);
	for(unsigned int i=0; i<parameters.size(); i++){
		std::string name;
		if(i < matrixElements.size()){
			name = "ME-"+std::to_string(i);
			min->SetLimitedVariable(i,name,parameters.at(i),0.01,par_LL.at(i),par_UL.at(i));
		}
		else{
			name = "Scaling-"+std::to_string(i-matrixElements.size());
			min->SetLowerLimitedVariable(i,name,parameters.at(i),0.0001,0);//,par_LL.at(i),par_UL.at(i));
		}
	}
	
	typedef std::chrono::high_resolution_clock Clock;
	typedef std::chrono::milliseconds milliseconds;
	Clock::time_point t0 = Clock::now();

	min->Minimize();

	Clock::time_point t1 = Clock::now();
	milliseconds ms = std::chrono::duration_cast<milliseconds>(t1-t0);

	std::cout << std::endl;

	min->PrintResults();

	std::cout << "Fitting time: " << ms.count() << " ms" <<	std::endl;
	
	const double	*res = min->X();
	for(unsigned int i=0;i<parameters.size();i++)
		parameters[i] = res[i];

	//for(unsigned int i=0;i<parameters.size();i++)
	//	parameters[i] = theFCN.GetParameter(i);

}

void CoulExFitter::CreateScalingParameter(std::vector<int> expnum, double scaling, double scaling_LL, double scaling_UL){

	ScalingParameter tmpScaling;
	tmpScaling.SetExperimentVector(expnum);
	tmpScaling.SetScalingValue(scaling,scaling_LL,scaling_UL);

	scalingParameters.push_back(tmpScaling);

}

void CoulExFitter::AddFittingMatrixElement(int lambda, int init, int fin, double ME, double LL, double UL){
	//parameters.push_back(ME);
	//par_LL.push_back(LL);
	//par_UL.push_back(UL);
	MatrixElement tmpME(matrixElements.size(),lambda,init,fin,ME,LL,UL);
	matrixElements.push_back(tmpME);
}

void CoulExFitter::AddCorrectionFactor(TVectorD corrFac){
	correctionFactors.push_back(corrFac);
}
void CoulExFitter::SetCorrectionFactor(int i, TVectorD corrFac){
	if((i < (int)correctionFactors.size()))
		correctionFactors.at(i) = corrFac;
	else
		std::cout << "Outside vector range" << std::endl;
}

void CoulExFitter::DefineExperiment(double thetacm){
	ExperimentData tmpExp;
	tmpExp.SetThetaCM(thetacm);
	exptData.push_back(tmpExp);		
}
void CoulExFitter::AddData(int nExpt, int init, int fin, double counts, double unc){
	exptData.at(nExpt).AddData(init,fin,counts,unc);
}  

void CoulExFitter::AddLifetime(int index, double lifetime, double unc){
	LitLifetime tmpLifetime(index,lifetime,unc);
	litLifetimes.push_back(tmpLifetime);	
}
void CoulExFitter::AddBranchingRatio(int index_I1, int index_F1, int index_F2, double br, double unc){
	LitBranchingRatio tmpBR(index_I1,index_F1,index_F2,br,unc);
	litBranchingRatios.push_back(tmpBR);	
}  
void CoulExFitter::AddMixingRatio(int index_I, int index_F, double delta, double unc){
	LitMixingRatio tmpMR(index_I,index_F,delta,unc);
	litMixingRatios.push_back(tmpMR);	
}    

void CoulExFitter::ClearAll(){

	index.clear();
	parameters.clear();			
	matrixElements.clear();			
	pointCalcs.clear();			
	exptData.clear();			
	litLifetimes.clear();			
	litBranchingRatios.clear();		
	litMixingRatios.clear();		
	EffectiveCrossSection.clear();

}

void CoulExFitter::Print() const{

	if(exptData.size()>0){
		std::cout 	<< "\n\n"
				<< "Experimental data:" << std::endl;
		
		std::cout 	<< exptData.size() << " experiments" << std::endl;
		
		for(unsigned int i=0;i<exptData.size();i++){
			std::cout	<< "Experiment " << i+1 << std::endl;
			std::cout	<< "Theta [CM]: " << exptData.at(i).GetThetaCM() << std::endl;
			std::cout 	<< std::setw(15) << std::left << "Init. index:" 
					<< std::setw(15) << std::left << "Final index:"
					<< std::setw(15) << std::left << "Init. J:" 
					<< std::setw(15) << std::left << "Final J:"
					<< std::setw(10) << std::left << "Counts:"
					<< std::setw(10) << std::left << "Unc:"
					<< std::endl;
			for(unsigned int t=0;t<exptData.at(i).GetData().size();t++){
				std::cout 	<< std::setw(15) << std::left << exptData.at(i).GetDataPoint(t).GetInitialIndex()
						<< std::setw(15) << std::left << exptData.at(i).GetDataPoint(t).GetFinalIndex()
						<< std::setw(15) << std::left << fNucleus.GetLevelJ().at(exptData.at(i).GetDataPoint(t).GetInitialIndex())
						<< std::setw(15) << std::left << fNucleus.GetLevelJ().at(exptData.at(i).GetDataPoint(t).GetFinalIndex())
						<< std::setw(10) << std::left << exptData.at(i).GetDataPoint(t).GetCounts()
						<< std::setw(10) << std::left << exptData.at(i).GetDataPoint(t).GetUpUnc()
						<< std::endl;
			}			
		}
	}
	else	
		std::cout << "No experimental data declared" << std::endl;

	if(litLifetimes.size()>0){
		std::cout	<< "\n\n"
				<< "Literature lifetimes:"
				<< std::endl;
		std::cout	<< std::setw(8)  << std::left << "Index"
				<< std::setw(6)  << std::left << "J:"
				<< std::setw(15) << std::left << "Lifetime (ps)" 
				<< std::setw(15) << std::left << "Uncertainty:"
				<< std::endl;
		for(unsigned int i=0;i<litLifetimes.size();i++){
			std::cout 	<< std::setw(8)  << std::left << litLifetimes.at(i).GetIndex()
					<< std::setw(6)  << std::left << fNucleus.GetLevelJ().at(litLifetimes.at(i).GetIndex())
					<< std::setw(15) << std::left << litLifetimes.at(i).GetLifetime()
					<< std::setw(15) << std::left << litLifetimes.at(i).GetUpUnc()
					<< std::endl;
		}
	}

	if(litBranchingRatios.size()>0){
		std::cout	<< "\n\n"
				<< "Literature Branching Ratios:"
				<< std::endl;
		std::cout	<< std::setw(15) << std::left << "Init. Index"
				<< std::setw(15) << std::left << "Final Index 1"
				<< std::setw(15) << std::left << "Final Index 2"
				<< std::setw(10) << std::left << "J init:"
				<< std::setw(12) << std::left << "J final 1:"
				<< std::setw(12) << std::left << "J final 2:"
				<< std::setw(17) << std::left << "Branching Ratio" 
				<< std::setw(15) << std::left << "Uncertainty:"
				<< std::endl;
		for(unsigned int i=0;i<litBranchingRatios.size();i++){
			std::cout 	<< std::setw(15) << std::left << litBranchingRatios.at(i).GetInitialIndex()
					<< std::setw(15) << std::left << litBranchingRatios.at(i).GetFinalIndex_1()
					<< std::setw(15) << std::left << litBranchingRatios.at(i).GetFinalIndex_2()
					<< std::setw(10) << std::left << fNucleus.GetLevelJ().at(litBranchingRatios.at(i).GetInitialIndex())
					<< std::setw(12) << std::left << fNucleus.GetLevelJ().at(litBranchingRatios.at(i).GetFinalIndex_1())
					<< std::setw(12) << std::left << fNucleus.GetLevelJ().at(litBranchingRatios.at(i).GetFinalIndex_2())
					<< std::setw(17) << std::left << litBranchingRatios.at(i).GetBranchingRatio()
					<< std::setw(15) << std::left << litBranchingRatios.at(i).GetUpUnc()
					<< std::endl;
		}
	}
	
	if(litMixingRatios.size()>0){
		std::cout	<< "\n\n"
				<< "Literature Mixing Ratios:"
				<< std::endl;
		std::cout	<< std::setw(15) << std::left << "Init. Index"
				<< std::setw(15) << std::left << "Final Index"
				<< std::setw(10) << std::left << "J init:"
				<< std::setw(10) << std::left << "J final:"
				<< std::setw(14) << std::left << "Mixing Ratio" 
				<< std::setw(15) << std::left << "Uncertainty:"
				<< std::endl;
		for(unsigned int i=0;i<litMixingRatios.size();i++){
			std::cout 	<< std::setw(15) << std::left << litMixingRatios.at(i).GetInitialIndex()
					<< std::setw(15) << std::left << litMixingRatios.at(i).GetFinalIndex()
					<< std::setw(10) << std::left << fNucleus.GetLevelJ().at(litMixingRatios.at(i).GetInitialIndex())
					<< std::setw(10) << std::left << fNucleus.GetLevelJ().at(litMixingRatios.at(i).GetFinalIndex())
					<< std::setw(14) << std::left << litMixingRatios.at(i).GetMixingRatio()
					<< std::setw(15) << std::left << litMixingRatios.at(i).GetUpUnc()
					<< std::endl;
		}
	}

	std::string	mult[7] = {"E1","E2","E3","E4","E5","E6","M1"};

	std::cout 	<< "\n\n"
			<< "Starting matrix elements:"
			<< std::endl;
	for(unsigned int l = 0; l < fNucleus.GetMatrixElements().size(); l++){
		if(MiscFunctions::GetMaxMatrix(fNucleus.GetMatrixElements().at(l)) > 0){
			std::cout	<< mult[l] << " matrix elements"
					<< std::endl;	
			TMatrixD	tmpMat;
			tmpMat.ResizeTo(fNucleus.GetMatrixElements().at(l).GetNrows(),fNucleus.GetMatrixElements().at(l).GetNcols());
			tmpMat =	fNucleus.GetMatrixElements().at(l);
			for(unsigned int s = 0; s < matrixElements.size(); s++){
				if(matrixElements.at(s).GetLambda() == (int)l){
					tmpMat[matrixElements.at(s).GetInitialState()][matrixElements.at(s).GetFinalState()] = matrixElements.at(s).GetMatrixElement();
					tmpMat[matrixElements.at(s).GetFinalState()][matrixElements.at(s).GetInitialState()] = matrixElements.at(s).GetMatrixElement();
				}
			}
			MiscFunctions::PrintMatrixNucleus(tmpMat,fNucleus);	

		}
	}	

}
