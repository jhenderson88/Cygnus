#include "CoulExSimFitter.h"

CoulExSimFitter::CoulExSimFitter()
{

	ClearAll();

	first		= true;

	fLikelihood	= false;

	fDoFullUnc	= false;

	maxIter		= 500;
	maxCalls	= 500;
	fitTolerance	= 0.001;
	nThreads	= 1;

	verbose		= false;

}

void CoulExSimFitter::DoFit(const char* method, const char *algorithm){

	CoulExSimMinFCN theFCN(exptData_Beam,exptData_Target);

	theFCN.SetBeamMatrixElements(matrixElements_Beam);
	theFCN.SetTargetMatrixElements(matrixElements_Target);
	theFCN.SetScalingParameters(scalingParameters);

	theFCN.SetBeamLitLifetimes(litLifetimes_Beam);
	theFCN.SetBeamLitBranching(litBranchingRatios_Beam);
	theFCN.SetBeamLitMixing(litMixingRatios_Beam);	
	theFCN.SetBeamLitMatrixElements(litMatrixElements_Beam);	
	theFCN.SetTargetLitLifetimes(litLifetimes_Target);
	theFCN.SetTargetLitBranching(litBranchingRatios_Target);
	theFCN.SetTargetLitMixing(litMixingRatios_Target);	
	theFCN.SetTargetLitMatrixElements(litMatrixElements_Target);	

	theFCN.SetBeamPointCalcs(pointCalcs_Beam);
	theFCN.SetTargetPointCalcs(pointCalcs_Target);

	theFCN.SetBeamNucleus(&fNucleus_Beam);
	theFCN.SetTargetNucleus(&fNucleus_Target);

	theFCN.SetBeamCorrectionFactors(correctionFactors_Beam);
	theFCN.SetTargetCorrectionFactors(correctionFactors_Target);

	theFCN.SetIter(maxIter);
	theFCN.SetCalls(maxCalls);

	theFCN.SetNthreads(nThreads);

	theFCN.SetVerbose(verbose);

	theFCN.SetupCalculation();

	theFCN.SetLikelihoodFit(fLikelihood);

	for(size_t m=0;m<matrixElements_Beam.size();m++)
		matrixElements_Beam.at(m).Print();
	for(size_t m=0;m<matrixElements_Target.size();m++)
		matrixElements_Target.at(m).Print();

	parameters.clear();
	par_LL.clear();
	par_UL.clear();
	for(unsigned int i=0;i<matrixElements_Beam.size();i++){
		parameters.push_back(matrixElements_Beam.at(i).GetMatrixElement());
		par_LL.push_back(matrixElements_Beam.at(i).GetMatrixElementLowerLimit());
		par_UL.push_back(matrixElements_Beam.at(i).GetMatrixElementUpperLimit());
	}
	for(unsigned int i=0;i<matrixElements_Target.size();i++){
		parameters.push_back(matrixElements_Target.at(i).GetMatrixElement());
		par_LL.push_back(matrixElements_Target.at(i).GetMatrixElementLowerLimit());
		par_UL.push_back(matrixElements_Target.at(i).GetMatrixElementUpperLimit());
	}
	for(unsigned int i=0;i<scalingParameters.size();i++){
		parameters.push_back(scalingParameters.at(i).GetScalingParameter());
		par_LL.push_back(scalingParameters.at(i).GetScalingLowerLimit());
		par_UL.push_back(scalingParameters.at(i).GetScalingUpperLimit());
	}	

	std::cout 	<< std::setw(12) << std::left << "Parameters:" 
			<< std::endl;
	for(unsigned int i=0;i<matrixElements_Beam.size();i++){
		std::cout	<< std::setw(11) << std::left << "Beam ME:" 
				<< std::setw(4) << std::left << i+1;
	}
	for(unsigned int i=0;i<matrixElements_Target.size();i++){
		std::cout	<< std::setw(11) << std::left << "Target ME:" 
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

	if(fLikelihood)
		min->SetErrorDef(0.5);

	std::cout << "Iterations: " << maxIter << std::endl;
	std::cout << "Calls: " << maxCalls << std::endl;

	min->SetMaxFunctionCalls(maxCalls);
	min->SetMaxIterations(maxIter);
	min->SetTolerance(fitTolerance);
	min->SetFunction(f);
	for(unsigned int i=0; i<parameters.size(); i++){
		std::string name;
		if(i < matrixElements_Beam.size()){
			name = "Beam-ME-"+std::to_string(i);
			min->SetLimitedVariable(i,name,parameters.at(i),0.01,par_LL.at(i),par_UL.at(i));
		}
		else if(i < matrixElements_Target.size() + matrixElements_Beam.size()){
			name = "Target-ME-"+std::to_string(i - matrixElements_Beam.size());
			min->SetLimitedVariable(i,name,parameters.at(i),0.01,par_LL.at(i),par_UL.at(i));
		}
		else{
			name = "Scaling-"+std::to_string(i-(matrixElements_Beam.size() + matrixElements_Target.size()));
			min->SetLowerLimitedVariable(i,name,parameters.at(i),0.01,0);
		}
	}

	if(!verbose && !fLikelihood){
		std::cout 	<< std::setw(12) << std::left << "Iteration:" 
				<< std::setw(13) << std::left << "Chi2 value:" 
				<< std::setw(7)  << std::left << "NDF:"
				<< std::setw(13) << std::left << "Red. Chi2:"
				<< std::setw(12) << std::left << "Lit. Chi2:" 
				<< std::setw(24) << std::left << "Processing time: (ms)" 
				<< std::endl;
	}
	else{
		std::cout 	<< std::setw(12) << std::left << "Iteration:" 
				<< std::setw(13) << std::left << "-Ln(L) value:" 
				<< std::setw(7)  << std::left << "NDF:"
				<< std::setw(13) << std::left << "Red. Ln(L):"
				<< std::setw(12) << std::left << "Lit. -Ln(L):" 
				<< std::setw(24) << std::left << "Processing time: (ms)" 
				<< std::endl;
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

	covMat.ResizeTo(parameters.size(),parameters.size());
	corMat.ResizeTo(parameters.size(),parameters.size());
	for(unsigned int i=0;i<parameters.size();i++){
		for(unsigned int j=i;j<parameters.size();j++){
			covMat[i][j] = min->CovMatrix(i,j);
			covMat[j][i] = min->CovMatrix(j,i);
			corMat[i][j] = min->Correlation(i,j);
			corMat[j][i] = min->Correlation(j,i);
		}
	}

	if(DoFullUncertainty()){
		std::cout	<< "MINOS uncertainties (asymmetric):"
				<< std::endl;
		std::vector<double> errLowVec, errUpVec;
		for(unsigned int i=0;i<parameters.size();i++){
			double errLow, errUp;
			min->GetMinosError(i,errLow,errUp);
			errLowVec.push_back(errLow);
			errUpVec.push_back(errUp);
		}
		std::cout	<< "Correlated uncertainty calculation completed"
				<< std::endl;
		std::cout	<< std::setw(14) << std::left << "Parameter" 
				<< std::setw(14) << std::left << "Value" 
				<< std::setw(14) << "+" 
				<< std::setw(3) << "/"
				<< std::setw(14) << "-" 
				<< std::endl;
		for(unsigned int i=0;i<parameters.size();i++){
			std::cout	<< std::setw(14) << std::left << min->VariableName(i)
					<< std::setw(14) << std::left << parameters[i]
					<< std::setw(14) << errUpVec[i] 
					<< std::setw(3) << ""
					<< std::setw(14) << errLowVec[i]
					<< std::endl;
		}
	}

}

void CoulExSimFitter::CreateScalingParameter(std::vector<int> expnum, double scaling, double scaling_LL, double scaling_UL){

	ScalingParameter tmpScaling;
	tmpScaling.SetExperimentVector(expnum);
	tmpScaling.SetScalingValue(scaling,scaling_LL,scaling_UL);

	scalingParameters.push_back(tmpScaling);

}

void CoulExSimFitter::AddBeamFittingMatrixElement(int lambda, int init, int fin, double ME, double LL, double UL){
	MatrixElement tmpME(matrixElements_Beam.size(),lambda,init,fin,ME,LL,UL);
	matrixElements_Beam.push_back(tmpME);
}
void CoulExSimFitter::AddTargetFittingMatrixElement(int lambda, int init, int fin, double ME, double LL, double UL){
	MatrixElement tmpME(matrixElements_Target.size(),lambda,init,fin,ME,LL,UL);
	matrixElements_Target.push_back(tmpME);
}

void CoulExSimFitter::AddBeamCorrectionFactor(TVectorD corrFac){
	correctionFactors_Beam.push_back(corrFac);
}
void CoulExSimFitter::SetBeamCorrectionFactor(int i, TVectorD corrFac){
	if((i < (int)correctionFactors_Beam.size()))
		correctionFactors_Beam.at(i) = corrFac;
	else
		std::cout << "Outside vector range" << std::endl;
}
void CoulExSimFitter::AddTargetCorrectionFactor(TVectorD corrFac){
	correctionFactors_Target.push_back(corrFac);
}
void CoulExSimFitter::SetTargetCorrectionFactor(int i, TVectorD corrFac){
	if((i < (int)correctionFactors_Target.size()))
		correctionFactors_Target.at(i) = corrFac;
	else
		std::cout << "Outside vector range" << std::endl;
}

void CoulExSimFitter::DefineExperiment(double thetacm){
	ExperimentData tmpExp;
	tmpExp.SetThetaCM(thetacm);
	exptData_Beam.push_back(tmpExp);		
	exptData_Target.push_back(tmpExp);		
}
void CoulExSimFitter::AddBeamData(int nExpt, int init, int fin, double counts, double unc){
	exptData_Beam.at(nExpt).AddData(init,fin,counts,unc);
} 
void CoulExSimFitter::AddTargetData(int nExpt, int init, int fin, double counts, double unc){
	exptData_Target.at(nExpt).AddData(init,fin,counts,unc);
}   

void CoulExSimFitter::AddBeamLifetime(int index, double lifetime, double unc){
	LitLifetime tmpLifetime(index,lifetime,unc);
	litLifetimes_Beam.push_back(tmpLifetime);	
}
void CoulExSimFitter::AddBeamBranchingRatio(int index_I1, int index_F1, int index_F2, double br, double unc){
	LitBranchingRatio tmpBR(index_I1,index_F1,index_F2,br,unc);
	litBranchingRatios_Beam.push_back(tmpBR);	
}  
void CoulExSimFitter::AddBeamMixingRatio(int index_I, int index_F, double delta, double unc){
	LitMixingRatio tmpMR(index_I,index_F,delta,unc);
	litMixingRatios_Beam.push_back(tmpMR);	
}
void CoulExSimFitter::AddBeamMatrixElement(int mult, int index_I, int index_F, double me, double unc){
	LitMatrixElement tmpME(mult,index_I,index_F,me,unc);
	litMatrixElements_Beam.push_back(tmpME);	
}
void CoulExSimFitter::AddTargetLifetime(int index, double lifetime, double unc){
	LitLifetime tmpLifetime(index,lifetime,unc);
	litLifetimes_Target.push_back(tmpLifetime);	
}
void CoulExSimFitter::AddTargetBranchingRatio(int index_I1, int index_F1, int index_F2, double br, double unc){
	LitBranchingRatio tmpBR(index_I1,index_F1,index_F2,br,unc);
	litBranchingRatios_Target.push_back(tmpBR);	
}  
void CoulExSimFitter::AddTargetMixingRatio(int index_I, int index_F, double delta, double unc){
	LitMixingRatio tmpMR(index_I,index_F,delta,unc);
	litMixingRatios_Target.push_back(tmpMR);	
}    
void CoulExSimFitter::AddTargetMatrixElement(int mult, int index_I, int index_F, double me, double unc){
	LitMatrixElement tmpME(mult, index_I,index_F,me,unc);
	litMatrixElements_Target.push_back(tmpME);	
}    

void CoulExSimFitter::ClearAll(){

	index.clear();
	parameters.clear();
			
	matrixElements_Beam.clear();			
	pointCalcs_Beam.clear();			
	exptData_Beam.clear();			
	litLifetimes_Beam.clear();			
	litBranchingRatios_Beam.clear();		
	litMixingRatios_Beam.clear();		
	litMatrixElements_Beam.clear();		
	EffectiveCrossSection_Beam.clear();
			
	matrixElements_Target.clear();			
	pointCalcs_Target.clear();			
	exptData_Target.clear();			
	litLifetimes_Target.clear();			
	litBranchingRatios_Target.clear();		
	litMixingRatios_Target.clear();		
	litMatrixElements_Target.clear();		
	EffectiveCrossSection_Target.clear();

}

void CoulExSimFitter::Print() const{

	if(exptData_Beam.size()>0){
		std::cout 	<< "\n\n"
				<< "Experimental data (beam):" << std::endl;
		
		std::cout 	<< exptData_Beam.size() << " experiments" << std::endl;
		
		for(unsigned int i=0;i<exptData_Beam.size();i++){
			std::cout	<< "Experiment " << i+1 << std::endl;
			std::cout	<< "Theta [CM]: " << exptData_Beam.at(i).GetThetaCM() << std::endl;
			std::cout 	<< std::setw(15) << std::left << "Init. index:" 
					<< std::setw(15) << std::left << "Final index:"
					<< std::setw(15) << std::left << "Init. J:" 
					<< std::setw(15) << std::left << "Final J:"
					<< std::setw(10) << std::left << "Counts:"
					<< std::setw(10) << std::left << "Unc:"
					<< std::endl;
			for(unsigned int t=0;t<exptData_Beam.at(i).GetData().size();t++){
				std::cout 	<< std::setw(15) << std::left << exptData_Beam.at(i).GetDataPoint(t).GetInitialIndex()
						<< std::setw(15) << std::left << exptData_Beam.at(i).GetDataPoint(t).GetFinalIndex()
						<< std::setw(15) << std::left << fNucleus_Beam.GetLevelJ().at(exptData_Beam.at(i).GetDataPoint(t).GetInitialIndex())
						<< std::setw(15) << std::left << fNucleus_Beam.GetLevelJ().at(exptData_Beam.at(i).GetDataPoint(t).GetFinalIndex())
						<< std::setw(10) << std::left << exptData_Beam.at(i).GetDataPoint(t).GetCounts()
						<< std::setw(10) << std::left << exptData_Beam.at(i).GetDataPoint(t).GetUpUnc()
						<< std::endl;
			}			
		}
	}
	else	
		std::cout << "No experimental data declared" << std::endl;

	if(litLifetimes_Beam.size()>0){
		std::cout	<< "\n\n"
				<< "Literature lifetimes (beam):"
				<< std::endl;
		std::cout	<< std::setw(8)  << std::left << "Index"
				<< std::setw(6)  << std::left << "J:"
				<< std::setw(15) << std::left << "Lifetime (ps)" 
				<< std::setw(15) << std::left << "Uncertainty:"
				<< std::endl;
		for(unsigned int i=0;i<litLifetimes_Beam.size();i++){
			std::cout 	<< std::setw(8)  << std::left << litLifetimes_Beam.at(i).GetIndex()
					<< std::setw(6)  << std::left << fNucleus_Beam.GetLevelJ().at(litLifetimes_Beam.at(i).GetIndex())
					<< std::setw(15) << std::left << litLifetimes_Beam.at(i).GetLifetime()
					<< std::setw(15) << std::left << litLifetimes_Beam.at(i).GetUpUnc()
					<< std::endl;
		}
	}

	if(litBranchingRatios_Beam.size()>0){
		std::cout	<< "\n\n"
				<< "Literature Branching Ratios (beam):"
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
		for(unsigned int i=0;i<litBranchingRatios_Beam.size();i++){
			std::cout 	<< std::setw(15) << std::left << litBranchingRatios_Beam.at(i).GetInitialIndex()
					<< std::setw(15) << std::left << litBranchingRatios_Beam.at(i).GetFinalIndex_1()
					<< std::setw(15) << std::left << litBranchingRatios_Beam.at(i).GetFinalIndex_2()
					<< std::setw(10) << std::left << fNucleus_Beam.GetLevelJ().at(litBranchingRatios_Beam.at(i).GetInitialIndex())
					<< std::setw(12) << std::left << fNucleus_Beam.GetLevelJ().at(litBranchingRatios_Beam.at(i).GetFinalIndex_1())
					<< std::setw(12) << std::left << fNucleus_Beam.GetLevelJ().at(litBranchingRatios_Beam.at(i).GetFinalIndex_2())
					<< std::setw(17) << std::left << litBranchingRatios_Beam.at(i).GetBranchingRatio()
					<< std::setw(15) << std::left << litBranchingRatios_Beam.at(i).GetUpUnc()
					<< std::endl;
		}
	}
	
	if(litMixingRatios_Beam.size()>0){
		std::cout	<< "\n\n"
				<< "Literature Mixing Ratios (beam):"
				<< std::endl;
		std::cout	<< std::setw(15) << std::left << "Init. Index"
				<< std::setw(15) << std::left << "Final Index"
				<< std::setw(10) << std::left << "J init:"
				<< std::setw(10) << std::left << "J final:"
				<< std::setw(14) << std::left << "Mixing Ratio" 
				<< std::setw(15) << std::left << "Uncertainty:"
				<< std::endl;
		for(unsigned int i=0;i<litMixingRatios_Beam.size();i++){
			std::cout 	<< std::setw(15) << std::left << litMixingRatios_Beam.at(i).GetInitialIndex()
					<< std::setw(15) << std::left << litMixingRatios_Beam.at(i).GetFinalIndex()
					<< std::setw(10) << std::left << fNucleus_Beam.GetLevelJ().at(litMixingRatios_Beam.at(i).GetInitialIndex())
					<< std::setw(10) << std::left << fNucleus_Beam.GetLevelJ().at(litMixingRatios_Beam.at(i).GetFinalIndex())
					<< std::setw(14) << std::left << litMixingRatios_Beam.at(i).GetMixingRatio()
					<< std::setw(15) << std::left << litMixingRatios_Beam.at(i).GetUpUnc()
					<< std::endl;
		}
	}

	std::string	mult[7] = {"E1","E2","E3","E4","E5","E6","M1"};

	std::cout 	<< "\n\n"
			<< "Starting matrix elements (beam):"
			<< std::endl;
	for(unsigned int l = 0; l < fNucleus_Beam.GetMatrixElements().size(); l++){
		if(MiscFunctions::GetMaxMatrix(fNucleus_Beam.GetMatrixElements().at(l)) > 0){
			std::cout	<< mult[l] << " matrix elements"
					<< std::endl;	
			TMatrixD	tmpMat;
			tmpMat.ResizeTo(fNucleus_Beam.GetMatrixElements().at(l).GetNrows(),fNucleus_Beam.GetMatrixElements().at(l).GetNcols());
			tmpMat =	fNucleus_Beam.GetMatrixElements().at(l);
			for(unsigned int s = 0; s < matrixElements_Beam.size(); s++){
				if(matrixElements_Beam.at(s).GetLambda() == (int)l){
					tmpMat[matrixElements_Beam.at(s).GetInitialState()][matrixElements_Beam.at(s).GetFinalState()] = matrixElements_Beam.at(s).GetMatrixElement();
					tmpMat[matrixElements_Beam.at(s).GetFinalState()][matrixElements_Beam.at(s).GetInitialState()] = matrixElements_Beam.at(s).GetMatrixElement();
				}
			}
			MiscFunctions::PrintMatrixNucleus(tmpMat,fNucleus_Beam);	

		}
	}	

	if(exptData_Target.size()>0){
		std::cout 	<< "\n\n"
				<< "Experimental data (target):" << std::endl;
		
		std::cout 	<< exptData_Target.size() << " experiments" << std::endl;
		
		for(unsigned int i=0;i<exptData_Target.size();i++){
			std::cout	<< "Experiment " << i+1 << std::endl;
			std::cout	<< "Theta [CM]: " << exptData_Target.at(i).GetThetaCM() << std::endl;
			std::cout 	<< std::setw(15) << std::left << "Init. index:" 
					<< std::setw(15) << std::left << "Final index:"
					<< std::setw(15) << std::left << "Init. J:" 
					<< std::setw(15) << std::left << "Final J:"
					<< std::setw(10) << std::left << "Counts:"
					<< std::setw(10) << std::left << "Unc:"
					<< std::endl;
			for(unsigned int t=0;t<exptData_Target.at(i).GetData().size();t++){
				std::cout 	<< std::setw(15) << std::left << exptData_Target.at(i).GetDataPoint(t).GetInitialIndex()
						<< std::setw(15) << std::left << exptData_Target.at(i).GetDataPoint(t).GetFinalIndex()
						<< std::setw(15) << std::left << fNucleus_Target.GetLevelJ().at(exptData_Target.at(i).GetDataPoint(t).GetInitialIndex())
						<< std::setw(15) << std::left << fNucleus_Target.GetLevelJ().at(exptData_Target.at(i).GetDataPoint(t).GetFinalIndex())
						<< std::setw(10) << std::left << exptData_Target.at(i).GetDataPoint(t).GetCounts()
						<< std::setw(10) << std::left << exptData_Target.at(i).GetDataPoint(t).GetUpUnc()
						<< std::endl;
			}			
		}
	}
	else	
		std::cout << "No experimental data declared" << std::endl;

	if(litLifetimes_Target.size()>0){
		std::cout	<< "\n\n"
				<< "Literature lifetimes (target):"
				<< std::endl;
		std::cout	<< std::setw(8)  << std::left << "Index"
				<< std::setw(6)  << std::left << "J:"
				<< std::setw(15) << std::left << "Lifetime (ps)" 
				<< std::setw(15) << std::left << "Uncertainty:"
				<< std::endl;
		for(unsigned int i=0;i<litLifetimes_Target.size();i++){
			std::cout 	<< std::setw(8)  << std::left << litLifetimes_Target.at(i).GetIndex()
					<< std::setw(6)  << std::left << fNucleus_Target.GetLevelJ().at(litLifetimes_Target.at(i).GetIndex())
					<< std::setw(15) << std::left << litLifetimes_Target.at(i).GetLifetime()
					<< std::setw(15) << std::left << litLifetimes_Target.at(i).GetUpUnc()
					<< std::endl;
		}
	}

	if(litBranchingRatios_Target.size()>0){
		std::cout	<< "\n\n"
				<< "Literature Branching Ratios (target):"
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
		for(unsigned int i=0;i<litBranchingRatios_Target.size();i++){
			std::cout 	<< std::setw(15) << std::left << litBranchingRatios_Target.at(i).GetInitialIndex()
					<< std::setw(15) << std::left << litBranchingRatios_Target.at(i).GetFinalIndex_1()
					<< std::setw(15) << std::left << litBranchingRatios_Target.at(i).GetFinalIndex_2()
					<< std::setw(10) << std::left << fNucleus_Target.GetLevelJ().at(litBranchingRatios_Target.at(i).GetInitialIndex())
					<< std::setw(12) << std::left << fNucleus_Target.GetLevelJ().at(litBranchingRatios_Target.at(i).GetFinalIndex_1())
					<< std::setw(12) << std::left << fNucleus_Target.GetLevelJ().at(litBranchingRatios_Target.at(i).GetFinalIndex_2())
					<< std::setw(17) << std::left << litBranchingRatios_Target.at(i).GetBranchingRatio()
					<< std::setw(15) << std::left << litBranchingRatios_Target.at(i).GetUpUnc()
					<< std::endl;
		}
	}
	
	if(litMixingRatios_Target.size()>0){
		std::cout	<< "\n\n"
				<< "Literature Mixing Ratios (target):"
				<< std::endl;
		std::cout	<< std::setw(15) << std::left << "Init. Index"
				<< std::setw(15) << std::left << "Final Index"
				<< std::setw(10) << std::left << "J init:"
				<< std::setw(10) << std::left << "J final:"
				<< std::setw(14) << std::left << "Mixing Ratio" 
				<< std::setw(15) << std::left << "Uncertainty:"
				<< std::endl;
		for(unsigned int i=0;i<litMixingRatios_Target.size();i++){
			std::cout 	<< std::setw(15) << std::left << litMixingRatios_Target.at(i).GetInitialIndex()
					<< std::setw(15) << std::left << litMixingRatios_Target.at(i).GetFinalIndex()
					<< std::setw(10) << std::left << fNucleus_Target.GetLevelJ().at(litMixingRatios_Target.at(i).GetInitialIndex())
					<< std::setw(10) << std::left << fNucleus_Target.GetLevelJ().at(litMixingRatios_Target.at(i).GetFinalIndex())
					<< std::setw(14) << std::left << litMixingRatios_Target.at(i).GetMixingRatio()
					<< std::setw(15) << std::left << litMixingRatios_Target.at(i).GetUpUnc()
					<< std::endl;
		}
	}

	std::cout 	<< "\n\n"
			<< "Starting matrix elements (target):"
			<< std::endl;
	for(unsigned int l = 0; l < fNucleus_Target.GetMatrixElements().size(); l++){
		if(MiscFunctions::GetMaxMatrix(fNucleus_Target.GetMatrixElements().at(l)) > 0){
			std::cout	<< mult[l] << " matrix elements"
					<< std::endl;	
			TMatrixD	tmpMat;
			tmpMat.ResizeTo(fNucleus_Target.GetMatrixElements().at(l).GetNrows(),fNucleus_Target.GetMatrixElements().at(l).GetNcols());
			tmpMat =	fNucleus_Target.GetMatrixElements().at(l);
			for(unsigned int s = 0; s < matrixElements_Target.size(); s++){
				if(matrixElements_Target.at(s).GetLambda() == (int)l){
					tmpMat[matrixElements_Target.at(s).GetInitialState()][matrixElements_Target.at(s).GetFinalState()] = matrixElements_Target.at(s).GetMatrixElement();
					tmpMat[matrixElements_Target.at(s).GetFinalState()][matrixElements_Target.at(s).GetInitialState()] = matrixElements_Target.at(s).GetMatrixElement();
				}
			}
			MiscFunctions::PrintMatrixNucleus(tmpMat,fNucleus_Target);	

		}
	}	

}
