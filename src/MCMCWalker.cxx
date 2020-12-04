#include "MCMCWalker.h"

MCMCWalker::MCMCWalker()
{

	ClearAll();

	SetPoissonUncertainties(false);
	SetLikelihoodFit(false);
	SetDoFullUncertainty(false);

	SetMaxIterations(500);
	SetMaxFunctionCalls(500);
	SetTolerance(0.001);
	SetNthreads(1);
	SetVerbose(false);

	ScalingSigma.clear();
	MatrixElementSigma.clear();

}

void MCMCWalker::CreateScalingParameter(std::vector<int> expnum){

	ScalingParameter tmpScaling;
	tmpScaling.SetExperimentVector(expnum);

}

void MCMCWalker::AddFittingMatrixElement(int lambda, int init, int fin, double ME, double LL, double UL, double sigma){
	MatrixElement tmpME(GetMatrixElements().size(),lambda,init,fin,ME,LL,UL);
	AddMatrixElement(tmpME);
	MatrixElementSigma.push_back(sigma);
}

void MCMCWalker::RandomizeStarting(){

	TRandom3 r(0);
	std::vector<MatrixElement> ME = GetMatrixElements();
	for(size_t m = 0; m < ME.size(); m++)
		ME.at(m).SetMatrixElement(r.Rndm()*(ME.at(m).GetMatrixElementUpperLimit() - ME.at(m).GetMatrixElementLowerLimit() - 2*MatrixElementSigma.at(m)) + ME.at(m).GetMatrixElementLowerLimit() + MatrixElementSigma.at(m));
	SetMatrixElements(ME);

}

void MCMCWalker::DoMCMCFit(){

	typedef std::chrono::high_resolution_clock Clock;
	typedef std::chrono::seconds seconds;

	CoulExMinFCN theFCN(GetData());

	theFCN.SetMatrixElements(GetMatrixElements());
	theFCN.SetScalingParameters(GetScalingParameters());

	theFCN.SetLitLifetimes(GetLitLifetimes());
	theFCN.SetLitBranching(GetLitBranching());
	theFCN.SetLitMixing(GetLitMixing());
	theFCN.SetLitMatrixElements(GetLitMatrixElements());	

	theFCN.SetPointCalcs(GetPointCalcs());

	Nucleus	tmpNucleus = GetNucleus();
	theFCN.SetNucleus(&tmpNucleus);

	theFCN.SetCorrectionFactors(GetCorrectionFactors());

	theFCN.SetIter(GetMaxIterations());
	theFCN.SetCalls(GetMaxFunctionCalls());

	theFCN.SetNthreads(GetNthreads());

	theFCN.SetVerbose(false);

	theFCN.SetupCalculation();

	theFCN.SetPoisson(false);

	theFCN.SetLikelihoodFit(true);

	std::vector<double>	fitPar;
	std::vector<double>	fitUL;
	std::vector<double>	fitLL;
	std::vector<double>	fitSigma;

	for(unsigned int i=0;i<GetMatrixElements().size();i++){
		fitPar.push_back(GetMatrixElements().at(i).GetMatrixElement());
		fitLL.push_back(GetMatrixElements().at(i).GetMatrixElementLowerLimit());
		fitUL.push_back(GetMatrixElements().at(i).GetMatrixElementUpperLimit());
		fitSigma.push_back(MatrixElementSigma.at(i));
	}

	SetFitParameters(fitPar);
	SetFitLL(fitLL);
	SetFitUL(fitUL);	

	theFCN.SetNpar(fitPar.size());

	parTracker.clear();
	parTracker.push_back(fitPar);

	std::vector<double>	curr_Par = fitPar;
	std::vector<double>	prop_Par;
	prop_Par.resize(fitPar.size());

	sampler = TRandom3(0);
	rand	= TRandom3(0);
		
	double	*curr_arr = curr_Par.data();	// C++11

	double	curr_Prob	= -theFCN(curr_arr);
	probTracker.push_back(curr_Prob);

	if(GetVerbose()){	
		std::cout 	<< std::setw(12) << std::left << "Iteration:" 
				<< std::setw(15) << std::left << "Accepted (%):"
				<< std::setw(15) << std::left << "Run time (s):"
				<< std::endl;
	}

	std::vector<int>	accepted;

	Clock::time_point t0 = Clock::now();
	Clock::time_point t1;

	for(int c = 1; c < GetMaxFunctionCalls(); c++){

		for(size_t s = 0; s<GetFitParameters().size(); s++){
			prop_Par[s] = SteppingKernel(curr_Par[s],fitSigma[s],fitUL[s],fitLL[s]);
		}

		double	*prop_arr = prop_Par.data();
		double	prop_Prob = -theFCN(prop_arr);
		double	A = exp(prop_Prob - curr_Prob);
		if(accepted.size() > 50)
			accepted.erase(accepted.begin());
		if(std::min(1.,A) > sampler.Rndm()){
			parTracker.push_back(prop_Par);
			curr_Par	= prop_Par;
			curr_Prob 	= prop_Prob;
			probTracker.push_back(prop_Prob);
			accepted.push_back(1);
		}
		else{
			parTracker.push_back(curr_Par);
			accepted.push_back(0);
			probTracker.push_back(curr_Prob);
		}
		if((c % 20) == 0){
			double sum = 0;
			for (auto& n : accepted)
				sum += n;
			sum /= (double)accepted.size();
			if(GetVerbose()){
				t1 = Clock::now();
				seconds s = std::chrono::duration_cast<seconds>(t1-t0);
				std::cout	<< std::setw(12) << std::left << c
						<< std::setw(15) << std::left << sum * 100
						<< std::setw(15) << std::left << s.count()
						<< "\r"
						<< std::flush; 
			}
		}
	}
	if(GetVerbose())
		std::cout << std::endl;	

}

double MCMCWalker::SteppingKernel(double x, double s, double UL, double LL){

	double val = x + rand.Gaus(0,s);
	while(val < LL || val > UL)
		val = x + rand.Gaus(0,s);	

	return val; 

}

void MCMCWalker::WriteTrackedParameters(const char *outfilename){

	std::ofstream outfile(outfilename);
	for(size_t s = 0; s<parTracker.size(); s++){
		for(size_t p = 0; p<parTracker.at(s).size();p++)
			outfile 	<< std::setw(12) << std::left << parTracker.at(s).at(p);
		outfile		<< std::setw(12) << std::left << probTracker.at(s); 
		outfile << "\n";
	}
	

}
