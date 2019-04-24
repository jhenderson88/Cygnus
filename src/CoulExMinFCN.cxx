#include "CoulExMinFCN.h"

// This is the function which Minuit2 will attempt to minimize

void ThreadTask(PointCoulEx &p, double theta){

	p.CalculatePointProbabilities(theta);	

	return;
	
}

void CoulExMinFCN::SetupCalculation(){

	exptIndex.resize(exptData.size());
	for(unsigned int i=0;i<scalingParameters.size();i++)
		for(unsigned int s=0;s<scalingParameters.at(i).GetExperimentNumbers().size();s++)
			exptIndex[scalingParameters.at(i).GetExperimentNumbers().at(s)] = (int)i;

	std::cout	<< std::setw(13) << std::left << "Experiment: "
			<< std::setw(14) << std::left << "Scaling index: "
			<< std::setw(14) << std::left << "Starting scaling:"
			<< std::endl;
	for(unsigned int i=0;i<exptIndex.size();i++){
		std::cout	<< std::setw(13) << std::left << i+1
				<< std::setw(14) << std::left << exptIndex.at(i)
				<< std::endl;
	}
		

}

double CoulExMinFCN::operator()(const double* par){

	double chisq = 0;
	int NDF = 0;

	typedef std::chrono::high_resolution_clock Clock;
	typedef std::chrono::milliseconds milliseconds;
	Clock::time_point t0 = Clock::now();

	Nucleus nucl = fNucleus;	

	for(unsigned int i=0;i<ME.size();i++){
		nucl.SetMatrixElement(ME.at(i).GetLambda(),ME.at(i).GetInitialState(),ME.at(i).GetFinalState(),par[i]);
		std::cout	<< par[i] 
				<< std::endl;
	}

	// 	COMPARE WITH LITERATURE CONSTRAINTS:
	TransitionRates rates(&nucl);
	double lifetime_chisq = 0;	
	for(unsigned int i=0;i<litLifetimes.size();i++){
		double 	tmp = 0;
		int	index		= litLifetimes.at(i).GetIndex();
		double	lifetime	= litLifetimes.at(i).GetLifetime();
		double	calcLifetime	= rates.GetLifetimes()[index];
		if(calcLifetime > lifetime)
			tmp = (calcLifetime - lifetime) / litLifetimes.at(i).GetUpUnc();
		else
			tmp = (calcLifetime - lifetime) / litLifetimes.at(i).GetDnUnc();
		chisq += tmp * tmp;
		lifetime_chisq += tmp*tmp;
		NDF++;
	}
	double br_chisq = 0;
	for(unsigned int i=0;i<litBranchingRatios.size();i++){
		double 	tmp 		= 0;
		int	index_init 	= litBranchingRatios.at(i).GetInitialIndex();
		int	index_final1	= litBranchingRatios.at(i).GetFinalIndex_1();
		int	index_final2	= litBranchingRatios.at(i).GetFinalIndex_2();
		double	BR		= litBranchingRatios.at(i).GetBranchingRatio();
		double  calcBR		= rates.GetBranchingRatios()[index_final1][index_init] / rates.GetBranchingRatios()[index_final2][index_init];
		if(calcBR > BR)
			tmp = (BR - calcBR) / litBranchingRatios.at(i).GetUpUnc();
		else
			tmp = (BR - calcBR) / litBranchingRatios.at(i).GetDnUnc();
		chisq += tmp * tmp;
		br_chisq += tmp*tmp;
		NDF++;
	}
	double mr_chisq = 0;
	for(unsigned int i=0;i<litMixingRatios.size();i++){
		double tmp;
		int 	index_init	= litMixingRatios.at(i).GetInitialIndex();
		int	index_final	= litMixingRatios.at(i).GetFinalIndex();
		double	delta		= litMixingRatios.at(i).GetMixingRatio();
		double	calcDelta	= rates.GetMixingRatios()[index_final][index_init];
		if(calcDelta > delta)
			tmp = (delta - calcDelta) / litMixingRatios.at(i).GetUpUnc();
		else
			tmp = (delta - calcDelta) / litMixingRatios.at(i).GetDnUnc();
		chisq += tmp * tmp;		
		mr_chisq += tmp*tmp;
		NDF++;
	}
	if(verbose)
		std::cout << "Literature chi-squared: " << chisq << std::endl;
	double litchisq = chisq;
	//	COULEX AND STUFF:
	
	//	First, new point calculations with new matrix elements
	for(unsigned int i=0;i<pointCalcs.size();i++){
		pointCalcs.at(i).SetNucleus(&nucl);
	}
	if(nThreads > 1){
		std::vector<std::thread> Threads;
		Threads.resize(nThreads-1);
		size_t exptcounter = 0;
		while(exptcounter < pointCalcs.size()){
			size_t activethreads = 0;
			for(size_t t=0; t < (size_t)(nThreads-1); t++){
				if(exptcounter < pointCalcs.size()){
					Threads[t] = std::thread(ThreadTask, std::ref(pointCalcs.at(exptcounter)), exptData.at(exptcounter).GetThetaCM());
					activethreads++;
				}		
				exptcounter++;
			}
			if(exptcounter < pointCalcs.size()){
				ThreadTask(std::ref(pointCalcs.at(exptcounter)), exptData.at(exptcounter).GetThetaCM());
				exptcounter++;
			}
			for(size_t t=0; t < activethreads; t++)
				Threads[t].join();
		}
	}
	else{
		for(unsigned int i=0;i<pointCalcs.size();i++)
			pointCalcs.at(i).CalculatePointProbabilities(exptData.at(i).GetThetaCM());
	}
	EffectiveCrossSection.clear();	// "Effective" cross section = cross section + feeding
	//	Now, determine the cross-sections to each state

	for(unsigned int i=0;i<pointCalcs.size();i++){
		TVectorD tmpVec;
		tmpVec.ResizeTo(pointCalcs.at(i).GetProbabilitiesVector().GetNrows());
		tmpVec	= pointCalcs.at(i).GetProbabilitiesVector();

		for(int s = 0; s<tmpVec.GetNrows(); s++){
			tmpVec[s] *= pointCalcs.at(i).GetReaction()->RutherfordCM(exptData.at(i).GetThetaCM()) * correctionFactors.at(i)[s] * 1e5;
		}
		TMatrixD tmpMat;
		tmpMat.ResizeTo(rates.GetBranchingRatios().GetNrows(),rates.GetBranchingRatios().GetNcols());
		tmpMat = GammaYield::GammaRayYield(tmpVec,rates.GetBranchingRatios());
		EffectiveCrossSection.push_back(tmpMat);
	}
	if(verbose)
		std::cout << std::endl;

	std::vector<double>	scaling;
	scaling.resize(exptData.size());
	for(unsigned int i=0;i<exptData.size();i++){
		scaling.at(i) = par[ME.size() + exptIndex.at(i)];
	}

	if(verbose)
		std::cout << "Experiment scaling: " << scaling.at(0) << std::endl;
	if(verbose){
		std::cout 	<< std::setw(7) << std::left << "Expt:";
		for(unsigned int t=0;t<exptData.at(0).GetData().size();t++){
			std::cout 	<< std::setw( 6) << std::left << "Init:"
					<< std::setw( 6) << std::left << "Finl:"
					<< std::setw(12) << std::left << "Calc:"
					<< std::setw(12) << std::left << "Expt:"
					<< std::setw(12) << std::left << "Calc/Expt:"
					<< std::setw(14) << std::left << "Chisq cont.:";
		}
		std::cout 	<< std::endl;
	}
	
	for(unsigned int i=0;i<exptData.size();i++){
		if(verbose)
			std::cout	<< std::setw(7) << std::left << i+1;
		for(unsigned int t=0;t<exptData.at(i).GetData().size();++t){
			double 	tmp 		= 0;
			int	index_init 	= exptData.at(i).GetData().at(t).GetInitialIndex();
			int	index_final 	= exptData.at(i).GetData().at(t).GetFinalIndex();
			double 	calcCounts 	= scaling.at(i) * EffectiveCrossSection.at(i)[index_final][index_init];
			double 	exptCounts 	= exptData.at(i).GetData().at(t).GetCounts();
			if(verbose){
				std::cout 	<< std::setw( 6) << std::left << index_init 
						<< std::setw( 6) << std::left << index_final 
						<< std::setw(12) << std::left << calcCounts 
						<< std::setw(12) << std::left << exptCounts 
						<< std::setw(12) << std::left << calcCounts/exptCounts
						<< std::setw(14) << std::left << TMath::Power((calcCounts - exptCounts)/exptData.at(i).GetData().at(t).GetUpUnc(),2);
			}
			if(calcCounts > exptCounts)
				tmp 		= (calcCounts - exptCounts) / exptData.at(i).GetData().at(t).GetUpUnc();
			else
				tmp 		= (calcCounts - exptCounts) / exptData.at(i).GetData().at(t).GetDnUnc();
			chisq		+= tmp * tmp;
			NDF++;
		}
		if(verbose)
			std::cout << std::endl;
	}
	if(verbose)
		std::cout << "Chisq: " << chisq << " scaling: " << par[nPar-1] << std::endl;

	iter++;

	Clock::time_point t1 = Clock::now();
	milliseconds ms = std::chrono::duration_cast<milliseconds>(t1-t0);

	if((iter % 10) == 0 && !verbose)
			std::cout 	<< "Iteration: " << std::setw(5) << iter 
					<< ", chi-squared value: " << std::setw(7) << chisq 
					<< " NDF: " << NDF 
					<< ", reduced chi-squared: " << std::setw(7) << chisq / (double)NDF
					<< ", literature component of chisq: " << std::setw(7) << litchisq 
					<< " processing time:  " << std::setw(5) << ms.count() << " ms\r" << std::flush;

	return chisq;	

}

void CoulExMinFCN::ClearAll(){

	parameters.clear();			
	ME.clear();			
	pointCalcs.clear();			
	exptData.clear();			
	litLifetimes.clear();			
	litBranchingRatios.clear();		
	litMixingRatios.clear();		
	EffectiveCrossSection.clear();

}
