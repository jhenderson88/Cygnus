#include "CoulExSimMinFCN.h"

// This is the function which Minuit2 will attempt to minimize
void SimMinThreadTask(PointCoulEx &p, double theta){

	p.CalculatePointProbabilities(theta);	

	return;
	
}



void CoulExSimMinFCN::SetupCalculation(){

	exptIndex.resize(exptData_Beam.size());
	for(unsigned int i=0;i<scalingParameters.size();i++)
		for(unsigned int s=0;s<scalingParameters.at(i).GetExperimentNumbers().size();s++)
			exptIndex[scalingParameters.at(i).GetExperimentNumbers().at(s)] = (int)i;

	std::cout	<< std::setw(13) << std::left << "Experiment: "
			<< std::setw(14) << std::left << "Scaling index: "
			<< std::endl;
	for(unsigned int i=0;i<exptIndex.size();i++){
		std::cout	<< std::setw(13) << std::left << i+1
				<< std::setw(14) << std::left << exptIndex.at(i)
				<< std::endl;
	}
		

}

double CoulExSimMinFCN::operator()(const double* par){

	double chisq = 0;
	int NDF = 0;

	typedef std::chrono::high_resolution_clock Clock;
	typedef std::chrono::milliseconds milliseconds;
	Clock::time_point t0 = Clock::now();

	Nucleus nucl_b = fNucleus_Beam;	
	Nucleus nucl_t = fNucleus_Target;

	for(unsigned int i=0;i<ME_Beam.size();i++)
		nucl_b.SetMatrixElement(ME_Beam.at(i).GetLambda(),ME_Beam.at(i).GetInitialState(),ME_Beam.at(i).GetFinalState(),par[i]);
	for(unsigned int i=0;i<ME_Target.size();i++)
		nucl_t.SetMatrixElement(ME_Target.at(i).GetLambda(),ME_Target.at(i).GetInitialState(),ME_Target.at(i).GetFinalState(),par[i + ME_Beam.size()]);

	// 	COMPARE WITH LITERATURE CONSTRAINTS:
	// 	First, compare with the literature for the beam:
	TransitionRates rates_b(&nucl_b);
	double lifetime_chisq = 0;	
	for(unsigned int i=0;i<litLifetimes_Beam.size();i++){
		double 	tmp = 0;
		int	index		= litLifetimes_Beam.at(i).GetIndex();
		double	lifetime	= litLifetimes_Beam.at(i).GetLifetime();
		double	calcLifetime	= rates_b.GetLifetimes()[index];
		if(calcLifetime > lifetime)
			tmp = (calcLifetime - lifetime) / litLifetimes_Beam.at(i).GetUpUnc();
		else
			tmp = (calcLifetime - lifetime) / litLifetimes_Beam.at(i).GetDnUnc();
		chisq += tmp * tmp;
		lifetime_chisq += tmp*tmp;
		NDF++;
	}
	double br_chisq = 0;
	for(unsigned int i=0;i<litBranchingRatios_Beam.size();i++){
		double 	tmp 		= 0;
		int	index_init 	= litBranchingRatios_Beam.at(i).GetInitialIndex();
		int	index_final1	= litBranchingRatios_Beam.at(i).GetFinalIndex_1();
		int	index_final2	= litBranchingRatios_Beam.at(i).GetFinalIndex_2();
		double	BR		= litBranchingRatios_Beam.at(i).GetBranchingRatio();
		double  calcBR		= rates_b.GetBranchingRatios()[index_final1][index_init] / rates_b.GetBranchingRatios()[index_final2][index_init];
		if(calcBR > BR)
			tmp = (BR - calcBR) / litBranchingRatios_Beam.at(i).GetUpUnc();
		else
			tmp = (BR - calcBR) / litBranchingRatios_Beam.at(i).GetDnUnc();
		chisq += tmp * tmp;
		br_chisq += tmp*tmp;
		NDF++;
	}
	double mr_chisq = 0;
	for(unsigned int i=0;i<litMixingRatios_Beam.size();i++){
		double tmp;
		int 	index_init	= litMixingRatios_Beam.at(i).GetInitialIndex();
		int	index_final	= litMixingRatios_Beam.at(i).GetFinalIndex();
		double	delta		= litMixingRatios_Beam.at(i).GetMixingRatio();
		double	calcDelta	= rates_b.GetMixingRatios()[index_final][index_init];
		if(calcDelta > delta)
			tmp = (delta - calcDelta) / litMixingRatios_Beam.at(i).GetUpUnc();
		else
			tmp = (delta - calcDelta) / litMixingRatios_Beam.at(i).GetDnUnc();
		chisq += tmp * tmp;		
		mr_chisq += tmp*tmp;
		NDF++;
	}
	// 	Now, compare with the literature for the target:
	TransitionRates rates_t(&nucl_b);
	for(unsigned int i=0;i<litLifetimes_Target.size();i++){
		double 	tmp = 0;
		int	index		= litLifetimes_Target.at(i).GetIndex();
		double	lifetime	= litLifetimes_Target.at(i).GetLifetime();
		double	calcLifetime	= rates_t.GetLifetimes()[index];
		if(calcLifetime > lifetime)
			tmp = (calcLifetime - lifetime) / litLifetimes_Target.at(i).GetUpUnc();
		else
			tmp = (calcLifetime - lifetime) / litLifetimes_Target.at(i).GetDnUnc();
		chisq += tmp * tmp;
		lifetime_chisq += tmp*tmp;
		NDF++;
	}
	for(unsigned int i=0;i<litBranchingRatios_Target.size();i++){
		double 	tmp 		= 0;
		int	index_init 	= litBranchingRatios_Target.at(i).GetInitialIndex();
		int	index_final1	= litBranchingRatios_Target.at(i).GetFinalIndex_1();
		int	index_final2	= litBranchingRatios_Target.at(i).GetFinalIndex_2();
		double	BR		= litBranchingRatios_Target.at(i).GetBranchingRatio();
		double  calcBR		= rates_t.GetBranchingRatios()[index_final1][index_init] / rates_t.GetBranchingRatios()[index_final2][index_init];
		if(calcBR > BR)
			tmp = (BR - calcBR) / litBranchingRatios_Target.at(i).GetUpUnc();
		else
			tmp = (BR - calcBR) / litBranchingRatios_Target.at(i).GetDnUnc();
		chisq += tmp * tmp;
		br_chisq += tmp*tmp;
		NDF++;
	}
	for(unsigned int i=0;i<litMixingRatios_Target.size();i++){
		double tmp;
		int 	index_init	= litMixingRatios_Target.at(i).GetInitialIndex();
		int	index_final	= litMixingRatios_Target.at(i).GetFinalIndex();
		double	delta		= litMixingRatios_Target.at(i).GetMixingRatio();
		double	calcDelta	= rates_t.GetMixingRatios()[index_final][index_init];
		if(calcDelta > delta)
			tmp = (delta - calcDelta) / litMixingRatios_Target.at(i).GetUpUnc();
		else
			tmp = (delta - calcDelta) / litMixingRatios_Target.at(i).GetDnUnc();
		chisq += tmp * tmp;		
		mr_chisq += tmp*tmp;
		NDF++;
	}
	if(verbose)
		std::cout << "Literature chi-squared: " << chisq << std::endl;
	double litchisq = chisq;
	//	COULEX AND STUFF:
	
	//	First, new point calculations with new matrix elements for the beam:
	for(unsigned int i=0;i<pointCalcs_Beam.size();i++){
		pointCalcs_Beam.at(i).SetNucleus(&nucl_b);
	}
	if(nThreads > 1){
		std::vector<std::thread> Threads;
		Threads.resize(nThreads-1);
		size_t exptcounter = 0;
		while(exptcounter < pointCalcs_Beam.size()){
			size_t activethreads = 0;
			for(size_t t=0; t < (size_t)(nThreads-1); t++){
				if(exptcounter < pointCalcs_Beam.size()){
					Threads[t] = std::thread(SimMinThreadTask, std::ref(pointCalcs_Beam.at(exptcounter)), exptData_Beam.at(exptcounter).GetThetaCM());
					activethreads++;
				}		
				exptcounter++;
			}
			if(exptcounter < pointCalcs_Beam.size()){
				SimMinThreadTask(std::ref(pointCalcs_Beam.at(exptcounter)), exptData_Beam.at(exptcounter).GetThetaCM());
				exptcounter++;
			}
			for(size_t t=0; t < activethreads; t++)
				Threads[t].join();
		}
	}
	else{
		for(unsigned int i=0;i<pointCalcs_Beam.size();i++)
			pointCalcs_Beam.at(i).CalculatePointProbabilities(exptData_Beam.at(i).GetThetaCM());
	}
	//	Now, new point calculations with new matrix elements for the target:
	for(unsigned int i=0;i<pointCalcs_Target.size();i++){
		pointCalcs_Target.at(i).SetNucleus(&nucl_t);
	}
	if(nThreads > 1){
		std::vector<std::thread> Threads;
		Threads.resize(nThreads-1);
		size_t exptcounter = 0;
		while(exptcounter < pointCalcs_Target.size()){
			size_t activethreads = 0;
			for(size_t t=0; t < (size_t)(nThreads-1); t++){
				if(exptcounter < pointCalcs_Target.size()){
					Threads[t] = std::thread(SimMinThreadTask, std::ref(pointCalcs_Target.at(exptcounter)), exptData_Target.at(exptcounter).GetThetaCM());
					activethreads++;
				}		
				exptcounter++;
			}
			if(exptcounter < pointCalcs_Target.size()){
				SimMinThreadTask(std::ref(pointCalcs_Target.at(exptcounter)), exptData_Target.at(exptcounter).GetThetaCM());
				exptcounter++;
			}
			for(size_t t=0; t < activethreads; t++)
				Threads[t].join();
		}
	}
	else{
		for(unsigned int i=0;i<pointCalcs_Target.size();i++)
			pointCalcs_Target.at(i).CalculatePointProbabilities(exptData_Target.at(i).GetThetaCM());
	}

	//	Prepare the calculations for comparison, first beam:
	EffectiveCrossSection_Beam.clear();	
	for(unsigned int i=0;i<pointCalcs_Beam.size();i++){
		TVectorD tmpVec;
		tmpVec.ResizeTo(pointCalcs_Beam.at(i).GetProbabilitiesVector().GetNrows());
		tmpVec	= pointCalcs_Beam.at(i).GetProbabilitiesVector();
		for(int s = 0; s<tmpVec.GetNrows(); s++){
			tmpVec[s] *= pointCalcs_Beam.at(i).GetReaction()->RutherfordCM(exptData_Beam.at(i).GetThetaCM()) * correctionFactors_Beam.at(i)[s] * 1e5;
		}
		TMatrixD tmpMat;
		tmpMat.ResizeTo(rates_b.GetBranchingRatios().GetNrows(),rates_b.GetBranchingRatios().GetNcols());
		tmpMat = GammaYield::GammaRayYield(tmpVec,rates_b.GetBranchingRatios());
		EffectiveCrossSection_Beam.push_back(tmpMat);
	}
	//	Now prepare the calculations for comparison with the target calculation:
	EffectiveCrossSection_Target.clear();	
	for(unsigned int i=0;i<pointCalcs_Target.size();i++){
		TVectorD tmpVec;
		tmpVec.ResizeTo(pointCalcs_Target.at(i).GetProbabilitiesVector().GetNrows());
		tmpVec	= pointCalcs_Target.at(i).GetProbabilitiesVector();
		for(int s = 0; s<tmpVec.GetNrows(); s++){
			tmpVec[s] *= pointCalcs_Target.at(i).GetReaction()->RutherfordCM(exptData_Target.at(i).GetThetaCM()) * correctionFactors_Target.at(i)[s] * 1e5;
		}
		TMatrixD tmpMat;
		tmpMat.ResizeTo(rates_t.GetBranchingRatios().GetNrows(),rates_t.GetBranchingRatios().GetNcols());
		tmpMat = GammaYield::GammaRayYield(tmpVec,rates_t.GetBranchingRatios());
		EffectiveCrossSection_Target.push_back(tmpMat);
	}
	if(verbose)
		std::cout << std::endl;

	//	Scaling is common for both beam and target, so no need to duplicate this:
	std::vector<double>	scaling;
	scaling.resize(exptData_Beam.size()); 
	// Number of experiments better be the same for target and beam, or we're in trouble
	for(unsigned int i=0;i<exptData_Beam.size();i++)
		scaling.at(i) = par[ME_Beam.size() + ME_Target.size() + exptIndex.at(i)];

	if(verbose)
		std::cout << "Experiment scaling: " << scaling.at(0) << std::endl;
	//	Everything needs printing for both beam and target...
	if(verbose){
		std::cout 	<< std::setw(7) << std::left << "Expt:";
		for(unsigned int t=0;t<exptData_Beam.at(0).GetData().size();t++){
			std::cout 	<< std::setw( 6) << std::left << "Init:"
					<< std::setw( 6) << std::left << "Finl:"
					<< std::setw(10) << std::left << "Calc:"
					<< std::setw(10) << std::left << "Expt:"
					<< std::setw(12) << std::left << "Calc/Expt:"
					<< std::setw(14) << std::left << "Chisq cont.:";
		}
		std::cout 	<< std::endl;
	}	
	for(unsigned int i=0;i<exptData_Beam.size();i++){
		if(verbose)
			std::cout	<< std::setw(7) << std::left << i+1;
		for(unsigned int t=0;t<exptData_Beam.at(i).GetData().size();++t){
			double 	tmp 		= 0;
			int	index_init 	= exptData_Beam.at(i).GetData().at(t).GetInitialIndex();
			int	index_final 	= exptData_Beam.at(i).GetData().at(t).GetFinalIndex();
			double 	calcCounts 	= scaling.at(i) * EffectiveCrossSection_Beam.at(i)[index_final][index_init];
			double 	exptCounts 	= exptData_Beam.at(i).GetData().at(t).GetCounts();
			if(verbose){
				std::cout 	<< std::setw( 6) << std::left << index_init 
						<< std::setw( 6) << std::left << index_final 
						<< std::setw(10) << std::left << calcCounts 
						<< std::setw(10) << std::left << exptCounts 
						<< std::setw(12) << std::left << calcCounts/exptCounts
						<< std::setw(14) << std::left << TMath::Power((calcCounts - exptCounts)/exptData_Beam.at(i).GetData().at(t).GetUpUnc(),2);
			}
			if(calcCounts > exptCounts)
				tmp 		= (calcCounts - exptCounts) / exptData_Beam.at(i).GetData().at(t).GetUpUnc();
			else
				tmp 		= (calcCounts - exptCounts) / exptData_Beam.at(i).GetData().at(t).GetDnUnc();
			chisq 			+= tmp * tmp;
			NDF++;
		}
		if(verbose)
			std::cout << std::endl;
	}
	if(verbose){
		std::cout 	<< std::setw(7) << std::left << "Expt:";
		for(unsigned int t=0;t<exptData_Target.at(0).GetData().size();t++){
			std::cout 	<< std::setw( 6) << std::left << "Init:"
					<< std::setw( 6) << std::left << "Finl:"
					<< std::setw(10) << std::left << "Calc:"
					<< std::setw(10) << std::left << "Expt:"
					<< std::setw(12) << std::left << "Calc/Expt:"
					<< std::setw(14) << std::left << "Chisq cont.:";
		}
		std::cout 	<< std::endl;
	}	
	for(unsigned int i=0;i<exptData_Target.size();i++){
		if(verbose)
			std::cout	<< std::setw(7) << std::left << i+1;
		for(unsigned int t=0;t<exptData_Target.at(i).GetData().size();++t){
			double 	tmp 		= 0;
			int	index_init 	= exptData_Target.at(i).GetData().at(t).GetInitialIndex();
			int	index_final 	= exptData_Target.at(i).GetData().at(t).GetFinalIndex();
			double 	calcCounts 	= scaling.at(i) * EffectiveCrossSection_Target.at(i)[index_final][index_init];
			double 	exptCounts 	= exptData_Target.at(i).GetData().at(t).GetCounts();
			if(verbose){
				std::cout 	<< std::setw( 6) << std::left << index_init 
						<< std::setw( 6) << std::left << index_final 
						<< std::setw(10) << std::left << calcCounts 
						<< std::setw(10) << std::left << exptCounts 
						<< std::setw(12) << std::left << calcCounts/exptCounts
						<< std::setw(14) << std::left << TMath::Power((calcCounts - exptCounts)/exptData_Target.at(i).GetData().at(t).GetUpUnc(),2);
			}
			if(calcCounts > exptCounts)
				tmp 		= (calcCounts - exptCounts) / exptData_Target.at(i).GetData().at(t).GetUpUnc();
			else
				tmp 		= (calcCounts - exptCounts) / exptData_Target.at(i).GetData().at(t).GetDnUnc();
			chisq 			+= tmp * tmp;
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
			std::cout 	<< std::setw(13) << std::left << "Iteration: " 
					<< std::setw(7)  << std::left << iter 
					<< std::setw(24) << std::left << " chi-squared value: " 
					<< std::setw(7)  << std::left << chisq 
					<< std::setw(7)  << std::left << " NDF: " 
					<< std::setw(4)  << std::left << NDF 
					<< std::setw(26) << std::left << " reduced chi-squared: " 
					<< std::setw(7)  << std::left << chisq / (double)NDF
					<< std::setw(36) << std::left << " literature component of chisq: " 
					<< std::setw(7)  << std::left << litchisq 
					<< std::setw(24) << std::left << " processing time:  " 
					<< std::setw(7)  << std::left << ms.count() 
					<< " ms\r" << std::flush;

	return chisq;	

}

void CoulExSimMinFCN::ClearAll(){

	parameters.clear();			
	
	ME_Beam.clear();			
	pointCalcs_Beam.clear();			
	exptData_Beam.clear();			
	litLifetimes_Beam.clear();			
	litBranchingRatios_Beam.clear();		
	litMixingRatios_Beam.clear();		
	EffectiveCrossSection_Beam.clear();

	ME_Target.clear();			
	pointCalcs_Target.clear();			
	exptData_Target.clear();			
	litLifetimes_Target.clear();			
	litBranchingRatios_Target.clear();		
	litMixingRatios_Target.clear();		
	EffectiveCrossSection_Target.clear();

}
