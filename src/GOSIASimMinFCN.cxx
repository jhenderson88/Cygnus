#include "GOSIASimMinFCN.h"

void GOSIASimMinFCN::SetupCalculation(){

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

double GOSIASimMinFCN::operator()(const double* par){

	double chisq = 0;
	int NDF = 0;

	typedef std::chrono::high_resolution_clock Clock;
	typedef std::chrono::milliseconds milliseconds;
	Clock::time_point t0 = Clock::now();


	Nucleus nucl_b = fNucleus_Beam;	
	Nucleus nucl_t = fNucleus_Target;

	parameters.clear();
	for(unsigned int i=0;i<ME_Beam.size();i++){
		nucl_b.SetMatrixElement(ME_Beam.at(i).GetLambda(),ME_Beam.at(i).GetInitialState(),ME_Beam.at(i).GetFinalState(),par[i]);
		parameters.push_back(par[i]);
	}
	for(unsigned int i=0;i<ME_Target.size();i++){
		nucl_t.SetMatrixElement(ME_Target.at(i).GetLambda(),ME_Target.at(i).GetInitialState(),ME_Target.at(i).GetFinalState(),par[i + ME_Beam.size()]);
		parameters.push_back(par[i + ME_Beam.size()]);
	}

	if(verbose){
		std::cout	<< std::endl;
		std::cout	<< "Parameters:"
				<< std::endl;
		for(unsigned int i=0;i<ME_Beam.size();i++)
			std::cout	<< std::setw(14) << std::left << par[i];
		std::cout	<< std::endl;
		for(unsigned int i=0;i<ME_Target.size();i++)
			std::cout	<< std::setw(14) << std::left << par[i + ME_Beam.size()];
		std::cout	<< std::endl;
	}

	// 	COMPARE WITH LITERATURE CONSTRAINTS:
	// 	First, compare with the literature for the beam:
	TransitionRates rates_b(&nucl_b);
	double lifetime_chisq = 0;	
	for(unsigned int i=0;i<litLifetimes_Beam.size();i++){
		double 	tmp = 0;
		int	index		= litLifetimes_Beam.at(i).GetIndex();
		double	lifetime	= litLifetimes_Beam.at(i).GetLifetime();
		double	calcLifetime	= rates_b.GetLifetimes()[index];
		if(fLikelihood){
			double	sigma		= litLifetimes_Beam.at(i).GetUpUnc() * litLifetimes_Beam.at(i).GetDnUnc();
			double	sigma_prime	= (litLifetimes_Beam.at(i).GetUpUnc() - litLifetimes_Beam.at(i).GetDnUnc());
			chisq 			+= 0.5 * TMath::Power((calcLifetime - lifetime),2)/(sigma + sigma_prime * (calcLifetime - lifetime));
		}
		else{
			if(calcLifetime > lifetime)
				tmp = (calcLifetime - lifetime) / litLifetimes_Beam.at(i).GetUpUnc();
			else
				tmp = (calcLifetime - lifetime) / litLifetimes_Beam.at(i).GetDnUnc();
			chisq += tmp * tmp;
			lifetime_chisq += tmp*tmp;
		}
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
		if(fLikelihood){
			double	sigma		= litBranchingRatios_Beam.at(i).GetUpUnc() * litBranchingRatios_Beam.at(i).GetDnUnc();
			double	sigma_prime	= (litBranchingRatios_Beam.at(i).GetUpUnc() - litBranchingRatios_Beam.at(i).GetDnUnc());
			chisq 			+= 0.5 * TMath::Power((calcBR - BR),2)/(sigma + sigma_prime * (calcBR - BR));
		}
		else{
			if(calcBR > BR)
				tmp = (BR - calcBR) / litBranchingRatios_Beam.at(i).GetUpUnc();
			else
				tmp = (BR - calcBR) / litBranchingRatios_Beam.at(i).GetDnUnc();
			chisq += tmp * tmp;
			br_chisq += tmp*tmp;
		}
		NDF++;
	}
	double mr_chisq = 0;
	for(unsigned int i=0;i<litMixingRatios_Beam.size();i++){
		double tmp;
		int 	index_init	= litMixingRatios_Beam.at(i).GetInitialIndex();
		int	index_final	= litMixingRatios_Beam.at(i).GetFinalIndex();
		double	delta		= litMixingRatios_Beam.at(i).GetMixingRatio();
		double	calcDelta	= rates_b.GetMixingRatios()[index_final][index_init];
		if(fLikelihood){
			double	sigma		= litMixingRatios_Beam.at(i).GetUpUnc() * litMixingRatios_Beam.at(i).GetDnUnc();
			double	sigma_prime	= (litMixingRatios_Beam.at(i).GetUpUnc() - litMixingRatios_Beam.at(i).GetDnUnc());
			chisq 			+= 0.5 * TMath::Power((calcDelta - delta),2)/(sigma + sigma_prime * (calcDelta - delta));
		}
		else{
			if(calcDelta > delta)
				tmp = (delta - calcDelta) / litMixingRatios_Beam.at(i).GetUpUnc();
			else
				tmp = (delta - calcDelta) / litMixingRatios_Beam.at(i).GetDnUnc();
			chisq += tmp * tmp;		
			mr_chisq += tmp*tmp;
		}
		NDF++;
	}
	double me_chisq = 0;
	for(unsigned int i=0;i<litMatrixElements_Beam.size();i++){
		double tmp;
		int	mult		= litMatrixElements_Beam.at(i).GetMultipolarity();
		int 	index_init	= litMatrixElements_Beam.at(i).GetInitialIndex();
		int	index_final	= litMatrixElements_Beam.at(i).GetFinalIndex();
		double	ME		= litMatrixElements_Beam.at(i).GetMatrixElement();
		double	calcME		= nucl_b.GetMatrixElements().at(mult)[index_init][index_final];
		if(fLikelihood){
			double	sigma		= litMatrixElements_Beam.at(i).GetUpUnc() * litMatrixElements_Beam.at(i).GetDnUnc();
			double	sigma_prime	= (litMatrixElements_Beam.at(i).GetUpUnc() - litMatrixElements_Beam.at(i).GetDnUnc());
			chisq 			+= 0.5 * TMath::Power((calcME - ME),2)/(sigma + sigma_prime * (calcME - ME));
		}
		else{
			if(calcME > ME)
				tmp = (ME - calcME) / litMatrixElements_Beam.at(i).GetUpUnc();
			else
				tmp = (ME - calcME) / litMatrixElements_Beam.at(i).GetDnUnc();
			chisq += tmp * tmp;		
			me_chisq += tmp*tmp;
		}
		NDF++;
	}
	// 	Now, compare with the literature for the target:
	TransitionRates rates_t(&nucl_t);
	for(unsigned int i=0;i<litLifetimes_Target.size();i++){
		double 	tmp = 0;
		int	index		= litLifetimes_Target.at(i).GetIndex();
		double	lifetime	= litLifetimes_Target.at(i).GetLifetime();
		double	calcLifetime	= rates_t.GetLifetimes()[index];
		if(fLikelihood){
			double	sigma		= litLifetimes_Target.at(i).GetUpUnc() * litLifetimes_Target.at(i).GetDnUnc();
			double	sigma_prime	= (litLifetimes_Target.at(i).GetUpUnc() - litLifetimes_Target.at(i).GetDnUnc());
			chisq 			+= 0.5 * TMath::Power((calcLifetime - lifetime),2)/(sigma + sigma_prime * (calcLifetime - lifetime));
		}
		else{
			if(calcLifetime > lifetime)
				tmp = (calcLifetime - lifetime) / litLifetimes_Target.at(i).GetUpUnc();
			else
				tmp = (calcLifetime - lifetime) / litLifetimes_Target.at(i).GetDnUnc();
			chisq += tmp * tmp;
			lifetime_chisq += tmp*tmp;
		}
		NDF++;
	}
	for(unsigned int i=0;i<litBranchingRatios_Target.size();i++){
		double 	tmp 		= 0;
		int	index_init 	= litBranchingRatios_Target.at(i).GetInitialIndex();
		int	index_final1	= litBranchingRatios_Target.at(i).GetFinalIndex_1();
		int	index_final2	= litBranchingRatios_Target.at(i).GetFinalIndex_2();
		double	BR		= litBranchingRatios_Target.at(i).GetBranchingRatio();
		double  calcBR		= rates_t.GetBranchingRatios()[index_final1][index_init] / rates_t.GetBranchingRatios()[index_final2][index_init];
		if(fLikelihood){
			double	sigma		= litBranchingRatios_Target.at(i).GetUpUnc() * litBranchingRatios_Target.at(i).GetDnUnc();
			double	sigma_prime	= (litBranchingRatios_Target.at(i).GetUpUnc() - litBranchingRatios_Target.at(i).GetDnUnc());
			chisq 			+= 0.5 * TMath::Power((calcBR - BR),2)/(sigma + sigma_prime * (calcBR - BR));
		}
		else{
			if(calcBR > BR)
				tmp = (BR - calcBR) / litBranchingRatios_Target.at(i).GetUpUnc();
			else
				tmp = (BR - calcBR) / litBranchingRatios_Target.at(i).GetDnUnc();
			chisq += tmp * tmp;
			br_chisq += tmp*tmp;
		}
		NDF++;
	}
	for(unsigned int i=0;i<litMixingRatios_Target.size();i++){
		double tmp;
		int 	index_init	= litMixingRatios_Target.at(i).GetInitialIndex();
		int	index_final	= litMixingRatios_Target.at(i).GetFinalIndex();
		double	delta		= litMixingRatios_Target.at(i).GetMixingRatio();
		double	calcDelta	= rates_t.GetMixingRatios()[index_final][index_init];
		if(fLikelihood){
			double	sigma		= litMixingRatios_Target.at(i).GetUpUnc() * litMixingRatios_Target.at(i).GetDnUnc();
			double	sigma_prime	= (litMixingRatios_Target.at(i).GetUpUnc() - litMixingRatios_Target.at(i).GetDnUnc());
			chisq 			+= 0.5 * TMath::Power((calcDelta - delta),2)/(sigma + sigma_prime * (calcDelta - delta));
		}
		else{
			if(calcDelta > delta)
				tmp = (delta - calcDelta) / litMixingRatios_Target.at(i).GetUpUnc();
			else
				tmp = (delta - calcDelta) / litMixingRatios_Target.at(i).GetDnUnc();
			chisq += tmp * tmp;		
			mr_chisq += tmp*tmp;
		}
		NDF++;
	}
	for(unsigned int i=0;i<litMatrixElements_Target.size();i++){
		double tmp;
		int	mult		= litMatrixElements_Target.at(i).GetMultipolarity();
		int 	index_init	= litMatrixElements_Target.at(i).GetInitialIndex();
		int	index_final	= litMatrixElements_Target.at(i).GetFinalIndex();
		double	ME		= litMatrixElements_Target.at(i).GetMatrixElement();
		double	calcME		= nucl_t.GetMatrixElements().at(mult)[index_init][index_final];
		if(fLikelihood){
			double	sigma		= litMatrixElements_Target.at(i).GetUpUnc() * litMatrixElements_Target.at(i).GetDnUnc();
			double	sigma_prime	= (litMatrixElements_Target.at(i).GetUpUnc() - litMatrixElements_Target.at(i).GetDnUnc());
			chisq 			+= 0.5 * TMath::Power((calcME - ME),2)/(sigma + sigma_prime * (calcME - ME));
		}
		else{
			if(calcME > ME)
				tmp = (ME - calcME) / litMatrixElements_Target.at(i).GetUpUnc();
			else
				tmp = (ME - calcME) / litMatrixElements_Target.at(i).GetDnUnc();
			chisq += tmp * tmp;		
			me_chisq += tmp*tmp;
		}
		NDF++;
	}
	double litchisq = chisq;
	//	COULEX AND STUFF:
	
	std::ofstream		beam_bst(beamBSTFile);
	for(size_t i=0;i<beamMapping_i.size();i++){
		beam_bst << nucl_b.GetMatrixElements().at(beamMapping_l.at(i))[beamMapping_f.at(i)][beamMapping_i.at(i)] << "\n";
	}
	beam_bst.close();
	
	std::ofstream		target_bst(targetBSTFile);
	for(size_t i=0;i<targetMapping_i.size();i++){
		target_bst 	<< nucl_t.GetMatrixElements().at(targetMapping_l.at(i))[targetMapping_f.at(i)][targetMapping_i.at(i)] << "\n";
	}
	target_bst.close();

	std::string	str;
	str = "./gosia < "+ beamGOSIAFile_inp+"> /dev/null";
	const char* 	c_b = str.c_str();
	system(c_b);		// 	Run the beam GOSIA file
	str = "./gosia < "+ targetGOSIAFile_inp+"> /dev/null";
	const char* 	c_t = str.c_str();
	system(c_t);	// 	Run the target GOSIA file

	const char	*b_out = beamGOSIAFile_out.c_str();
	const char	*t_out = targetGOSIAFile_out.c_str();
	GOSIAReader	beam_gosiaReader(&nucl_b,b_out);	//	Grab the GOSIA yields
	GOSIAReader	target_gosiaReader(&nucl_t,t_out);	//	Grab the GOSIA yields

	std::vector<ExperimentData>	beamCalc	= beam_gosiaReader.GetGOSIAData();
	std::vector<ExperimentData>	targetCalc	= target_gosiaReader.GetGOSIAData();
	EffectiveCrossSection_Beam.clear();	
	EffectiveCrossSection_Target.clear();	

	for(size_t i=0; i<beamCalc.size(); i++){
		TMatrixD	tmpMat;
		tmpMat.ResizeTo(rates_b.GetBranchingRatios().GetNrows(),rates_b.GetBranchingRatios().GetNcols());
		size_t	nRows = beamCalc.at(i).GetData().size();
		for(size_t j=0; j<nRows; j++){
			int	init		= beamCalc.at(i).GetData().at(j).GetInitialIndex();
			int	fina		= beamCalc.at(i).GetData().at(j).GetFinalIndex();
			double 	counts 		= beamCalc.at(i).GetData().at(j).GetCounts();
			if(false && i==0){
				std::cout	<< std::setw( 6) << std::left << init
						<< std::setw( 6) << std::left << fina
						<< std::setw(12) << std::left << counts	
						<< std::setw(12) << std::left << correctionFactors_Beam.at(i)[j]
						<< std::setw(12) << std::left << correctionFactors_Beam.at(i)[j] * counts
						<< std::endl;
			}
			tmpMat[fina][init]	= counts * correctionFactors_Beam.at(i)[j];
			tmpMat[init][fina]	= counts * correctionFactors_Beam.at(i)[j];
		}
		EffectiveCrossSection_Beam.push_back(tmpMat);
	}

	for(size_t i=0; i<targetCalc.size(); i++){
		TMatrixD	tmpMat;
		tmpMat.ResizeTo(rates_t.GetBranchingRatios().GetNrows(),rates_t.GetBranchingRatios().GetNcols());
		size_t	nRows = targetCalc.at(i).GetData().size();
		for(size_t j=0; j<nRows; j++){
			int	init		= targetCalc.at(i).GetData().at(j).GetInitialIndex();
			int	fina		= targetCalc.at(i).GetData().at(j).GetFinalIndex();
			double 	counts 		= targetCalc.at(i).GetData().at(j).GetCounts();
			if(false && i==0){
				std::cout	<< std::setw( 6) << std::left << init
						<< std::setw( 6) << std::left << fina
						<< std::setw(12) << std::left << counts	
						<< std::setw(12) << std::left << correctionFactors_Target.at(i)[j]
						<< std::setw(12) << std::left << correctionFactors_Target.at(i)[j] * counts
						<< std::endl;
			}
			tmpMat[fina][init]	= counts * correctionFactors_Target.at(i)[j];
			tmpMat[init][fina]	= counts * correctionFactors_Target.at(i)[j];
		}
		EffectiveCrossSection_Target.push_back(tmpMat);
	}

	if(verbose)
		std::cout << std::endl;

	//	Everything needs printing for both beam and target...
	if(verbose && !fLikelihood){
		std::cout 	<< std::setw(7) << std::left << "Beam:"
				<< std::endl;
		std::cout 	<< std::setw( 4) << std::left << " ";
		for(unsigned int t=0;t<exptData_Beam.at(exptData_Beam.size()-1).GetData().size();t++){
			std::cout 	<< std::setw( 6) << std::left << "Init:"
					<< std::setw( 6) << std::left << "Finl:"
					<< std::setw(10) << std::left << "Calc:"
					<< std::setw(10) << std::left << "Expt:"
					<< std::setw(10) << std::left << "C/E:"
					<< std::setw(12) << std::left << "Chisq:";
		}
		std::cout 	<< std::endl;
	}
	else if(verbose){
		std::cout 	<< std::setw(7) << std::left << "Beam:"
				<< std::endl;
		std::cout 	<< std::setw( 4) << std::left << " ";
		for(unsigned int t=0;t<exptData_Beam.at(exptData_Beam.size()-1).GetData().size();t++){
			std::cout	<< std::setw( 6) << std::left << "Init:"
					<< std::setw( 6) << std::left << "Finl:"
					<< std::setw(10) << std::left << "Calc:"
					<< std::setw(10) << std::left << "Expt:"
					<< std::setw(10) << std::left << "Calc/Expt:"
					<< std::setw(12) << std::left << "-Ln(L) cont.:";
		}
		std::cout 	<< std::endl;
	}

	std::vector<double>	scaling;
	scaling.resize(exptData_Beam.size());
	for(unsigned int s=0;s<scalingParameters.size();s++){
		double	calScaling 	= 0;
		double	weightSum	= 0;
		for(unsigned int ss=0;ss<scalingParameters.at(s).GetExperimentNumbers().size();ss++){
			unsigned int i = scalingParameters.at(s).GetExperimentNumbers().at(ss);
			for(unsigned int t=0;t<exptData_Beam.at(i).GetData().size();++t){
				int	index_init 	= exptData_Beam.at(i).GetData().at(t).GetInitialIndex();
				int	index_final 	= exptData_Beam.at(i).GetData().at(t).GetFinalIndex();
				double 	calcCounts 	= EffectiveCrossSection_Beam.at(i)[index_final][index_init];
				double 	exptCounts 	= exptData_Beam.at(i).GetData().at(t).GetCounts();
				if(exptCounts > 0 && calcCounts > 0){
					double	sigma		= (exptData_Beam.at(i).GetData().at(t).GetUpUnc() + exptData_Beam.at(i).GetData().at(t).GetDnUnc())/2.;  // Average uncertainty
					double	ratio		= exptCounts / calcCounts;
					double	r_sigma		= ratio * sigma / exptCounts;
					double	weight		= 1. / pow(r_sigma,2);
					calScaling 		+= weight * ratio;
					weightSum		+= weight; 
				}
			}
			for(unsigned int t=0;t<exptData_Target.at(i).GetData().size();++t){
				int	index_init 	= exptData_Target.at(i).GetData().at(t).GetInitialIndex();
				int	index_final 	= exptData_Target.at(i).GetData().at(t).GetFinalIndex();
				double 	calcCounts 	= EffectiveCrossSection_Target.at(i)[index_final][index_init];
				double 	exptCounts 	= exptData_Target.at(i).GetData().at(t).GetCounts();
				if(exptCounts > 0 && calcCounts > 0){
					double	sigma		= (exptData_Target.at(i).GetData().at(t).GetUpUnc() + exptData_Target.at(i).GetData().at(t).GetDnUnc())/2.;  // Average uncertainty
					double	ratio		= exptCounts / calcCounts;
					double	r_sigma		= ratio * sigma / exptCounts;
					double	weight		= 1. / pow(r_sigma,2);
					calScaling 		+= weight * ratio;
					weightSum		+= weight; 
				}
			}
		}
		calScaling /= weightSum;
		for(unsigned int ss=0;ss<scalingParameters.at(s).GetExperimentNumbers().size();ss++){
			unsigned int i 	= scalingParameters.at(s).GetExperimentNumbers().at(ss);
			scaling[i]	= calScaling;
		}
	}

	for(unsigned int i=0;i<exptData_Beam.size();i++){
		if(verbose)
			std::cout	<< std::setw(4) << std::left << i+1;
		for(unsigned int t=0;t<exptData_Beam.at(i).GetData().size();++t){
			double 	tmp 		= 0;
			int	index_init 	= exptData_Beam.at(i).GetData().at(t).GetInitialIndex();
			int	index_final 	= exptData_Beam.at(i).GetData().at(t).GetFinalIndex();
			double 	calcCounts 	= scaling.at(i) * EffectiveCrossSection_Beam.at(i)[index_final][index_init];
			double 	exptCounts 	= exptData_Beam.at(i).GetData().at(t).GetCounts();
			double	sigma		= exptData_Beam.at(i).GetData().at(t).GetUpUnc() * exptData_Beam.at(i).GetData().at(t).GetDnUnc();
			double	sigma_prime	= (exptData_Beam.at(i).GetData().at(t).GetUpUnc() - exptData_Beam.at(i).GetData().at(t).GetDnUnc());
			if(calcCounts > 0 && sigma > 0){
				if(verbose){
					if(fLikelihood){
						std::cout 	<< std::setw( 6) << std::left << index_init 
								<< std::setw( 6) << std::left << index_final 
								<< std::setw(10) << std::left << calcCounts 
								<< std::setw(10) << std::left << exptCounts 
								<< std::setw(10) << std::left << calcCounts/exptCounts
								<< std::setw(12) << std::left << 0.5 * TMath::Power((exptCounts - calcCounts),2)/(sigma + sigma_prime * (exptCounts - calcCounts));
					}
					else{
						std::cout 	<< std::setw( 6) << std::left << index_init 
								<< std::setw( 6) << std::left << index_final 
								<< std::setw(10) << std::left << calcCounts 
								<< std::setw(10) << std::left << exptCounts 
								<< std::setw(10) << std::left << calcCounts/exptCounts
								<< std::setw(12) << std::left << TMath::Power((calcCounts - exptCounts)/exptData_Beam.at(i).GetData().at(t).GetUpUnc(),2);
					}
				}
				if(fLikelihood){
					chisq		+= 0.5 * TMath::Power((exptCounts - calcCounts),2)/(sigma + sigma_prime * (exptCounts - calcCounts));
				}
				else{
					if(calcCounts > exptCounts)
						tmp 		= (calcCounts - exptCounts) / exptData_Beam.at(i).GetData().at(t).GetUpUnc();
					else
						tmp 		= (calcCounts - exptCounts) / exptData_Beam.at(i).GetData().at(t).GetDnUnc();
					chisq		+= tmp * tmp;
				}
				NDF++;
			}
		}
		if(verbose)
			std::cout << std::endl;
	}
	if(verbose && !fLikelihood){
		std::cout 	<< std::setw(7) << std::left << "Target:"
				<< std::endl;
		std::cout 	<< std::setw( 4) << std::left << "";
		for(unsigned int t=0;t<exptData_Target.at(exptData_Target.size()-1).GetData().size();t++){
			std::cout	<< std::setw( 6) << std::left << "Init:"
					<< std::setw( 6) << std::left << "Finl:"
					<< std::setw(10) << std::left << "Calc:"
					<< std::setw(10) << std::left << "Expt:"
					<< std::setw(10) << std::left << "C/E:"
					<< std::setw(12) << std::left << "Chisq:";
		}
		std::cout 	<< std::endl;
	}
	else if(verbose){
		std::cout 	<< std::setw(7) << std::left << "Target:"
				<< std::endl;
		std::cout 	<< std::setw( 4) << std::left << " ";
		for(unsigned int t=0;t<exptData_Target.at(exptData_Target.size()-1).GetData().size();t++){
			std::cout	<< std::setw( 6) << std::left << "Init:"
					<< std::setw( 6) << std::left << "Finl:"
					<< std::setw(10) << std::left << "Calc:"
					<< std::setw(10) << std::left << "Expt:"
					<< std::setw(10) << std::left << "Calc/Expt:"
					<< std::setw(12) << std::left << "-Ln(L) cont.:";
		}
		std::cout 	<< std::endl;
	}

	for(unsigned int i=0;i<exptData_Target.size();i++){
		if(verbose)
			std::cout	<< std::setw(4) << std::left << i+1;
		for(unsigned int t=0;t<exptData_Target.at(i).GetData().size();++t){
			double 	tmp 		= 0;
			int	index_init 	= exptData_Target.at(i).GetData().at(t).GetInitialIndex();
			int	index_final 	= exptData_Target.at(i).GetData().at(t).GetFinalIndex();
			double 	calcCounts 	= scaling.at(i) * EffectiveCrossSection_Target.at(i)[index_final][index_init];
			double 	exptCounts 	= exptData_Target.at(i).GetData().at(t).GetCounts();
			double	sigma		= exptData_Target.at(i).GetData().at(t).GetUpUnc() * exptData_Target.at(i).GetData().at(t).GetDnUnc();
			double	sigma_prime	= (exptData_Target.at(i).GetData().at(t).GetUpUnc() - exptData_Target.at(i).GetData().at(t).GetDnUnc());
			if(calcCounts > 0 && sigma > 0){
				if(verbose){
					if(fLikelihood){
						std::cout 	<< std::setw( 6) << std::left << index_init 
								<< std::setw( 6) << std::left << index_final 
								<< std::setw(10) << std::left << calcCounts 
								<< std::setw(10) << std::left << exptCounts 
								<< std::setw(10) << std::left << calcCounts/exptCounts
								<< std::setw(12) << std::left << 0.5 * TMath::Power((exptCounts - calcCounts),2)/(sigma + sigma_prime * (exptCounts - calcCounts));
					}
					else{
						std::cout 	<< std::setw( 6) << std::left << index_init 
								<< std::setw( 6) << std::left << index_final 
								<< std::setw(10) << std::left << calcCounts 
								<< std::setw(10) << std::left << exptCounts 
								<< std::setw(10) << std::left << calcCounts/exptCounts
								<< std::setw(12) << std::left << TMath::Power((calcCounts - exptCounts)/exptData_Target.at(i).GetData().at(t).GetUpUnc(),2);
					}
				}
				if(fLikelihood){
					chisq		+= 0.5 * TMath::Power((exptCounts - calcCounts),2)/(sigma + sigma_prime * (exptCounts - calcCounts));
				}
				else{
					if(calcCounts > exptCounts)
						tmp 		= (calcCounts - exptCounts) / exptData_Target.at(i).GetData().at(t).GetUpUnc();
					else
						tmp 		= (calcCounts - exptCounts) / exptData_Target.at(i).GetData().at(t).GetDnUnc();
					chisq		+= tmp * tmp;
				}
				NDF++;
			}
		}
		if(verbose)
			std::cout << std::endl;
	}
	if(verbose){
		std::cout 	<< std::setw(26) << std::left << "Literature chi-squared: " 
				<< std::setw(12) << std::left << std::setprecision(6) << litchisq
				<< std::setw(26) << std::left << "Total chi-squared: " 
				<< std::setw(12) << std::left << std::setprecision(6) << chisq 
				<< std::setw(26) << std::left << "NDF: " 
				<< std::setw(12) << std::left << std::setprecision(6) << NDF
				<< std::setw(26) << std::left << "Chisq / NDF: " 
				<< std::setw(12) << std::left << std::setprecision(6) << chisq / (float)NDF
				<< std::endl;
	}

	iter++;

	Clock::time_point t1 = Clock::now();
	milliseconds ms = std::chrono::duration_cast<milliseconds>(t1-t0);

	if((iter % 100) == 0 && !verbose)
		std::cout 	<< std::setw(12) << std::left << iter 
				<< std::setw(13) << std::left << chisq 
				<< std::setw(7)  << std::left << NDF
				<< std::setw(13) << std::left << chisq/(double)NDF 
				<< std::setw(12) << std::left << litchisq 
				<< std::setw(24) << std::left << ms.count() 
				<< "\r" << std::flush;

	return chisq;	

}

void GOSIASimMinFCN::ClearAll(){

	parameters.clear();			
	
	ME_Beam.clear();			
	exptData_Beam.clear();			
	litLifetimes_Beam.clear();			
	litBranchingRatios_Beam.clear();		
	litMixingRatios_Beam.clear();		
	EffectiveCrossSection_Beam.clear();

	ME_Target.clear();			
	exptData_Target.clear();			
	litLifetimes_Target.clear();			
	litBranchingRatios_Target.clear();		
	litMixingRatios_Target.clear();		
	EffectiveCrossSection_Target.clear();

}
