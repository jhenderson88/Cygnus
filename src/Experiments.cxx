#include "Experiments.h"

Experiments::Experiments() 
{ 

	fAccuracy	= 1e-5;

	fProjectileExcitation 	= true;
	fUseEfficiency	= false;
	verbose		= false; 
	nThreads	= 1;

}

Experiments::Experiments(Nucleus* nucl, Reaction* reac) 
{

	fAccuracy	= 1e-5;

	fProjectileExcitation 	= true;
	fUseEfficiency	= false;
	verbose		= false;
	nThreads	= 1;

	fNucleus	= *nucl;
	fReaction	= *reac;

	experimentRanges.clear();
	pointCalculation.clear();

}

Experiments::Experiments(const Experiments &e) : 
experimentRanges(e.experimentRanges.size()), pointCalculation(e.pointCalculation.size())
{

	fProjectileExcitation 	= e.fProjectileExcitation;

	fUseEfficiency	= e.fUseEfficiency;

	fStopping	= e.fStopping;	

	fAccuracy	= e.fAccuracy;

	nThreads	= e.nThreads;

	verbose		= e.verbose;

	fNucleus	= e.fNucleus;
	fReaction	= e.fReaction;	

	integratedCrossSections.resize(e.integratedCrossSections.size());
	pointCrossSections.resize(e.pointCrossSections.size());
	correctionFactors.resize(e.correctionFactors.size());

	for(size_t i=0; i<e.experimentRanges.size() ;i++)
		experimentRanges.at(i) = e.experimentRanges.at(i); 
	for(size_t i=0; i<e.pointCalculation.size() ;i++)
		pointCalculation.at(i) = e.pointCalculation.at(i);

	for(size_t i=0; i<e.integratedCrossSections.size(); i++){
		integratedCrossSections.at(i).ResizeTo(e.integratedCrossSections.at(i).GetNrows());
		integratedCrossSections.at(i) = e.integratedCrossSections.at(i);
	}	
	for(size_t i=0; i<e.pointCrossSections.size(); i++){
		pointCrossSections.at(i).ResizeTo(e.pointCrossSections.at(i).GetNrows());
		pointCrossSections.at(i) = e.pointCrossSections.at(i);
	}	
	for(size_t i=0; i<e.correctionFactors.size(); i++){
		correctionFactors.at(i).ResizeTo(e.correctionFactors.at(i).GetNrows());
		correctionFactors.at(i) = e.correctionFactors.at(i);
	}	


}
Experiments& Experiments::operator = (const Experiments& e){

	fProjectileExcitation 	= e.fProjectileExcitation;

	fUseEfficiency	= e.fUseEfficiency;

	fStopping	= e.fStopping;
	
	fAccuracy	= e.fAccuracy;

	nThreads	= e.nThreads;

	verbose		= e.verbose;

	fNucleus	= e.fNucleus;
	fReaction	= e.fReaction;

	experimentRanges.resize(e.experimentRanges.size());
	pointCalculation.resize(e.pointCalculation.size());
	integratedCrossSections.resize(e.integratedCrossSections.size());
	pointCrossSections.resize(e.pointCrossSections.size());
	correctionFactors.resize(e.correctionFactors.size());

	for(size_t i=0; i<e.experimentRanges.size() ;i++)
		experimentRanges.at(i) = e.experimentRanges.at(i); 
	for(size_t i=0; i<e.pointCalculation.size() ;i++)
		pointCalculation.at(i) = e.pointCalculation.at(i);

	for(size_t i=0; i<e.integratedCrossSections.size(); i++){
		integratedCrossSections.at(i).ResizeTo(e.integratedCrossSections.at(i).GetNrows());
		integratedCrossSections.at(i) = e.integratedCrossSections.at(i);
	}	
	for(size_t i=0; i<e.pointCrossSections.size(); i++){
		pointCrossSections.at(i).ResizeTo(e.pointCrossSections.at(i).GetNrows());
		pointCrossSections.at(i) = e.pointCrossSections.at(i);
	}	
	for(size_t i=0; i<e.correctionFactors.size(); i++){
		correctionFactors.at(i).ResizeTo(e.correctionFactors.at(i).GetNrows());
		correctionFactors.at(i) = e.correctionFactors.at(i);
	}	

	return *this;

} 

void Experiments::NewExperimentRange(double thetamin, double thetamax, int nT, double energymin, double energymax, int nE, bool targetDetection){

	ExperimentRange tmp(&fNucleus,&fReaction);
	tmp.SetRanges(thetamin,thetamax,nT,energymin,energymax,nE);
	tmp.SetTargetDetection(targetDetection);
	experimentRanges.push_back(tmp);
	
	Reaction tmpReac = fReaction;
	Nucleus tmpNucl = fNucleus;
	tmpReac.SetLabEnergy(tmp.GetMeanEnergy());
	PointCoulEx tmpPoin(&tmpNucl,&tmpReac);
	pointCalculation.push_back(tmpPoin);

}

void Experiments::SetStopping(){

	if(fStopping.StoppingSize()==0){
		std::cout	<< "Stopping powers not set!"
				<< std::endl;
		return;
	}
	for(size_t c = 0;c<experimentRanges.size();c++)
		experimentRanges.at(c).SetStopping(fStopping);
}

void Experiments::PointCorrections(){

	if(fStopping.StoppingSize()==0){
		std::cout	<< "Stopping powers not set, cannot complete integration"
				<< std::endl;
		return;
	}

	correctionFactors.resize(experimentRanges.size());

	std::cout	<< std::setw(18) << std::left << "Theta CM min:" 
			<< std::setw(18) << std::left << "Theta CM max:"
			<< std::setw(18) << std::left << "Theta Lab min:" 
			<< std::setw(18) << std::left << "Theta Lab max:"
			<< std::setw(14) << std::left << "N-threads:"
			<< std::setw(20) << std::left << "Processing time [ms]:"
			<< std::endl;

	for(size_t c = 0;c<experimentRanges.size();c++){

		experimentRanges.at(c).SetProjectileExcitation(fProjectileExcitation);
		experimentRanges.at(c).SetNthreads(nThreads);
		experimentRanges.at(c).SetAccuracy(fAccuracy);
		experimentRanges.at(c).SetStopping(fStopping);
		//experimentRanges.at(c).UseEfficiency(fUseEfficiency);
		experimentRanges.at(c).IntegrateRange();
		TVectorD tmpVec_int;
		tmpVec_int.ResizeTo(experimentRanges.at(c).GetIntegratedCrossSection_TVec().GetNrows());
		tmpVec_int = experimentRanges.at(c).GetIntegratedCrossSection_TVec();

		pointCalculation.at(c).SetProjectileExcitation(fProjectileExcitation);
		pointCalculation.at(c).SetVerbose(verbose);
		pointCalculation.at(c).SetAccuracy(fAccuracy);
		pointCalculation.at(c).CalculatePointProbabilities(experimentRanges.at(c).GetMeanThetaCM());
		TVectorD tmpVec;
		tmpVec.ResizeTo(pointCalculation.at(c).GetProbabilitiesVector().GetNrows());
		tmpVec = pointCalculation.at(c).GetProbabilitiesVector();
		tmpVec *= pointCalculation.at(c).GetReaction()->RutherfordCM(experimentRanges.at(c).GetMeanThetaCM()); 

		integratedCrossSections.push_back(tmpVec_int);
		pointCrossSections.push_back(tmpVec);

		correctionFactors.at(c).ResizeTo(tmpVec.GetNrows());
		for(int s = 0; s<pointCrossSections.at(c).GetNrows(); s++){
			correctionFactors.at(c)[s] = integratedCrossSections.at(c)[s] / pointCrossSections.at(c)[s];
		}

	}

}

void Experiments::PrintDetails() const{

	std::cout << experimentRanges.size() << " experiments defined" << std::endl;

	for(unsigned int i=0;i<experimentRanges.size();i++){
		experimentRanges.at(i).PrintDetails();
		pointCrossSections.at(i).Print();
	}

}

void Experiments::PrintPointCorrections() {

	for(size_t i = 0; i < experimentRanges.size(); i++){
		std::cout 	<< std::setw(10) << std::left << "Experiment " << i+1 << std::endl;
		std::cout	<< std::setw(16) << std::left << "Mean Theta CM:" << experimentRanges.at(i).GetMeanThetaCM() << std::endl;
		std::cout 	<< std::setw(8)  << std::left << "State:"
				<< std::setw(10) << std::left << "Corr'n:" 
				<< std::setw(14) << std::left << "Integral:"
				<< std::setw(14) << std::left << "Point:"
				<< std::setw(20) << std::left << "Corrected Point" 
				<< std::endl;
		for(int s =0; s<correctionFactors.at(i).GetNrows(); s++){
			std::cout 	<< std::setw(8)  << std::left << s+1
					<< std::setw(10) << std::left << correctionFactors.at(i)[s]
					<< std::setw(14) << std::left << integratedCrossSections.at(i)[s]
					<< std::setw(14) << std::left << pointCrossSections.at(i)[s] 
					<< std::setw(20) << std::left << pointCalculation.at(i).GetProbabilitiesVector()[s] * correctionFactors.at(i)[s] * pointCalculation.at(i).GetReaction()->RutherfordCM(experimentRanges.at(i).GetMeanThetaCM())
					<< std::endl;
		}
	}
}

void Experiments::WriteIntegralFits(const char* outfilename, const char* option) {

	TFile *outfile = new TFile(outfilename,option);
	TDirectory *dir;
	if(outfile->GetDirectory("Cross_Section_Distributions"))
		dir = outfile->GetDirectory("Cross_Section_Distributions");
	else
		dir  = outfile->mkdir("Cross_Section_Distributions");
	TDirectory *sub_dir[experimentRanges.size()];
	dir->cd();
	for(unsigned int i=0;i < experimentRanges.size(); i++){
		char dname[32];
		sprintf(dname,"Experiment_%i",i+1);
		sub_dir[i] = dir->mkdir(dname);
		sub_dir[i]->cd();
		TGraph2D gRuth = *experimentRanges.at(i).GetRutherfordThetaEnergy();	
		char rname[64];
		sprintf(rname,"Rutherford_Theta_Energy_Experiment_%i",i+1);
		gRuth.SetName(rname);
		gRuth.Write();
		for(int s=0; s < experimentRanges.at(i).GetIntegratedCrossSection_TVec().GetNrows(); s++){
			TGraph2D g2 = *experimentRanges.at(i).InterpolatedEnergyTheta(s);
			char gname[64];
			sprintf(gname,"Fine_Theta_Energy_CS_Expt_%i_State_%i",i+1,s);
			g2.SetName(gname);
			TGraph2D g2_2 = *experimentRanges.at(i).GetThetaEnergyGraph(s);
			sprintf(gname,"Theta_Energy_CS_Expt_%i_State_%i",i+1,s);
			g2_2.SetName(gname);
			TGraph2D g2_3 = *experimentRanges.at(i).InterpolatedEnergyTheta(s,true);
			sprintf(gname,"Efficiency_Modified_CS_Expt_%i_State_%i",i+1,s);
			g2_3.SetName(gname);
			g2.Write();
			g2_2.Write();
			g2_3.Write();

		}
		dir->cd();
	}
	outfile->cd();
	outfile->Close();

}
