#include "TMatrixD.h"
#include "TVectorD.h"

#include "ParticleDetector.h"
#include "ParticleDetectorS3.h"
#include "ScalingParameter.h"
#include "CoulExFitter.h"
#include "GammaYield.h"
#include "PointCoulEx.h"
#include "Experiments.h"
#include "TransitionRates.h"
#include "DataReader.h"
#include "Nucleus.h"
#include "NucleusReader.h"
#include "StoppingPower.h"

#include <iostream>

void RunFitter(const char* beamfile, const char* targfile, int threads = 1){

	NucleusReader *nuclreader_b = new NucleusReader(beamfile);
	NucleusReader *nuclreader_t = new NucleusReader(targfile);
	Nucleus *nucl_b = nuclreader_b->GetNucleus();
	Nucleus *nucl_t = nuclreader_t->GetNucleus();

	nucl_b->PrintNucleus();
	nucl_t->PrintNucleus();

	int    rmin[6] = {1,9,17,1,9,17};
	int    rmax[6] = {8,16,24,8,16,24};
	double tmin[6] = {10,10,10,10,10,10};
	double tmax[6] = {60,60,60,60,60,60};
	ParticleDetectorS3 *detectors[6];
	for(int i=0;i<6;i++){
		detectors[i] = new ParticleDetectorS3(i,30,0.,0.,rmin[i],rmax[i],tmin[i],tmax[i]);
		if(i==0)
			detectors[i]->WriteParticleDetector("OutputFitting.root","RECREATE");
		else
			detectors[i]->WriteParticleDetector("OutputFitting.root");
		std::cout << detectors[i]->GetThetaMin() << "\t" << detectors[i]->GetThetaMax() << std::endl;
	}
	
	Reaction *reac = new Reaction(nucl_b->GetA(),nucl_b->GetZ(),nucl_t->GetA(),nucl_t->GetZ(),83.);
	//reac->SetMass(21.9996,109.90952); 
	reac->SetExcitationEnergy(nucl_b->GetLevelEnergies().at(1));
	reac->SetGOSIAKinematics(false);
	Experiments *expts_b = new Experiments(nucl_b,reac);
	reac->SetExcitationEnergy(nucl_t->GetLevelEnergies().at(1));
	Experiments *expts_t = new Experiments(nucl_t,reac);
	expts_b->SetProjectileExcitation(true);
	expts_t->SetProjectileExcitation(false);
	expts_b->SetAccuracy(1e-8);
	expts_t->SetAccuracy(1e-8);
	bool tarDet[6] = {false,false,false,true,true,true};
	double thetamin[6] = {20.15,32.35,42.,20.15,32.35,42.};
	double thetamax[6] = {32.325,41.975,49.375,32.325,41.975,49.375};
	for(int i=0;i<6;i++){
		//expts_b->NewExperimentRange(thetamin[i],thetamax[i],6,73.,83.,5,tarDet[i]);
		//expts_t->NewExperimentRange(thetamin[i],thetamax[i],6,73.,83.,5,tarDet[i]);
		expts_b->NewExperimentRange(detectors[i]->GetThetaMin()-3,detectors[i]->GetThetaMax()+3,6,73.,83.,5,tarDet[i]);
		expts_b->SetParticleDetectorEff(i,detectors[i]->GetThetaEfficiencyGraph());
		expts_t->NewExperimentRange(detectors[i]->GetThetaMin()-3,detectors[i]->GetThetaMax()+3,6,73.,83.,5,tarDet[i]);
		expts_t->SetParticleDetectorEff(i,detectors[i]->GetThetaEfficiencyGraph());
	}

	double E[6] = {73,75,77,79,81,83};
	double ELoss[6] = {10,10,10,10,10,10};
	StoppingPower dEdX;
	for(int e=0;e<6;e++)
		dEdX.AddStoppingPower(E[e],ELoss[e]);
	dEdX.FitStoppingPowers();

	expts_b->SetStopping(dEdX);
	expts_b->SetVerbose(false);		
	expts_b->SetNthreads(threads);
	expts_b->UseEfficiency(false);
	expts_b->PointCorrections();
	expts_t->SetStopping(dEdX);
	expts_t->SetVerbose(false);		
	expts_t->SetNthreads(threads);
	expts_t->UseEfficiency(false);
	expts_t->PointCorrections();
	expts_t->WriteIntegralFits("OutputFitting_Beam.root","RECREATE");
	expts_t->WriteIntegralFits("OutputFitting_Target.root","RECREATE");

	TransitionRates *rates_b = new TransitionRates(nucl_b);
	TransitionRates *rates_t = new TransitionRates(nucl_t);

	std::vector<PointCoulEx> poinVec_b;
	for(unsigned int e=0;e<6;e++){
		GammaYield::PrintYields(expts_b->GetExperimentRange(e),*rates_b,*nucl_b);
		PointCoulEx tmpPoin = expts_b->GetPointCalculation(e);
		char fname[64];
		sprintf(fname,"CoulExDetails_Beam_%i.txt",e+1);
		tmpPoin.WriteDetailsToFile(fname);
		poinVec_b.push_back(tmpPoin);
	}
	std::vector<PointCoulEx> poinVec_t;
	for(unsigned int e=0;e<6;e++){
		GammaYield::PrintYields(expts_t->GetExperimentRange(e),*rates_t,*nucl_t);
		PointCoulEx tmpPoin = expts_t->GetPointCalculation(e);
		char fname[64];
		sprintf(fname,"CoulExDetails_Target_%i.txt",e+1);
		tmpPoin.WriteDetailsToFile(fname);
		poinVec_t.push_back(tmpPoin);
	}
	std::cout	<< "Rutherford cross sections:"
			<< std::endl;
	std::cout	<< std::setw(15) << std::left << "Experiment:"
			<< std::setw(15) << std::left << "Beam:"
			<< std::setw(15) << std::left << "Target:"
			<< std::endl;
	for(unsigned int e=0;e<6;e++){
		std::cout	<< std::setw(15) << std::left << e+1
				<< std::setw(15) << std::left << expts_b->GetExperimentRange(e).IntegrateRutherford()
				<< std::setw(15) << std::left << expts_t->GetExperimentRange(e).IntegrateRutherford()
				<< std::endl;
	}

	return;

	//expts->PrintPointCorrections();
	/*expts->WriteIntegralFits("OutputFitting.root");

	CoulExFitter *fitter = new CoulExFitter();
	fitter->SetNucleus(nucl);
	for(unsigned int e=0;e<6;e++)
		fitter->AddCorrectionFactor(expts->GetCorrectionFactors().at(e));
	fitter->SetData(exptVec);
	for(size_t i = 0;i<poinVec.size();i++)
		poinVec.at(i).SetVerbose(false);
	fitter->SetPointCalcs(poinVec);
	fitter->AddLifetime(1,4.07,0.289);
	fitter->AddFittingMatrixElement(1,0,1,0.46,0.44,0.48);
	fitter->AddFittingMatrixElement(1,1,1,-0.500,-1,-0.1);
	std::vector<int> tmpVec;
	for(unsigned int e=0;e<6;e++)
		tmpVec.push_back((int)e);
	fitter->CreateScalingParameter(tmpVec,8,0.01,10);
	fitter->SetNthreads(threads);
	//fitter->Print();
	fitter->SetVerbose(true);

	fitter->SetMaxIterations(1000);
	fitter->SetMaxFunctionCalls(1000);
	fitter->SetTolerance(0.0001);
	std::cout << "First minimization" << std::endl;
	//fitter->DoFit("Genetic","");
	fitter->DoFit("Minuit2","Migrad");
	//std::cout << "Second minimization" << std::endl;
	//fitter->DoFit("Genetic","ConjugateFR");
	//fitter->DoFit("GSLMultiMin","ConjugatePR");
	*/
}


int main(int argc, char** argv){

	RunFitter(argv[1],argv[2],std::atoi(argv[3]));


}
