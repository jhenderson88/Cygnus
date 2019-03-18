#include "TMatrixD.h"
#include "TVectorD.h"

#include "ParticleDetector.h"
#include "ParticleDetectorS3.h"
#include "ScalingParameter.h"
#include "CoulExFitter.h"
#include "CoulExSimFitter.h"
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
	
	Reaction *reac_b = new Reaction(nucl_b->GetA(),nucl_b->GetZ(),nucl_t->GetA(),nucl_t->GetZ(),83.);
	reac_b->SetGOSIAKinematics(true);
	reac_b->SetExcitationEnergy(nucl_b->GetLevelEnergies().at(1));
	Experiments *expts_b = new Experiments(nucl_b,reac_b);
	Reaction *reac_t = new Reaction(nucl_b->GetA(),nucl_b->GetZ(),nucl_t->GetA(),nucl_t->GetZ(),83.);
	reac_t->SetGOSIAKinematics(true);
	reac_t->SetExcitationEnergy(nucl_t->GetLevelEnergies().at(1));
	Experiments *expts_t = new Experiments(nucl_t,reac_t);
	expts_b->SetProjectileExcitation(true);
	expts_t->SetProjectileExcitation(false);
	expts_b->SetAccuracy(1e-7);
	expts_t->SetAccuracy(1e-7);
	bool tarDet[6] = {false,false,false,true,true,true};
	double thetamin[6] = {20.15,32.35,42.,20.15,32.35,42.};
	double thetamax[6] = {32.325,41.975,49.375,32.325,41.975,49.375};
	for(int i=0;i<6;i++){
		expts_b->NewExperimentRange(thetamin[i],thetamax[i],6,73.,83.,5,tarDet[i]);
		expts_t->NewExperimentRange(thetamin[i],thetamax[i],6,73.,83.,5,tarDet[i]);
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

	rates_t->Print();

	DataReader *dataReader_b = new DataReader(nucl_b,"Data_Mg22.txt");
	DataReader *dataReader_t = new DataReader(nucl_b,"Data_Pd110.txt"); 
	std::vector<ExperimentData> exptVec_b;
	std::vector<ExperimentData> exptVec_t;

	std::vector<PointCoulEx> poinVec_b;
	for(unsigned int e=0;e<6;e++){
		//GammaYield::PrintYields(expts_b->GetExperimentRange(e),*rates_b,*nucl_b);
		ExperimentData tmpExpt = dataReader_b->GetExperimentData().at(e);
		tmpExpt.SetThetaCM(expts_b->GetExperimentRange(e).GetMeanThetaCM());
		tmpExpt.Print();
		exptVec_b.push_back(tmpExpt);
		PointCoulEx tmpPoin = expts_b->GetPointCalculation(e);
		char fname[64];
		sprintf(fname,"CoulExDetails_Beam_%i.txt",e+1);
		tmpPoin.WriteDetailsToFile(fname);
		poinVec_b.push_back(tmpPoin);
	}
	std::vector<PointCoulEx> poinVec_t;
	for(unsigned int e=0;e<6;e++){
		//GammaYield::PrintYields(expts_t->GetExperimentRange(e),*rates_t,*nucl_t);
		ExperimentData tmpExpt = dataReader_t->GetExperimentData().at(e);
		tmpExpt.SetThetaCM(expts_t->GetExperimentRange(e).GetMeanThetaCM());
		tmpExpt.Print();
		exptVec_t.push_back(tmpExpt);
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

	std::cout	<< "\n Beam Point Corrections:"
			<< std::endl;
	expts_b->PrintPointCorrections();
	std::cout	<< "\n Target Point Corrections:"
			<< std::endl;
	expts_t->PrintPointCorrections();

	CoulExSimFitter *fitter = new CoulExSimFitter();
	fitter->SetBeamNucleus(nucl_b);
	fitter->SetTargetNucleus(nucl_t);
	for(unsigned int e=0;e<6;e++)
		fitter->AddBeamCorrectionFactor(expts_b->GetCorrectionFactors().at(e));
	for(unsigned int e=0;e<6;e++)
		fitter->AddTargetCorrectionFactor(expts_t->GetCorrectionFactors().at(e));
	fitter->SetBeamData(exptVec_b);
	fitter->SetTargetData(exptVec_t);
	for(size_t i = 0;i<poinVec_t.size();i++)
		poinVec_t.at(i).SetVerbose(false);
	for(size_t i = 0;i<poinVec_b.size();i++)
		poinVec_b.at(i).SetVerbose(false);
	fitter->SetBeamPointCalcs(poinVec_b);
	fitter->SetTargetPointCalcs(poinVec_t);
	fitter->AddTargetLifetime(1,63.48,1.01);
	fitter->AddBeamFittingMatrixElement(1,0,1,0.46,0.3,0.8);
	fitter->AddBeamFittingMatrixElement(1,1,1,-0.500,-1,1);
	fitter->AddTargetFittingMatrixElement(1,0,1,0.9,0.8,1.05);
	fitter->AddTargetFittingMatrixElement(1,1,1,-0.5,-1.2,0.);
	std::vector<int> tmpVec;
	for(unsigned int e=0;e<6;e++)
		tmpVec.push_back((int)e);
	fitter->CreateScalingParameter(tmpVec,10.1,0.01,1e8);
	fitter->SetNthreads(threads);
	fitter->SetVerbose(true);

	//fitter->Print();

	//return;

	fitter->SetMaxIterations(5000);
	fitter->SetMaxFunctionCalls(5000);
	fitter->SetTolerance(0.01);
	std::cout << "First minimization" << std::endl;
	fitter->DoFit("Minuit2","Migrad");
	
}


int main(int argc, char** argv){

	RunFitter(argv[1],argv[2],std::atoi(argv[3]));


}
