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

void RunFitter(const char* nuclfile = "NucleusFile.txt", const char* datafile = "Data.txt", int threads = 1){

	NucleusReader *nuclreader = new NucleusReader(nuclfile);
	Nucleus *nucl = nuclreader->GetNucleus();

	int    rmin[6] = {1,7,13,19,1,7};
	int    rmax[6] = {6,12,18,24,6,12};
	double tmin[6] = {10,10,10,10,10,10};
	double tmax[6] = {60,60,60,60,60,60};
	ParticleDetectorS3 *detectors[6];
	for(int i=0;i<6;i++){
		detectors[i] 	= new ParticleDetectorS3(i,26,-0.2,0.4,rmin[i],rmax[i],tmin[i],tmax[i]);
		if(i==0)
			detectors[i]->WriteParticleDetector("OutputFitting.root","RECREATE");
		else
			detectors[i]->WriteParticleDetector("OutputFitting.root");
	}
	
	Reaction *reac = new Reaction(76,34,208,82,287.71);
	reac->SetMass(71.927,207.9767); 
	Experiments *expts = new Experiments(nucl,reac);
	expts->SetAccuracy(1e-5);
	bool tarDet[6] = {false,false,false,false,true,true};
	for(int i=0;i<6;i++){
		expts->NewExperimentRange(detectors[i]->GetThetaMin()-2,detectors[i]->GetThetaMax()+2,11,272.148,287.71,5,tarDet[i]);
		expts->SetParticleDetectorEff(i,detectors[i]->GetThetaEfficiencyHist());
	}


	double E[7] = {270,274,278,282,286,290,294};
	double ELoss[7] = {16.969,16.965,16.960,16.954,16.947,16.939,16.929};
	StoppingPower dEdX;
	for(int e=0;e<7;e++)
		dEdX.AddStoppingPower(E[e],ELoss[e]);
	dEdX.FitStoppingPowers();
	expts->SetStopping(dEdX);
	expts->SetVerbose(false);		
	expts->SetNthreads(threads);
	expts->UseEfficiency();
	expts->PointCorrections();

	TransitionRates *rates = new TransitionRates(nucl);

	DataReader *dataReader = new DataReader(nucl,datafile);
	std::vector<ExperimentData> exptVec;
	std::vector<PointCoulEx> poinVec;
	for(unsigned int e=0;e<6;e++){
		//GammaYield::PrintYields(expts->GetExperimentRange(e),*rates,*nucl);
		ExperimentData tmpExpt = dataReader->GetExperimentData().at(e);
		tmpExpt.SetThetaCM(expts->GetExperimentRange(e).GetMeanThetaCM());
		exptVec.push_back(tmpExpt);
		PointCoulEx tmpPoin = expts->GetPointCalculation(e);
		char fname[64];
		sprintf(fname,"CoulExDetails_%i.txt",e+1);
		tmpPoin.WriteDetailsToFile(fname);
		poinVec.push_back(tmpPoin);
	}

	//expts->PrintPointCorrections();
	expts->WriteIntegralFits("OutputFitting.root");

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

}


int main(int argc, char** argv){

	if(argc == 1)
		RunFitter();
	else if(argc == 2)
		RunFitter(argv[1]);
	else if(argc == 3)
		RunFitter(argv[1],argv[2]);
	else if(argc == 4)
		RunFitter(argv[1],argv[2],std::atoi(argv[3]));
	else if(argc > 4){
		RunFitter(argv[1],argv[2]);
		std::cout << "Too many inputs, only inputs 1 & 2 used. Nucleus file: " << argv[1] << " and data file " << argv[2] << std::endl;
	}

}
