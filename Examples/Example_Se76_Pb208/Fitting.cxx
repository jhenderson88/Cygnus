#include "TMatrixD.h"
#include "TVectorD.h"

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

	nucl->PrintNucleus();
	
	Reaction *reac = new Reaction(76,34,208,82,309.7);
	reac->SetMass(75.9192,207.9767); 
	Experiments *expts = new Experiments(nucl,reac);
	expts->SetAccuracy(1e-5);
	double tmin[6] = {20,30,40,20,30,40};
	double tmax[6] = {30,40,50,30,40,50};
	bool tarDet[6] = {false,false,false,true,true,true};
	for(int i=0;i<6;i++){
		expts->NewExperimentRange(tmin[i],tmax[i],11,289.7,309.7,5,tarDet[i]);
	}
	double E[8] = {280.,285.,290.,295.,300.,305.,310.,315.};
	double ELoss[8] = {16.973,16.964,16.954,16.942,16.928,16.911,16.892,16.870};
	StoppingPower dEdX;
	for(int e=0;e<8;e++)
		dEdX.AddStoppingPower(E[e],ELoss[e]);
	dEdX.FitStoppingPowers();
	expts->SetStopping(dEdX);
	expts->SetVerbose(false);		
	expts->SetNthreads(threads);
	expts->PointCorrections();

	TransitionRates *rates = new TransitionRates(nucl);
	rates->Print();

	DataReader *dataReader = new DataReader(nucl,datafile);
	std::cout << "\n" << dataReader->GetExperimentData().size() << " experiments read in\n" << std::endl;
	std::vector<ExperimentData> exptVec;
	std::vector<PointCoulEx> poinVec;
	for(unsigned int e=0;e<6;e++){
		GammaYield::PrintYields(expts->GetExperimentRange(e),*rates,*nucl);
		ExperimentData tmpExpt = dataReader->GetExperimentData().at(e);
		tmpExpt.SetThetaCM(expts->GetExperimentRange(e).GetMeanThetaCM());
		exptVec.push_back(tmpExpt);
		PointCoulEx tmpPoin = expts->GetPointCalculation(e);
		char fname[64];
		sprintf(fname,"CoulExDetails_%i.txt",e+1);
		tmpPoin.WriteDetailsToFile(fname);
		poinVec.push_back(tmpPoin);
	}
	std::cout << std::endl;
	expts->PrintPointCorrections();
	expts->WriteIntegralFits("CrossSections.root","RECREATE");


	CoulExFitter *fitter = new CoulExFitter();
	fitter->SetNucleus(nucl);
	for(unsigned int e=0;e<6;e++)
		fitter->AddCorrectionFactor(expts->GetCorrectionFactors().at(e));
	fitter->SetData(exptVec);
	for(size_t i = 0;i<poinVec.size();i++)
		poinVec.at(i).SetVerbose(false);
	fitter->SetPointCalcs(poinVec);
	fitter->AddLifetime(1,17.745,0.289);
	fitter->AddLifetime(3,4.934,0.289);
	fitter->AddLifetime(4,2.193,0.072);
	fitter->AddBranchingRatio(3,0,1,0.584,0.032);
	fitter->AddFittingMatrixElement(1,0,1,0.486,0.320,1.270);
	fitter->AddFittingMatrixElement(1,0,3,0.571,0.05,0.2);
	fitter->AddFittingMatrixElement(1,1,1,-0.200,-1,-0.1);
	fitter->AddFittingMatrixElement(1,1,2,0.6999,0.1,1);
	fitter->AddFittingMatrixElement(1,1,3,0.412,0.1,1.8);
	fitter->AddFittingMatrixElement(1,1,4,1.7054,0.2,2.2);
	std::vector<int> tmpVec;
	for(unsigned int e=0;e<6;e++)
		tmpVec.push_back((int)e);
	fitter->CreateScalingParameter(tmpVec,0.1375,0.01,1e6);
	fitter->SetNthreads(threads);
	fitter->Print();
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
