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

void RunFitter(const char* nuclfile = "NucleusFile_Mg22.txt", int threads = 1){

	NucleusReader *nuclreader = new NucleusReader(nuclfile);
	Nucleus *nucl = nuclreader->GetNucleus();
	
	Reaction *reac = new Reaction(22,12,110,46,78);
	reac->SetGOSIAKinematics(true);
	Experiments *expts = new Experiments(nucl,reac);
	expts->UseEfficiency(false);
	expts->SetAccuracy(1e-5);
	bool tarDet[6] = {false,false,false,true,true,true};
	double theta_in[6]	= {20.15,32.35,42.,20.15,32.35,42.};
	double theta_out[6]	= {32.325,41.975,49.375,32.325,41.975,49.375};
	for(int i=0;i<6;i++){
		expts->NewExperimentRange(theta_in[i],theta_out[i],11,73,83,5,tarDet[i]);
	}

	double E[6] = {73,75,77,79,81,83};
	double ELoss[6] = {10,10,10,10,10,10};
	StoppingPower dEdX;
	for(int e=0;e<6;e++)
		dEdX.AddStoppingPower(E[e],ELoss[e]);
	dEdX.FitStoppingPowers();

	std::cout	<< "Point corrections:"
			<< std::endl;
	expts->SetStopping(dEdX);
	expts->SetVerbose(false);		
	expts->SetNthreads(threads);
	expts->PointCorrections();

	std::cout	<< "Rutherford cross-sections (integrated):"
			<< std::endl;
	for(unsigned int e=0;e<6;e++)
		std::cout	<< std::setw(12) << std::left << e+1
				<< std::setw(16) << std::left << expts->GetExperimentRange(e).GetIntegratedRutherford()
				<< std::endl;
	expts->WriteIntegralFits("OutputFitting.root","RECREATE");

	TransitionRates *rates = new TransitionRates(nucl);
	rates->Print();

	expts->PrintPointCorrections();

	
	for(unsigned int e=0;e<6;e++){
		PointCoulEx tmpPoin = expts->GetPointCalculation(e);
		char fname[64];
		sprintf(fname,"CoulExDetails_%i.txt",e+1);
		tmpPoin.WriteDetailsToFile(fname);
	}


}


int main(int argc, char** argv){

	if(argc == 1)
		RunFitter();
	else if(argc == 2)
		RunFitter(argv[1]);
	else if(argc == 3)
		RunFitter(argv[1],std::atoi(argv[2]));
	else if(argc > 3){
		RunFitter(argv[1],std::atoi(argv[2]));
		std::cout << "Too many inputs, only inputs 1 & 2 used. Nucleus file: " << argv[1] << " and data file " << argv[2] << std::endl;
	}

}
