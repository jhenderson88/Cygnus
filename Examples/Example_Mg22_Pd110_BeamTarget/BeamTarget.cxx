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

#include <stdio.h>
#include <string.h>
#include <iostream>

void RunFitter(const char* beamfile, const char* targfile, int threads = 1){
	
	//******************************************************//
	//		Set up and perform integration		//
	//******************************************************//

	//	Read in the beam and target nucleus files
	NucleusReader *nuclreader_b = new NucleusReader(beamfile);
	NucleusReader *nuclreader_t = new NucleusReader(targfile);
	//	Create the beam and target nuclei from the input files
	Nucleus *nucl_b = nuclreader_b->GetNucleus();
	Nucleus *nucl_t = nuclreader_t->GetNucleus();
	
	// 	Define the beam reaction kinematics
	Reaction *reac_b = new Reaction(nucl_b->GetA(),nucl_b->GetZ(),nucl_t->GetA(),nucl_t->GetZ(),83.);
	reac_b->SetGOSIAKinematics(true);					//	Use GOSIA-style kinematics
	reac_b->SetExcitationEnergy(nucl_b->GetLevelEnergies().at(1));		//	Use the first-excited state energy of the beam for kinematics
	//	Define the experiments for the beam-like nucleus based on the nucleus file and the kinematics
	Experiments *expts_b = new Experiments(nucl_b,reac_b);
	//	Define the target reaction kinematics
	Reaction *reac_t = new Reaction(nucl_b->GetA(),nucl_b->GetZ(),nucl_t->GetA(),nucl_t->GetZ(),83.);
	reac_t->SetGOSIAKinematics(true);					//	Use GOSIA-style kinematics
	reac_t->SetExcitationEnergy(nucl_t->GetLevelEnergies().at(1));		//	Use the first-excited state energy of the target for kinematics
	//	Define the experiments for the target-like nucleus based on the nucleus file and kinematics
	Experiments *expts_t = new Experiments(nucl_t,reac_t);	

	expts_b->SetProjectileExcitation(true);		//	Set beam excitation for the beam
	expts_t->SetProjectileExcitation(false);	//	... and target excitation for the target
	expts_b->SetAccuracy(1e-7);			
	expts_t->SetAccuracy(1e-7);

	//	Define whether the scattered beam- or target-like nucleus is detected (false == beam, true == target)
	bool tarDet[6] = {false,false,false,true,true,true};
	//	Define the theta ranges (in the lab) for the experiments
	double thetamin[6] = {20.15,32.35,42.,20.15,32.35,42.};
	double thetamax[6] = {32.325,41.975,49.375,32.325,41.975,49.375};
	
	//	Set up a new ExperimentRange for each experiment, one each for beam and target
	for(int i=0;i<6;i++){
		expts_b->NewExperimentRange(thetamin[i],thetamax[i],6,73.,83.,5,tarDet[i]);
		expts_t->NewExperimentRange(thetamin[i],thetamax[i],6,73.,83.,5,tarDet[i]);
	}

	//	Define energy loss in the target
	double E[6] = {73,75,77,79,81,83};
	double ELoss[6] = {10,10,10,10,10,10};
	StoppingPower dEdX;
	for(int e=0;e<6;e++)
		dEdX.AddStoppingPower(E[e],ELoss[e]);
	dEdX.FitStoppingPowers();	// 	Fit the stopping powers to allow extrapolation

	//	Set up the Experiments for the beam and target
	//	For the beam:
	expts_b->SetStopping(dEdX);	//	Energy loss
	expts_b->SetVerbose(false);	//	Make the command-line output quiet
	expts_b->SetNthreads(threads);	//	Define the number of threads the code can use
	expts_b->UseEfficiency(false);	//	We will not be using a theta-phi efficiency distribution
	expts_b->PointCorrections();	//	Perform the point corrections
	//	... and the target:
	expts_t->SetStopping(dEdX);	//	Energy loss
	expts_t->SetVerbose(false);	//	Make the command-line output quiet                    	
	expts_t->SetNthreads(threads);  //	Define the number of threads the code can use
	expts_t->UseEfficiency(false);  //	We will not be using a theta-phi efficiency distribution
	expts_t->PointCorrections();    //	Perform the point corrections
	//	Write the full integration to file as TGraph2Ds
	expts_t->WriteIntegralFits("OutputFitting_Beam.root","RECREATE");
	expts_t->WriteIntegralFits("OutputFitting_Target.root","RECREATE");
	
	//******************************************************//
	//		Integration step complete		//
	//******************************************************//

	//	Create transition rates based on the input nuclei
	TransitionRates *rates_b = new TransitionRates(nucl_b);
	TransitionRates *rates_t = new TransitionRates(nucl_t);

	//	Read the experimental data in for the beam and target
	DataReader *dataReader_b = new DataReader(nucl_b,"Data_Mg22.txt");	// 	Beam
	DataReader *dataReader_t = new DataReader(nucl_b,"Data_Pd110.txt"); 	//	Target
	//	Create vectors of experimental data to be filled
	std::vector<ExperimentData> exptVec_b;	//	Beam
	std::vector<ExperimentData> exptVec_t;	// 	Target	
	//	Create vectors of point calculations 
	std::vector<PointCoulEx> poinVec_b;	//	Beam
	std::vector<PointCoulEx> poinVec_t;	//	Target

	//	Fill experimental data for the beam	
	for(unsigned int e=0;e<6;e++){
		ExperimentData tmpExpt = dataReader_b->GetExperimentData().at(e);	//	Temporary object to hold the data for experiment e
		tmpExpt.SetThetaCM(expts_b->GetExperimentRange(e).GetMeanThetaCM());	//	Set the CM angle (presently unused)
		exptVec_b.push_back(tmpExpt);						//	Add the temporary data object to the beam vector
		PointCoulEx tmpPoin = expts_b->GetPointCalculation(e);			//	Grab the point calculation used to calculate the corrections
		char fname[64];
		sprintf(fname,"CoulExDetails_Beam_%i.txt",e+1);
		tmpPoin.WriteDetailsToFile(fname);					//	Write the point calculation details to file
		poinVec_b.push_back(tmpPoin);						//	Push the point calculation to the beam vector
	}
	//	Fill experimental data for the target
	for(unsigned int e=0;e<6;e++){
		ExperimentData tmpExpt = dataReader_t->GetExperimentData().at(e);	//	Temporary object to hold the data for experiment e
		tmpExpt.SetThetaCM(expts_t->GetExperimentRange(e).GetMeanThetaCM());    //	Set the CM angle (presently unused)
		tmpExpt.Print();                                                        //	Add the temporary data object to the target vector
		exptVec_t.push_back(tmpExpt);                                           //	Grab the point calculation used to calculate the corrections
		PointCoulEx tmpPoin = expts_t->GetPointCalculation(e);                                                                                       
		char fname[64];                                                                                                                              
		sprintf(fname,"CoulExDetails_Target_%i.txt",e+1);                       //	Write the point calculation details to file
		tmpPoin.WriteDetailsToFile(fname);                                      //	Push the point calculation to the target vector
		poinVec_t.push_back(tmpPoin);
	}
	//	Print the Rutherford cross sections (useful for debugging vs GOSIA)
	std::cout	<< "Rutherford cross sections:"
			<< std::endl;
	std::cout	<< std::setw(15) << std::left << "Experiment:"
			<< std::setw(15) << std::left << "Beam:"
			<< std::setw(15) << std::left << "Target:"
			<< std::endl;
	for(unsigned int e=0;e<6;e++){
		std::cout	<< std::setw(15) << std::left << e+1
				<< std::setw(15) << std::left << expts_b->GetExperimentRange(e).IntegrateRutherford()	//	Beam integrated Rutherford
				<< std::setw(15) << std::left << expts_t->GetExperimentRange(e).IntegrateRutherford()	//	Target integration Rutherford
				<< std::endl;
	}
	
	//	Print the point corrections calculated previously
	std::cout	<< "\n Beam Point Corrections:"
			<< std::endl;
	expts_b->PrintPointCorrections();
	std::cout	<< "\n Target Point Corrections:"
			<< std::endl;
	expts_t->PrintPointCorrections();

	//******************************************************//
	//		Begin fitting experimental data		//
	//******************************************************//
	
	//	Create a simultaneous fitter (similar to CoulExFitter, but allows for beam and target)
	CoulExSimFitter *fitter = new CoulExSimFitter();
	fitter->SetBeamNucleus(nucl_b);		//	Pass the fitter the beam nucleus
	fitter->SetTargetNucleus(nucl_t);	//	... and the target nucleus
	for(unsigned int e=0;e<6;e++)
		fitter->AddBeamCorrectionFactor(expts_b->GetCorrectionFactors().at(e));		//	Give the fitter the beam correction factors
	for(unsigned int e=0;e<6;e++)
		fitter->AddTargetCorrectionFactor(expts_t->GetCorrectionFactors().at(e));	//	... and the target correction factors
	fitter->SetBeamData(exptVec_b);		//	Give the fitter the beam data
	fitter->SetTargetData(exptVec_t);	//	... and the target data
	for(size_t i = 0;i<poinVec_t.size();i++)
		poinVec_t.at(i).SetVerbose(false);	//	Turn off the verbose option (keep things quiet)
	for(size_t i = 0;i<poinVec_b.size();i++)
		poinVec_b.at(i).SetVerbose(false);	//	... and for the target
	fitter->SetBeamPointCalcs(poinVec_b);		//	Give the fitter the beam point calculations
	fitter->SetTargetPointCalcs(poinVec_t);		//	... and the target point calculations

	//	Define literature limits for some of the target data (first excited state lifetime)
	fitter->AddTargetLifetime(1,63.48,1.01);			

	//	Tell the fitter which matrix elements we will be fitting in the beam
	fitter->AddBeamFittingMatrixElement(1,0,1,0.46,0.3,0.8);
	fitter->AddBeamFittingMatrixElement(1,1,1,-0.500,-1,1);
	//	... and in the target
	fitter->AddTargetFittingMatrixElement(1,0,1,0.9,0.8,1.05);
	fitter->AddTargetFittingMatrixElement(1,1,1,-0.5,-1.2,0.);
	//	Define a vector of integers to define common scaling
	std::vector<int> tmpVec;
	for(unsigned int e=0;e<6;e++)
		tmpVec.push_back((int)e);
	fitter->CreateScalingParameter(tmpVec,1000,0.01,1e8);	//	Set the scaling parameters and their allowed range
	fitter->SetNthreads(threads);				//	Define the number of threads to be used
	fitter->SetVerbose(true);				//	Make the fit verbose

	//	Do the Fit:
	fitter->SetMaxIterations(5000);
	fitter->SetMaxFunctionCalls(5000);
	fitter->SetTolerance(0.01);
	std::cout << "First minimization" << std::endl;
	fitter->DoFit("Minuit2","Migrad");
	
}


int main(int argc, char** argv){

	if(argc > 0 && strcmp(argv[1],"--h") == 0 || strcmp(argv[1],"--help") == 0){
		std::cout	<< "Arguments: beam nucleus file, target nucleus file, nCores"
				<< std::endl;
		return 0;
	}

	RunFitter(argv[1],argv[2],std::atoi(argv[3]));

}
