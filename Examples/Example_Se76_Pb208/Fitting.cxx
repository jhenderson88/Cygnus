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

	// 	Define the nucleus of interest
	NucleusReader *nuclreader = new NucleusReader(nuclfile);	//	Read from the nucleus file
	Nucleus *nucl = nuclreader->GetNucleus();			//	Construct a nucleus based on the input file
	nucl->PrintNucleus();						//	Print the output to the terminal
	
	//	Define the reaction (beam energy, etc.)	
	Reaction *reac = new Reaction(76,34,208,82,309.7);		//	Construct a Reaction. Arguments: Beam A, Beam Z, Target A, Target Z, Beam energy (MeV)
	reac->SetMass(75.9192,207.9767); 				//	Use the correct nuclear masses. Without this, Cygnus uses the A values to determine kinematics.

	//	Define the experiment.
	//	Note: ExperimentRanges in Cygnus correspond to any individual detector region, like a subset of rings on an annular silicon for which yields will be provided	
	Experiments *expts = new Experiments(nucl,reac);		//	Construct the Experiment with the previously defined Nucleus and Reaction.
	expts->SetAccuracy(1e-5);					//	Set the accuracy of the calculation. 1e-5 is rather low resolution. Default is 1e-8. Reduced accuracy leads to increased speed.
	double tmin[6] = {20,30,40,20,30,40};				//	Minimum theta (lab, degrees) of each of the six defined experimental ranges
	double tmax[6] = {30,40,50,30,40,50};				//	Maximum theta (lab, degrees) of each of the six defined experimental ranges
	bool tarDet[6] = {false,false,false,true,true,true};		//	Is the target detected in the particle detector? Defines the kinematics of the reaction.
	//	Create the ExperimentRange objects which are stored within the Experiment
	// 	Arguments: minimum theta (lab), maximum theta (lab), number of theta meshpoints, minimum energy (MeV), maximum energy (MeV), target detection bool
	for(int i=0;i<6;i++)
		expts->NewExperimentRange(tmin[i],tmax[i],11,289.7,309.7,5,tarDet[i]);

	//	Define the energy loss in the target - required for integration
	double E[8] = {280.,285.,290.,295.,300.,305.,310.,315.};
	double ELoss[8] = {16.973,16.964,16.954,16.942,16.928,16.911,16.892,16.870};
	StoppingPower dEdX;
	for(int e=0;e<8;e++)
		dEdX.AddStoppingPower(E[e],ELoss[e]);
	//	User provided stopping powers are fit with a cubic polynomial 
	dEdX.FitStoppingPowers();

	expts->SetStopping(dEdX);	//	Add the stopping powers to the experiment
	expts->SetVerbose(false);	//	Set verbocity
	expts->SetNthreads(threads);	//	Define the number of threads to be used in the integration process

	//	This performs the point corrections, in which a lot happens:
	//	1)	For each ExperimentRange, full calculations are performed at each meshpoint
	//	2) 	These meshpoints are used to fit a full theta-energy cross-section surface
	//	3)	This surface is numerically integrated in theta with dE/dX to yield integrated cross sections
	//	4)	Point calculations are performed at the mean theta and energy
	//	5)	Point calculations are compared with the full, integrated cross sections and correction factors are calculated
	//		which can be used to modify point calculations in order to capture the theta/energy range of the experiment.
	expts->PointCorrections();	

	//	TransitionRates defines the gamma-decay properties (lifetimes, branching ratios, mixing ratios) of the nucleus
	TransitionRates *rates = new TransitionRates(nucl);
	rates->Print();		//	Print the starting values to the terminal

	//	Now our starting calculations have been performed, we read in the data
	//	DataReader requires the Nucleus so it can associate yields with transitions
	DataReader *dataReader = new DataReader(nucl,datafile);
	std::cout << "\n" << dataReader->GetExperimentData().size() << " experiments read in\n" << std::endl;

	//	Vectors of ExperimentData and PointCalculations are constructed for use in the fitting
	std::vector<ExperimentData> exptVec;
	std::vector<PointCoulEx> poinVec;
	//	Loop over the experiments
	for(unsigned int e=0;e<6;e++){
		//	GammaYield::PrintYields is a static function that converts cross sections into gamma-ray yields. It is compatible with a number of arguments (see documentation).
		GammaYield::PrintYields(expts->GetExperimentRange(e),*rates,*nucl);
		//	A temporary ExperimentData object is created, to push back into the vector
		ExperimentData tmpExpt = dataReader->GetExperimentData().at(e);		//	Get the data from the DataReader object created from file
		tmpExpt.SetThetaCM(expts->GetExperimentRange(e).GetMeanThetaCM());	//	The mean center of mass angle is also stored here 
		exptVec.push_back(tmpExpt);						//	Push the temporary object into the vector
		//	A temporary PointCoulEx object is created for putting in the vector
		PointCoulEx tmpPoin = expts->GetPointCalculation(e);			//	Get the PointCalculation from the Experiment object. This is at the average theta/energy.
		char fname[64];
		sprintf(fname,"CoulExDetails_%i.txt",e+1);
		tmpPoin.WriteDetailsToFile(fname);					//	Write the details of the calculation
		poinVec.push_back(tmpPoin);						// 	Push it into the vector
	}
	std::cout << std::endl;
	expts->PrintPointCorrections();							//	Print the corrections
	expts->WriteIntegralFits("CrossSections.root","RECREATE");			//	Write the TGraph2D objects containing the cross-section surfaces to file for user inspection

	//******************************************************************************//
	//                         BEGIN THE FITTING PROCEDURE                          //
	//******************************************************************************//

	//	Create the fitter
	CoulExFitter *fitter = new CoulExFitter();
	fitter->SetNucleus(nucl);							//	Pass it a Nucleus
	for(unsigned int e=0;e<6;e++)
		fitter->AddCorrectionFactor(expts->GetCorrectionFactors().at(e));	//	Tell the fitter what the calculated correction factors are for each experiment
	fitter->SetData(exptVec);							//	Give the fitter the experimental data
	for(size_t i = 0;i<poinVec.size();i++)
		poinVec.at(i).SetVerbose(false);					//	Set the point calculation verbocity (recommended false, unless you like spam)
	fitter->SetPointCalcs(poinVec);							//	Pass the PointCoulEx objects to the fitter

	//	Define literature data for use in constraining the fits
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

	//	This section tells the fitter to use common scaling for all ExperimentRanges
	std::vector<int> tmpVec;
	for(unsigned int e=0;e<6;e++)
		tmpVec.push_back((int)e);
	// 	Arguments: vector containing experiments to be scaled with this parameter, starting value, lower limit, upper limit
	fitter->CreateScalingParameter(tmpVec,0.1375,0.01,1e6);	
	
	//	Define the number of threads the fitter can use
	fitter->SetNthreads(threads);
	// 	Print the starting values
	fitter->Print();					
	//	Set verbocity:
	//	If verbose, fitter will output detailed overview of agreement every 10 steps
	//	If not, fitter will update chi-squared value only	
	fitter->SetVerbose(true);							

	//	MaxIterations/MaxFunctionCalls are for Minuit2 and GSL respectively
	fitter->SetMaxIterations(1000);
	fitter->SetMaxFunctionCalls(1000);
	fitter->SetTolerance(0.0001);

	//	Perform the fit
	fitter->DoFit("Minuit2","Migrad");


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
