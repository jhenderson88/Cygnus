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
#include "TRandom3.h"
#include <fstream>

#include <iostream>

///
///	The goal of this script is to first simulate a nucleus based on an input file, and then to calculate expected yields
///
///	Then, based on the central (exact) values and using SQRT(N) uncertainties with randomized "experimental" variations
///	to create a new experimental "data" file with realisitic Poisson uncertainties.
///
///	Based on the simulated data, the code will then perform a fit, giving an indication of the sensitivity that can be
///	achieved given certain experimental conditions.
///


void RunFitter(const char* nuclfile = "NucleusFile_Se76.txt", int threads = 1){

	// 	Define the nucleus of interest
	NucleusReader *nuclreader = new NucleusReader(nuclfile);	//	Read from the nucleus file
	Nucleus *nucl = nuclreader->GetNucleus();			//	Construct a nucleus based on the input file
	nucl->PrintNucleus();						//	Print the output to the terminal

	//	Random number generator used to 
	//	create realistic yields
	TRandom3 rand(0);

	//	Define the reaction (beam energy, etc.)	
	Reaction *reac = new Reaction(76,34,208,82,309.7);		//	Construct a Reaction. Arguments: Beam A, Beam Z, Target A, Target Z, Beam energy (MeV)
	reac->SetMass(75.9192,207.9767); 				//	Use the correct nuclear masses. Without this, Cygnus uses the A values to determine kinematics.

	//	Define the experiment.
	//	Note: ExperimentRanges in Cygnus correspond to any individual detector region, like a subset of rings on an annular silicon for which yields will be provided	
	Experiments *expts = new Experiments(nucl,reac);		//	Construct the Experiment with the previously defined Nucleus and Reaction.
	expts->FixStep(true);
	expts->SetUseSymmetry(true);
	//expts->SetAccuracy(1e-6);					//	Set the accuracy of the calculation. 1e-5 is rather low resolution. Default is 1e-8. Reduced accuracy leads to increased speed.
	double tmin[6] = {20,30,40,20,30,40};				//	Minimum theta (lab, degrees) of each of the six defined experimental ranges
	double tmax[6] = {30,40,50,30,40,50};				//	Maximum theta (lab, degrees) of each of the six defined experimental ranges
	bool tarDet[6] = {false,false,false,true,true,true};		//	Is the target detected in the particle detector? Defines the kinematics of the reaction.
	//	Create the ExperimentRange objects which are stored within the Experiment
	// 	Arguments: minimum theta (lab), maximum theta (lab), number of theta meshpoints, minimum energy (MeV), maximum energy (MeV), target detection bool
	for(int i=0;i<6;i++)
		expts->NewExperimentRange(tmin[i],tmax[i],6,289.7,309.7,5,tarDet[i]);

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

	//	Write the integration distributions to files
	expts->WriteIntegralFits("OutputFitting_Simulated.root","RECREATE");	

	//	Constants for turning calculated yields into numbers of events		
	double 	YieldFactor = 0.00000756746 / 208.;	//	6.022E23 * 1E-30 * 4pi (NA * barns * solid angle)
	double 	runningSeconds = 3 * 24 * 60 * 60;	//	Length of the experiment in seconds
	double 	beamIntensity = 1e5;			//	Beam intensity in pps
	double	detectionEff = 0.1;			//	Detection efficiency, for simplicity this is assumed to be constant

	//	Output text file where we will write our simulated yields
	std::ofstream outfile("SimulatedData.txt");

	//	Syntax for the yield file
	outfile << "! Syntax: EXPT (flag for reader), nExpt, beam energy (MeV), theta_min (deg), theta_max (deg)"
		<< "\n";

	//	Object to convert cross sections to yields (using branching ratios)
	TransitionRates *rates = new TransitionRates(nucl);
	rates->Print();
	std::cout	<< std::endl;

	//	Prepare the simulated yields
	for(unsigned int e=0;e<6;e++){
		//	Syntax for each experiment for the yield file
		outfile	<< "EXPT\t" 
			<< e+1 	<< "\t"
			<< reac->GetLabEnergy() << "\t"
			<< tmin[e] << "\t"
			<< tmax[e] << "\n";
		//	Print the yields as determined from the TransitionRates object and the experimental cross
		//	section. Note that this is not yet in counting units.
		GammaYield::PrintYields(expts->GetExperimentRange(e),*rates,*nucl);
		//	Create a matrix containing the experimental "yields" (again, not yet in counting units)
		TMatrixD tmpMat = GammaYield::GammaRayYield(expts->GetExperimentRange(e),*rates);
		//	Loop over every state...
		for(int f = 0; f<nucl->GetNstates(); f++){
			//	... and evert other state
			for(int i = f+1; i<nucl->GetNstates(); i++){
				//	Calculate the "true" yield, now in correct units (events/experiment)
				double tmpYld = tmpMat[f][i] * YieldFactor * runningSeconds * beamIntensity * detectionEff;
				if(tmpYld > 0){
					//	Determine a simuated yield, based on random sampling of a Poisson distributions
					double tmpYldPoiss = rand.PoissonD(tmpYld);
					//	... and determine an uncertainty (for simplicity: sqrt(N))
					double tmpErr = TMath::Sqrt(tmpYldPoiss);
					//	Write the yield and uncertainty to file
					outfile	<< i << "\t"
						<< f << "\t"
						<< tmpYldPoiss << "\t"
						<< tmpErr << "\n";	
				}
			}
		}
	}
	outfile << "END\n";
	outfile.close();
	std::cout	<< std::endl;
	//	Yield simulation is now finished, and the yields are written to file

	//	Print the point corrections to the screen
	expts->PrintPointCorrections();

	//	Now, proceed to fit the simulated data
	DataReader *dataReader = new DataReader(nucl,"SimulatedData.txt");
	std::vector<ExperimentData> exptVec;
	std::vector<PointCoulEx> poinVec;
	for(unsigned int e=0;e<6;e++){
		ExperimentData tmpExpt = dataReader->GetExperimentData().at(e);
		tmpExpt.SetThetaCM(expts->GetExperimentRange(e).GetMeanThetaCM());
		exptVec.push_back(tmpExpt);
		PointCoulEx tmpPoin = expts->GetPointCalculation(e);
		poinVec.push_back(tmpPoin);
	}

	std::cout	<< "Begin fitting routines"
			<< std::endl;

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

	//	Define the matrix elements to be varied in the fit. Any not defined here will be fixed to the values in the Nucleus.
	fitter->AddFittingMatrixElement(1,0,1,0.61,0.55,0.70);
	fitter->AddFittingMatrixElement(1,0,3,0.271,0.05,0.4);
	fitter->AddFittingMatrixElement(1,1,1,-0.200,-0.8,-0.1);
	fitter->AddFittingMatrixElement(1,1,2,0.6999,0.2,1);
	fitter->AddFittingMatrixElement(1,1,3,0.412,0.1,1.2);
	fitter->AddFittingMatrixElement(1,1,4,1.1,0.8,1.2);

	//	This section tells the fitter to use common scaling for all ExperimentRanges
	std::vector<int> tmpVec;
	for(unsigned int e=0;e<6;e++)
		tmpVec.push_back((int)e);
	// 	Arguments: vector containing experiments to be scaled with this parameter, starting value, lower limit, upper limit
	fitter->CreateScalingParameter(tmpVec,0.001,1e-6,1e6);	
	
	//	Define the number of threads the fitter can use
	fitter->SetNthreads(threads);
	// 	Print the starting values
	fitter->Print();					
	//	Set verbocity:
	//	If verbose, fitter will output detailed overview of agreement every 10 steps
	//	If not, fitter will update chi-squared value only	
	fitter->SetVerbose(false);							
	fitter->SetLikelihoodFit();

	//	MaxIterations/MaxFunctionCalls are for Minuit2 and GSL respectively
	fitter->SetMaxIterations(5000);
	fitter->SetMaxFunctionCalls(5000);
	fitter->SetTolerance(0.01);

	//	Perform the fit
	fitter->DoFit("Minuit2","Combined");

	//******************************************************************************//
	//                     REPEAT INTEGRAL AND FIT FINAL YIELDS                     //
	//******************************************************************************//

	std::vector<double> fittedPar = fitter->GetFitParameters();
	Nucleus *nuclFitted = nuclreader->GetNucleus();
	for(size_t me = 0; me < fitter->GetMatrixElements().size(); me++){
		MatrixElement matrixElement = fitter->GetMatrixElements().at(me);
		nuclFitted->SetMatrixElement(matrixElement.GetLambda(),
						matrixElement.GetInitialState(),		
						matrixElement.GetFinalState(),
						fittedPar[me]);
	}
	for(size_t i = 0; i < fittedPar.size(); i++)
		std::cout	<< std::setw(12) << std::left << fittedPar[i];
	std::cout	<< std::endl;
	Experiments *exptsFitted = new Experiments(nuclFitted,reac);	
	exptsFitted->FixStep(true);
	exptsFitted->SetUseSymmetry(true);
	for(int i=0;i<6;i++)
		exptsFitted->NewExperimentRange(tmin[i],tmax[i],11,289.7,309.7,5,tarDet[i]);

	exptsFitted->SetStopping(dEdX);		//	Add the stopping powers to the experiment
	exptsFitted->SetVerbose(false);		//	Set verbocity
	exptsFitted->SetNthreads(threads);	//	Define the number of threads to be used in the integration process

	exptsFitted->PointCorrections();
	std::ofstream fittedYields;
	fittedYields.open("FittedYields.txt");

	for(unsigned int e=0;e<6;e++)
		GammaYield::WriteYields(exptsFitted->GetExperimentRange(e),*rates,*nuclFitted,fittedYields,fittedPar.at(fittedPar.size()-1));

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
