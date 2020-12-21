#include "TMatrixD.h"
#include "TVectorD.h"

#include <stdio.h>
#include <stdlib.h>
#include "ScalingParameter.h"
#include "CoulExFitter.h"
#include "GammaYield.h"
#include "PointCoulEx.h"
#include "Experiments.h"
#include "TransitionRates.h"
#include "GOSIAReader.h"
#include "DataReader.h"
#include "Nucleus.h"
#include "NucleusReader.h"
#include "StoppingPower.h"
#include "TRandom3.h"
#include <fstream>

#include "GOSIASimFitter.h"

#include <iostream>

///
///	The goal of this script is to use the traditional GOSIA code to perform the Coulomb
///	excitation analysis, but with a Minuit-based minimization algorithm.
///
///	In order for the minimizer to work, GOSIA needs to print its calculated yields with much
///	higher precision than is standard. This can be done by modifying the GOSIA write statement 99049 to:  
///	"99049 FORMAT (4X,1I3,4X,1I3,3X,1F4.1,3X,1F4.1,3X,F60.50,8X,1E11.5)"
///
///	This will make the formatting look pretty hideous, but gives the minimizer access to all of the precision
///	it needs.
///
///	The example below has comments describing the process. This code is a GOSIA in, GOSIA out
///	example. It uses GOSIA to calculate some yields, applies a typical poisson smearing to the data,
///	and then fits it. So everything should be ideal.
///

void Run(){

	//	Set up the Krypton GOSIA calculation
	const char* krfile = "NucleusFile_Kr80.txt";
	
	// 	Define the nucleus of interest
	NucleusReader *kr_nuclreader = new NucleusReader(krfile);	//	Read from the nucleus file
	Nucleus *kr_nucl = kr_nuclreader->GetNucleus();			//	Construct a nucleus based on the input file

	///
	///	We use the GOSIA OP,REST option to reset the matrix elements easily
	///	This requires us to know the ordering of the matrix elements, as provided in the 
	///	GOSIA file. The vectors below are used to specify this ordering so that the correct
	///	matrix elements are modified.
	///
	std::ofstream		kr_bst("kr80_pt.bst");
	std::vector<int>	kr_i, kr_f, kr_l;	// Initial, final, lambda
	kr_i.resize(8); kr_f.resize(8); kr_l.resize(8);
	kr_i[0]	= 1;	kr_f[0] = 0; kr_l[0] = 1;
	kr_i[1]	= 2;	kr_f[1] = 0; kr_l[1] = 1;
	kr_i[2]	= 1;	kr_f[2] = 1; kr_l[2] = 1;
	kr_i[3]	= 2;	kr_f[3] = 1; kr_l[3] = 1;
	kr_i[4]	= 3;	kr_f[4] = 1; kr_l[4] = 1;
	kr_i[5]	= 4;	kr_f[5] = 1; kr_l[5] = 1;
	kr_i[6]	= 4;	kr_f[6] = 4; kr_l[6] = 1;
	kr_i[7]	= 2;	kr_f[7] = 1; kr_l[7] = 6;

	///	Write the initial matrix elements to the .bst file, as specified in the Nucleus file and
	///	with an ordering defined by the above vectors
	for(size_t i=0;i<kr_i.size();i++)
		kr_bst << kr_nucl->GetMatrixElements().at(kr_l.at(i))[kr_f.at(i)][kr_i.at(i)] << "\n";
	kr_bst.close();

	//	Repeat for the Platinum calculation

	const char* ptfile = "NucleusFile_Pt196.txt";
	
	// 	Define the nucleus of interest
	NucleusReader *pt_nuclreader = new NucleusReader(ptfile);	//	Read from the nucleus file
	Nucleus *pt_nucl = pt_nuclreader->GetNucleus();			//	Construct a nucleus based on the input file

	///
	///	We use the GOSIA OP,REST option to reset the matrix elements easily
	///	This requires us to know the ordering of the matrix elements, as provided in the 
	///	GOSIA file. The vectors below are used to specify this ordering so that the correct
	///	matrix elements are modified.
	///
	std::ofstream		pt_bst("pt196_kr.bst");
	std::vector<int>	pt_i, pt_f, pt_l;	// Initial, final, lambda
	pt_i.resize(8); pt_f.resize(8); pt_l.resize(8);
	pt_i[0]	= 1;	pt_f[0] = 0; pt_l[0] = 1;
	pt_i[1]	= 1;	pt_f[1] = 1; pt_l[1] = 1;
	pt_i[2]	= 2;	pt_f[2] = 1; pt_l[2] = 1;
	pt_i[3]	= 3;	pt_f[3] = 1; pt_l[3] = 1;
	pt_i[4]	= 5;	pt_f[4] = 1; pt_l[4] = 1;
	pt_i[5]	= 2;	pt_f[5] = 2; pt_l[5] = 1;
	pt_i[6]	= 5;	pt_f[6] = 2; pt_l[6] = 1;
	pt_i[7]	= 2;	pt_f[7] = 1; pt_l[7] = 6;

	///	Write the initial matrix elements to the .bst file, as specified in the Nucleus file and
	///	with an ordering defined by the above vectors
	for(size_t i=0;i<pt_i.size();i++)
		pt_bst << pt_nucl->GetMatrixElements().at(pt_l.at(i))[pt_f.at(i)][pt_i.at(i)] << "\n";
	pt_bst.close();

	///
	///	If we want to fake up some yields using GOSIA first.
	///	Otherwise, we can just format existing yields correctly and feed them
	///	directly into the minimizer.
	///
	if(true){	

		system("./gosia < kr80_pt.inp.INTI > /dev/null");	// 	Run the beam-like GOSIA integration	
		system("./gosia < pt196_kr.inp.INTI > /dev/null");	//	Run the target-like GOSIA integration
		
		GOSIAReader	kr_gosiaReader_yield(kr_nucl,"kr80_pt.out.INTI");	//	Read in the results
		GOSIAReader	pt_gosiaReader_yield(pt_nucl,"pt196_kr.out.INTI");	//	Read in the results
	
		//	Here we define the "experiment":
		double 	YieldFactor = 0.00000756746 / 196.;	//	6.022E23 * 1E-30 * 4pi (NA * barns * solid angle)
		double 	runningSeconds = 6 * 60 * 60;		//	Length of the experiment in seconds
		double 	beamIntensity = 1e5;			//	Beam intensity in pps
		double	detectionEff = 0.1;			//	Detection efficiency, for simplicity this is assumed to be constant

		TRandom3 rand;

		//	Output text file where we will write our simulated yields
		std::ofstream outfile_kr("SimulatedData_Kr80.txt");

		//	Syntax for the yield file
		outfile_kr 	<< "! Syntax: EXPT (flag for reader), nExpt, beam energy (MeV), theta_min (deg), theta_max (deg)"
				<< "\n";

		//	Prepare the simulated yields
		for(size_t e=0;e<kr_gosiaReader_yield.GetGOSIAData().size();e++){
			//	Syntax for each experiment for the yield file
			outfile_kr	<< std::setw(12) << std::left << "EXPT" 
					<< std::setw(12) << std::left<< e+1 
					<< std::setw(12) << std::left<< 320
					<< std::setw(12) << std::left<< 0
					<< std::setw(12) << std::left<< 0 
					<< "\n";
			for(size_t ee=0;ee<kr_gosiaReader_yield.GetGOSIAData().at(e).GetData().size();ee++){
				//	Calculate the "true" yield, now in correct units (events/experiment)
				int init = kr_gosiaReader_yield.GetGOSIAData().at(e).GetDataPoint(ee).GetInitialIndex();
				int fina = kr_gosiaReader_yield.GetGOSIAData().at(e).GetDataPoint(ee).GetFinalIndex();
				double tmpYld = kr_gosiaReader_yield.GetGOSIAData().at(e).GetDataPoint(ee).GetCounts();
				tmpYld *= YieldFactor;
				tmpYld *= runningSeconds;
				tmpYld *= beamIntensity;
				tmpYld *= detectionEff;
				if(tmpYld > 0){
					//	Determine a simuated yield, based on random sampling of a Poisson distributions
					double tmpYldPoiss = rand.PoissonD(tmpYld);
					//	... and determine an uncertainty (for simplicity: sqrt(N))
					double tmpErr = TMath::Sqrt(tmpYldPoiss);
					//	Write the yield and uncertainty to file
					outfile_kr	<< std::setw(12) << std::left<< init
							<< std::setw(12) << std::left<< fina
							<< std::setw(12) << std::left<< tmpYldPoiss / detectionEff 
							<< std::setw(12) << std::left<< tmpErr / detectionEff 
							<< "\n";	
				}
			}
		}
		outfile_kr << "END\n";
		outfile_kr.close();

		detectionEff = 0.13;			//	Detection efficiency, for simplicity this is assumed to be constant

		//	Output text file where we will write our simulated yields
		std::ofstream outfile_pt("SimulatedData_Pt196.txt");

		//	Syntax for the yield file
		outfile_pt 	<< "! Syntax: EXPT (flag for reader), nExpt, beam energy (MeV), theta_min (deg), theta_max (deg)"
				<< "\n";

		//	Prepare the simulated yields
		for(size_t e=0;e<pt_gosiaReader_yield.GetGOSIAData().size();e++){
			//	Syntax for each experiment for the yield file
			outfile_pt	<< std::setw(12) << std::left<< "EXPT" 
					<< std::setw(12) << std::left<< e+1
					<< std::setw(12) << std::left<< 320
					<< std::setw(12) << std::left<< 0
					<< std::setw(12) << std::left<< 0 
					<< "\n";
			for(size_t ee=0;ee<pt_gosiaReader_yield.GetGOSIAData().at(e).GetData().size();ee++){
				//	Calculate the "true" yield, now in correct units (events/experiment)
				int init = pt_gosiaReader_yield.GetGOSIAData().at(e).GetDataPoint(ee).GetInitialIndex();
				int fina = pt_gosiaReader_yield.GetGOSIAData().at(e).GetDataPoint(ee).GetFinalIndex();
				double tmpYld = pt_gosiaReader_yield.GetGOSIAData().at(e).GetDataPoint(ee).GetCounts();
				tmpYld *= YieldFactor;
				tmpYld *= runningSeconds;
				tmpYld *= beamIntensity;
				tmpYld *= detectionEff;
				if(tmpYld > 0){
					//	Determine a simuated yield, based on random sampling of a Poisson distributions
					double tmpYldPoiss = rand.PoissonD(tmpYld);
					//	... and determine an uncertainty (for simplicity: sqrt(N))
					double tmpErr = TMath::Sqrt(tmpYldPoiss);
					//	Write the yield and uncertainty to file
					outfile_pt	<< std::setw(12) << std::left<< init
							<< std::setw(12) << std::left<< fina
							<< std::setw(12) << std::left<< tmpYldPoiss / detectionEff
							<< std::setw(12) << std::left<< tmpErr / detectionEff 
							<< "\n";	
				}
			}
		}
		outfile_pt << "END\n";
		outfile_pt.close();

	}

	///
	///	Now we have our yields (simulated or real) we need to calculate the point
	///	to integrated correction factors
	///


	///
	///	First, the beam
	///

	system("./gosia < kr80_pt.inp > /dev/null");		//	Run the point, beam-like calculations
	system("./gosia < kr80_pt.inp.INTI > /dev/null");	//	Run the integrated, beam-like calculations
	
	GOSIAReader	kr_gosiaReader_point(kr_nucl,"kr80_pt.out");
	GOSIAReader	kr_gosiaReader_inti(kr_nucl,"kr80_pt.out.INTI");

	//	Copy the output, for later comparison
	system("cp kr80_pt.out kr80_pt.out.init");

	///
	///	Here we loop over the yields and calculate the correction factor, before passing it
	///	to a vector to be used in the minimization
	///
	std::vector<TVectorD>	kr_Corr;
	for(size_t e=0; e < kr_gosiaReader_point.GetGOSIAData().size(); e++){	// Loop over experiments
		size_t	len = kr_gosiaReader_point.GetGOSIAData().at(e).GetData().size();
		TVectorD	tmpVec;
		tmpVec.ResizeTo(len);	
		for(size_t ee=0; ee<kr_gosiaReader_point.GetGOSIAData().at(e).GetData().size();ee++){
			double	point = kr_gosiaReader_point.GetGOSIAData().at(e).GetDataPoint(ee).GetCounts();
			double	inti = kr_gosiaReader_inti.GetGOSIAData().at(e).GetDataPoint(ee).GetCounts();
			if(point > 1e-8)
				tmpVec[ee] = inti/point;
			else
				tmpVec[ee] = 0;
		}
		kr_Corr.push_back(tmpVec);
	}

	///
	///	Repeat the process for the target
	///

	system("./gosia < pt196_kr.inp > /dev/null");		//	Run the point, target-like calculations
	system("./gosia < pt196_kr.inp.INTI > /dev/null");	//	Run the integrated, beam-like calculations

	GOSIAReader	pt_gosiaReader_point(pt_nucl,"pt196_kr.out");
	GOSIAReader	pt_gosiaReader_inti(pt_nucl,"pt196_kr.out.INTI");

	//	Copy the output, for later comparison
	system("cp pt196_kr.out pt196_kr.out.init");

	///
	///	Here we loop over the yields and calculate the correction factor, before passing it
	///	to a vector to be used in the minimization
	///
	std::vector<TVectorD>	pt_Corr;
	for(size_t e=0; e < pt_gosiaReader_point.GetGOSIAData().size(); e++){	// Loop over experiments
		size_t	len = pt_gosiaReader_point.GetGOSIAData().at(e).GetData().size();
		TVectorD	tmpVec;
		tmpVec.ResizeTo(len);	
		for(size_t ee=0; ee<pt_gosiaReader_point.GetGOSIAData().at(e).GetData().size();ee++){
			double	point = pt_gosiaReader_point.GetGOSIAData().at(e).GetDataPoint(ee).GetCounts();
			double	inti = pt_gosiaReader_inti.GetGOSIAData().at(e).GetDataPoint(ee).GetCounts();
			if(point > 1e-8)
				tmpVec[ee] = inti/point;
			else
				tmpVec[ee] = 0;
		}
		pt_Corr.push_back(tmpVec);
	}

	///
	///	Now we read in our experimental data
	///
	///	If you want to use other data than simulated, here's where you include it
	///
	DataReader *dataReader_t = new DataReader(pt_nucl,"SimulatedData_Pt196.txt");
	DataReader *dataReader_b = new DataReader(kr_nucl,"SimulatedData_Kr80.txt");

	///	
	///	First read in the target data
	///
	std::vector<ExperimentData> exptVec_t;
	for(unsigned int e=0;e<dataReader_t->GetExperimentData().size();e++){
		ExperimentData tmpExpt = dataReader_t->GetExperimentData().at(e);
		exptVec_t.push_back(tmpExpt);
	}

	///
	///	Then the beam
	///
	std::vector<ExperimentData> exptVec_b;
	for(unsigned int e=0;e<dataReader_b->GetExperimentData().size();e++){
		ExperimentData tmpExpt = dataReader_b->GetExperimentData().at(e);
		exptVec_b.push_back(tmpExpt);
	}

	///
	///	Now we have everything we need to begin the fitting procedure
	///

	std::cout	<< "Begin fitting routines"
			<< std::endl;

	//******************************************************************************//
	//                         BEGIN THE FITTING PROCEDURE                          //
	//******************************************************************************//

	///	Create a simultaneous fitter, designed to use external GOSIA yields
	GOSIASimFitter *fitter = new GOSIASimFitter();
	fitter->SetBeamNucleus(kr_nucl);		//	Pass the fitter the beam nucleus...
	fitter->SetTargetNucleus(pt_nucl);		//	... and the target nucleus
	///
	///	Now we need to tell the fitter about the correction factors for the beam...
	///
	for(unsigned int e=0;e<13;e++)
		fitter->AddBeamCorrectionFactor(kr_Corr.at(e));
	///
	///	... and the target
	///
	for(unsigned int e=0;e<13;e++)
		fitter->AddTargetCorrectionFactor(pt_Corr.at(e));
	fitter->SetBeamData(exptVec_b);			//	Give the fitter the beam data...
	fitter->SetTargetData(exptVec_t);		//	... and the target data

	///
	///	We also need to tell the fitter about the matrix-element ordering in the GOSIA
	///	file. So we pass the vector objects here.
	///
	fitter->SetBeamMapping(kr_i,kr_f,kr_l);
	fitter->SetTargetMapping(pt_i,pt_f,pt_l);

	///
	///	... and the fitter needs to know the various filenames to run, read and modify with
	///	new matrix elements.
	///
	fitter->SetBeamGOSIAInput("kr80_pt.inp");	// 	Input file
	fitter->SetTargetGOSIAInput("pt196_kr.inp");	// 	Input file
	fitter->SetBeamGOSIAOutput("kr80_pt.out");	//	Output file
	fitter->SetTargetGOSIAOutput("pt196_kr.out");	// 	Output file
	fitter->SetBeamBST("kr80_pt.bst");		//	Matrix element file
	fitter->SetTargetBST("pt196_kr.bst");		//	Matrix element file

	///
	///	Give the fitter some literature data to play with
	///
	fitter->AddTargetMatrixElement(1,0,1,1.172,0.005);
	fitter->AddTargetMatrixElement(1,1,1,0.83,0.09);

	///
	///	Tell the fitter which matrix elements we're fitting in the beam...
	///
	fitter->AddBeamFittingMatrixElement(1,0,1,0.675428,0.1,3);
	fitter->AddBeamFittingMatrixElement(1,0,2,0.0860115,-3,3);
	fitter->AddBeamFittingMatrixElement(1,1,1,-0.648803,-3,3);
	fitter->AddBeamFittingMatrixElement(1,1,2,0.766345,-3,3);
	fitter->AddBeamFittingMatrixElement(1,1,3,0.207623,-3,3);
	fitter->AddBeamFittingMatrixElement(1,1,4,0.932593,-3,3);
	fitter->AddBeamFittingMatrixElement(1,4,4,0.431781,-3,3);
	fitter->AddBeamFittingMatrixElement(6,1,2,0.166004,-3,3);
	///
	///	... and in the target
	///
	fitter->AddTargetFittingMatrixElement(1,0,1,1.17241,0.1,2);
	fitter->AddTargetFittingMatrixElement(1,1,1,0.83,0,3);

	///	
	///	We're going to scale each angular range independently in this case. So we create
	///	N vectors, containing a single integer n, where N is the number of angular ranges
	///	(experiments, in the GOSIA nomenclature) and n is the number of the experiment.
	///
	std::vector<int> tmpVec;
	for(int i=0;i<13;i++){
		tmpVec.clear();
		tmpVec.push_back(i);
		fitter->CreateScalingParameter(tmpVec);
	}
	fitter->SetLikelihoodFit(false);

	///
	///	Pass some fitter formatting information.
	///
	///	MaxIterations/MaxFunctionCalls are for Minuit2 and GSL respectively.
	///
	fitter->SetMaxIterations(50000);
	fitter->SetMaxFunctionCalls(50000);
	fitter->SetVerbose(true);

	///
	///	... and let's do a full uncertainty analysis, allowing for assymmetric errors
	///
	fitter->SetDoFullUncertainty(true);

	///
	///	Nowe we're all set up
	///
	///	Perform the fit
	///
	fitter->DoFit("Minuit2","Migrad");


}


int main(int argc, char** argv){

	Run();

}
