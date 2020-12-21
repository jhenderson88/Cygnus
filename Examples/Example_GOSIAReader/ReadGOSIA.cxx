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


void Run(){

	//	Set up the Krypton GOSIA calculation

	const char* krfile = "NucleusFile_Kr80.txt";
	
	// 	Define the nucleus of interest
	NucleusReader *kr_nuclreader = new NucleusReader(krfile);	//	Read from the nucleus file
	Nucleus *kr_nucl = kr_nuclreader->GetNucleus();			//	Construct a nucleus based on the input file
	//kr_nucl->PrintNucleus();
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
	for(size_t i=0;i<kr_i.size();i++)
		kr_bst << kr_nucl->GetMatrixElements().at(kr_l.at(i))[kr_f.at(i)][kr_i.at(i)] << "\n";
	kr_bst.close();

	//	Repeat for the Platinum calculation

	const char* ptfile = "NucleusFile_Pt196.txt";
	
	// 	Define the nucleus of interest
	NucleusReader *pt_nuclreader = new NucleusReader(ptfile);	//	Read from the nucleus file
	Nucleus *pt_nucl = pt_nuclreader->GetNucleus();			//	Construct a nucleus based on the input file
	//pt_nucl->PrintNucleus();
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
	for(size_t i=0;i<pt_i.size();i++)
		pt_bst << pt_nucl->GetMatrixElements().at(pt_l.at(i))[pt_f.at(i)][pt_i.at(i)] << "\n";
	pt_bst.close();

	if(true){	// If we want to simulate some fake yields firs

		system("./gosia < kr80_pt.inp.INTI > /dev/null");	
		system("./gosia < pt196_kr.inp.INTI > /dev/null");
		
		GOSIAReader	kr_gosiaReader_yield(kr_nucl,"kr80_pt.out.INTI");
		GOSIAReader	pt_gosiaReader_yield(pt_nucl,"pt196_kr.out.INTI");
	
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

	system("./gosia < kr80_pt.inp > /dev/null");	
	system("./gosia < kr80_pt.inp.INTI > /dev/null");	

	GOSIAReader	kr_gosiaReader_point(kr_nucl,"kr80_pt.out");
	GOSIAReader	kr_gosiaReader_inti(kr_nucl,"kr80_pt.out.INTI");

	system("cp kr80_pt.out kr80_pt.out.init");

	//std::vector<std::vector<double>>	kr_Corr;	// Krykron point corrections
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

	system("./gosia < pt196_kr.inp > /dev/null");	
	system("./gosia < pt196_kr.inp.INTI > /dev/null");	

	GOSIAReader	pt_gosiaReader_point(pt_nucl,"pt196_kr.out");
	GOSIAReader	pt_gosiaReader_inti(pt_nucl,"pt196_kr.out.INTI");

	system("cp pt196_kr.out pt196_kr.out.init");

	//std::vector<std::vector<double>>	pt_Corr;	// Krypton point corrections
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

	//	Now, proceed to fit the data
	DataReader *dataReader_t = new DataReader(pt_nucl,"SimulatedData_Pt196.txt");
	DataReader *dataReader_b = new DataReader(kr_nucl,"SimulatedData_Kr80.txt");

	//DataReader *dataReader_t = new DataReader(pt_nucl,"PtYields.yld");
	//DataReader *dataReader_b = new DataReader(kr_nucl,"KrYields.yld");

	std::vector<ExperimentData> exptVec_t;
		for(unsigned int e=0;e<dataReader_t->GetExperimentData().size();e++){
		ExperimentData tmpExpt = dataReader_t->GetExperimentData().at(e);
		exptVec_t.push_back(tmpExpt);
	}

	//	Now, proceed to fit the simulated data
	std::vector<ExperimentData> exptVec_b;
	for(unsigned int e=0;e<dataReader_b->GetExperimentData().size();e++){
		ExperimentData tmpExpt = dataReader_b->GetExperimentData().at(e);
		exptVec_b.push_back(tmpExpt);
	}

	std::cout	<< "Begin fitting routines"
			<< std::endl;

	//******************************************************************************//
	//                         BEGIN THE FITTING PROCEDURE                          //
	//******************************************************************************//

	//	Create a simultaneous fitter (similar to CoulExFitter, but allows for beam and target)
	GOSIASimFitter *fitter = new GOSIASimFitter();
	fitter->SetBeamNucleus(kr_nucl);		//	Pass the fitter the beam nucleus
	fitter->SetTargetNucleus(pt_nucl);		//	... and the target nucleus
	for(unsigned int e=0;e<13;e++)
		fitter->AddBeamCorrectionFactor(kr_Corr.at(e));
	for(unsigned int e=0;e<13;e++)
		fitter->AddTargetCorrectionFactor(pt_Corr.at(e));
	fitter->SetBeamData(exptVec_b);			//	Give the fitter the beam data
	fitter->SetTargetData(exptVec_t);		//	... and the target data

	fitter->SetBeamMapping(kr_i,kr_f,kr_l);
	fitter->SetTargetMapping(pt_i,pt_f,pt_l);

	fitter->SetBeamGOSIAInput("kr80_pt.inp");
	fitter->SetTargetGOSIAInput("pt196_kr.inp");

	fitter->SetBeamGOSIAOutput("kr80_pt.out");
	fitter->SetTargetGOSIAOutput("pt196_kr.out");

	fitter->SetBeamBST("kr80_pt.bst");
	fitter->SetTargetBST("pt196_kr.bst");

	fitter->AddBeamMixingRatio(2,1,6,1);
	fitter->AddBeamBranchingRatio(2,1,0,3.096,0.245);

	//	Define literature limits for some of the target data (first excited state lifetime)
	fitter->AddTargetMatrixElement(1,0,1,1.172,0.005);
	fitter->AddTargetMatrixElement(1,1,1,0.83,0.09);
	fitter->AddTargetMatrixElement(1,1,2,1.36,0.03);
	fitter->AddTargetMatrixElement(1,1,3,1.91,0.03);
	fitter->AddTargetMatrixElement(1,2,2,-0.51,0.21);
	fitter->AddTargetMatrixElement(6,1,2,0.376,0.03);

	//	Define literature data for use in constraining the fits
	// 	fitter->AddBeamLifetime(1,0.909,0.245);

	//	Tell the fitter which matrix elements we will be fitting in the beam
	fitter->AddBeamFittingMatrixElement(1,0,1,0.675428,0.1,3);
	fitter->AddBeamFittingMatrixElement(1,0,2,0.0860115,-3,3);
	fitter->AddBeamFittingMatrixElement(1,1,1,-0.648803,-3,3);
	fitter->AddBeamFittingMatrixElement(1,1,2,0.766345,-3,3);
	fitter->AddBeamFittingMatrixElement(1,1,3,0.207623,-3,3);
	fitter->AddBeamFittingMatrixElement(1,1,4,0.932593,-3,3);
	fitter->AddBeamFittingMatrixElement(1,4,4,0.431781,-3,3);
	fitter->AddBeamFittingMatrixElement(6,1,2,0.166004,-3,3);
	//	fitter->AddBeamFittingMatrixElement(1,0,1,0.65,0.4,1.2);
	//	fitter->AddBeamFittingMatrixElement(1,0,2,0.16305,0,2);
	//	fitter->AddBeamFittingMatrixElement(1,1,1,0.1000,-2,2);
	//	fitter->AddBeamFittingMatrixElement(1,1,2,0.49529,-2,2);
	//	fitter->AddBeamFittingMatrixElement(1,1,3,0.27565,0,2);
	//	fitter->AddBeamFittingMatrixElement(1,1,4,0.66796,-2,2);
	//	fitter->AddBeamFittingMatrixElement(1,4,4,-0.3052,-2,2);
	//	fitter->AddBeamFittingMatrixElement(6,1,2,0.01342,-2,2);
	//	... and in the target
	fitter->AddTargetFittingMatrixElement(1,0,1,1.17241,0.1,2);
	fitter->AddTargetFittingMatrixElement(1,1,1,0.83,0,3);
	//fitter->AddTargetFittingMatrixElement(1,1,2,1.35997,-3,3);
	//fitter->AddTargetFittingMatrixElement(1,1,3,1.90993,-3,3);
	//fitter->AddTargetFittingMatrixElement(1,1,5,0.288,-3,3);
	//fitter->AddTargetFittingMatrixElement(1,2,2,-0.51,-3,3);
	//fitter->AddTargetFittingMatrixElement(1,2,5,1.13,-2,3);
	//fitter->AddTargetFittingMatrixElement(6,1,2,0.376,-2,3);
	//	Define a vector of integers to define common scaling
	/*std::vector<int> tmpVec;
	tmpVec.push_back(0);
	tmpVec.push_back(6);
	fitter->CreateScalingParameter(tmpVec);	//	Set the scaling parameters and their allowed range
	tmpVec.push_back(1);	
	tmpVec.push_back(2);
	tmpVec.push_back(3);
	tmpVec.push_back(4);
	tmpVec.push_back(7);	
	tmpVec.push_back(8);
	tmpVec.push_back(9);
	tmpVec.push_back(10);
	fitter->CreateScalingParameter(tmpVec);	//	Set the scaling parameters and their allowed range
	tmpVec.clear();
	tmpVec.push_back(5);
	tmpVec.push_back(11);
	fitter->CreateScalingParameter(tmpVec);	//	Set the scaling parameters and their allowed range
	tmpVec.clear();
	tmpVec.push_back(12);
	fitter->CreateScalingParameter(tmpVec);*/
	std::vector<int> tmpVec;
	for(int i=0;i<13;i++){
		tmpVec.clear();
		tmpVec.push_back(i);
//	}
		fitter->CreateScalingParameter(tmpVec);
	}
	fitter->SetLikelihoodFit(false);

	//	MaxIterations/MaxFunctionCalls are for Minuit2 and GSL respectively
	fitter->SetMaxIterations(50000);
	fitter->SetMaxFunctionCalls(50000);
	fitter->SetVerbose(true);
	//fitter->SetTolerance(0.001);

	//fitter->SetDoFullUncertainty(true);

	//	Perform the fit
	fitter->DoFit("Minuit2","Migrad");
	//fitter->DoFit("Minuit2","Combined");
	//fitter->DoFit("GSLMultiMin", "ConjugateFR");


}


int main(int argc, char** argv){

	Run();

}
