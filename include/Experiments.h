#ifndef Experiments_h
#define Experiments_h

#include "ParticleDetector.h"
#include "StoppingPower.h"
#include "PointCoulEx.h"
#include "ExperimentRange.h"
#include "Nucleus.h"
#include "Reaction.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"

#include <vector>
#include <iostream>

class PointCoulEx;
class ExperimentRange;
class Nucleus;
class Reaction;

///
///	\class Experiments
///
///	\brief The Experiments class is the umbrella class into 
///	which all defined experiments are stored.
///
///	Each "Experiment" does not need to be a different
///	physical experiment, but might be a different theta
///	range on a detector, for example.
///
///	The Experiment class then contains a vector of 
///	ExperimentRange objects, in which theta-enegry
///	integrals are evaluated to determined an integrated
///	cross section. PointCoulEx calculations are also
///	performed, and a correction factor can thereby be
///	determined: effectively correcting point calcula-
///	tions to reproduce integrated yields. This varies
///	state-by-state, depending upon the angular distr-
///	ibution of the cross section, and is therefore
///	reported accordingly.
///

class Experiments
{

	public:
		Experiments();
		Experiments(Nucleus*,Reaction*);		/*!< Constructor passing a Nucleus and Reaction kinematics to the object */
		Experiments(const Experiments&);		/*!< Copy constructor */
		Experiments& operator = (const Experiments&);	/*!< Assignment operator */
		~Experiments() {;}

		///
		///	Define a new experimental range (for which data will be provided,
		///	or yields are required). 
		///
		///	thetaMin and thetaMax define the range in theta over which the 
		///	integration will be performed, with nTheta defining the number 
		///	of theta meshpoints to be used.
		///
		///	energyMin and energyMax define the energy range over which the
		///	integration will be performed, with nEnergy defining the number
		///	of energy meshpoints to be used.
		///
		///	targetDetection (default = false) determines whether the beam-
		///	or target-like nucleus was detected, allowing for the correct
		///	kinematic conversion.
		///
		void	NewExperimentRange(double thetaMin, double thetaMax, int nTheta, double energyMin, double energyMax, int nEnergy,bool targetDetection = false);		// Thetamin, thetamax, nTheta, energymin, energymax, nEnergy

		void	ClearExperiments()	{ experimentRanges.clear();	}	/*!<	Delete all previously defined experimental ranges */

		void	PointCorrections();						/*!<	Determine the correction factor required to convert a point to an energy and angle integrated calculation */

		void	PrintDetails() const;						/*!<	Print the details of the experimental ranges stored within the class */
		void	PrintPointCorrections();					/*!<	Print the correction factors as determined by PointCorrections */

		ExperimentRange			GetExperimentRange(int i)	const	{ return experimentRanges.at((unsigned int)i); 			}	/*!< Returns an experimental range by index */
		std::vector<TVectorD>		GetCorrectionFactors()		const	{ return correctionFactors;					}	/*!< Returns the vector of correction factors (TVectorD) */
		std::vector<TVectorD>		GetPointCrossSections()		const	{ return pointCrossSections;					}	/*!< Returns the vector of point cross sections (TVectorD) */
		std::vector<TVectorD>		GetIntegratedCrossSections()	const	{ return integratedCrossSections;				}	/*!< Returns the vector of integrated cross sections (TVectorD)*/

		PointCoulEx			GetPointCalculation(int i)	const	{ return pointCalculation.at((unsigned int)i);			}	/*!< Returns a point calculation by index */
		std::vector<PointCoulEx>	GetPointCalculations()		const	{ return pointCalculation;					}	/*!< Returns the vector of point calculations */

		int				GetNexpts()			const	{ return experimentRanges.size();				}	/*!< Returns the number of experimental ranges defined */

		void	SetVerbose(bool b = true)	{ verbose = b;		}										/*!< Set the verbocity of the calculation (for debugging) */

		void	WriteIntegralFits(const char*, const char* opt = "UPDATE");										/*!< Writes the theta-energy cross section surfaces (TGraph & TGraph2D) to a ROOT file */

		void	SetNthreads(int n)						{ nThreads = n;							}	/*!< Set the number of cores the code can use for integration */
		int	GetNthreads()						const	{ return nThreads;						}	/*!< Returns the number of cores the code is allowed to use*/

		void	SetAccuracy(double acc)						{ fAccuracy = acc;						}	/*!< Define the accuracy of the CoulEx calculations to be performed */

		void	SetStopping(StoppingPower s)					{ fStopping = s;						}	/*!< Define the stopping powers to be used in the integration process */

		void	SetParticleDetectorEff(int i, TGraph* p)			{ experimentRanges.at((unsigned int)i).SetDetectorEff(p);	}	/*!< Define the detection efficiency of the particle detector vs theta (lab) */

		void	UseEfficiency(bool b = true)					{ fUseEfficiency = b;						}	/*!< Use the efficiencies of the particle detector */

		void	SetProjectileExcitation(bool b = true)				{ fProjectileExcitation = b;					}	/*!< Determines whether the beam or target is being excited */

	private:

		double		fAccuracy;

		Nucleus		fNucleus;
		Reaction	fReaction;

		StoppingPower	fStopping;

		std::vector<ExperimentRange>	experimentRanges;
		std::vector<PointCoulEx>	pointCalculation;

		std::vector<TVectorD>		integratedCrossSections;
		std::vector<TVectorD>		pointCrossSections;
		std::vector<TVectorD>		correctionFactors;

		bool		verbose;
		bool		fUseEfficiency;
		bool		fProjectileExcitation;

		int 		nThreads;


};
#endif
