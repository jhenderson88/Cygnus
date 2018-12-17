#ifndef ExperimentRange_h
#define ExperimentRange_h

#include "ParticleDetector.h"
#include "ParticleDetectorS3.h"
#include "StoppingPower.h"
#include "Nucleus.h"
#include "Reaction.h"
#include "Integral.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1F.h"
#include "TFile.h"
#include <chrono>
#include <algorithm>

class Nucleus;
class Reaction;
class Integral;

///
///	\class ExperimentRange
///
///	\brief Define an experimental range, for determining integrated
///	yields.
///
///	The ExperimentRange class uses Integral to calculate the cross-
///	section meshpoints, state-by-state. It then fits the theta-energy 
///	cross-section dependence and integrates the surface, along with
///	a correction for energy loss in the target, in order to calculate 
///	an integrated cross section.
///

class ExperimentRange {

	public:
		ExperimentRange();
		ExperimentRange(Nucleus*,Reaction*);				/*!< Constructor with a Nucleus and Reaction class */
		ExperimentRange(const ExperimentRange &e);			/*!< Copy constructor */
		ExperimentRange& operator = (const ExperimentRange& e);		/*!< Assignment operator */
		~ExperimentRange() {;}
	
		///
		///	Define an experimental range:\n
		///	tmin - minimum theta (lab)
		///	tmax - maximum theta (lab)
		///	nt   - number of theta meshpoints for integration
		///	emin - minimum energy (MeV)
		///	emax - maximum energy (MeV)
		///	ne   - number of energy meshpoints for integration
		///
		void SetRanges(double tmin, double tmax, int nt, double emin, double emax, int ne) {
			SetThetaMin(tmin); 
			SetThetaMax(tmax); 
			SetNtheta(nt); 
			SetEnergyMin(emin); 
			SetEnergyMax(emax); 
			SetNenergy(ne);
			SetMeanValues();
		}

		void		SetMeanValues();	/*!< Define the mean theta and mean energy, for use in the point calculations */

		void		SetThetaMin(double t)	{ thetamin = t;	}	/*!< Define the minimum theta for integration (lab) */
		void		SetThetaMax(double t)	{ thetamax = t;	}	/*!< Define the maximum theta for integration (lab) */

		void		SetEnergyMin(double e)	{ energymin = e; }	/*!< Define the minimum energy for integration (MeV, lab) */
		void		SetEnergyMax(double e)	{ energymax = e; }	/*!< Define the maximum energy for integration (MeV, lab) */

		void		SetNtheta(int n)	{ nTheta = n;	}	/*!< Define the number of theta meshpoints for integration */
		void		SetNenergy(int n)	{ nEnergy = n;	}	/*!< Define the number of energy meshpoints for integration */

		void		IntegrateRange();				/*!< Perform the integration over the theta and energy ranges */
		TGraph2D*	GetThetaEnergyGraph(int s);			/*!< Return the cross-section graph vs theta and energy for a given state, s */
		TGraph2D*	InterpolatedEnergyTheta(int,bool useDetector = false);		/*!< Perform an interpolation between the meshpoints determined in the integration */

		void		SetTargetDetection(bool b)	{ targetDetection = b;	}	/*!< Define whether the target is being detected (default: beam-like detection */

		TVectorD	GetIntegratedCrossSection_TVec()		const	{ return IntegratedCrossSection_TVec;	}	/*!< Return a TVectorD containing the integrated cross section over theta and energy */

		double		IntegrateRutherford();				/*!< Integrate the Rutherford cross section over the angle and energy range */
		TGraph2D*	GetRutherfordThetaEnergy();			/*!< Create a TGraph2D of the Rutherford cross section over theta (cm) and energy */

		double		GetMeanThetaLab()				const	{ return meanThetaLab;			}	/*!< Return the mean theta value (lab) for point calculation  */
		double		GetMeanThetaCM()				const	{ return meanThetaCM;			}	/*!< Return the mean theta value (CoM) for point calculation  */
		double		GetMeanEnergy()					const	{ return meanEnergy;			}	/*!< Return the mean energy (MeV, lab) for point calculation */
	
		double		GetThetaMin() 					const	{ return thetamin;			}	/*!< Return the minimum theta value (lab) for integration */
		double		GetThetaMax()					const	{ return thetamax;			}	/*!< Return the maximum theta value (lab) for integration */
		double		GetEMin()					const	{ return energymin;			}	/*!< Return the minimum energy value (MeV, lab) for integration */
		double		GetEMax()					const	{ return energymax;			}	/*!< Return the maximum energy valie (MeV, lab) for integration */

		void		PrintDetails() const;	/*!< Print the details of the experimental range */
	
		void		SetVerbose(bool b = true)				{ verbose = b;				} 	/*!< Define the verbocity of the calculations */

		void		SetNthreads(int n)					{ nThreads = n;				}	/*!< Define the number of cores to be used in the integration */	
		int		GetNthreads()					const	{ return nThreads;			}	/*!< Return the number of cores to be used in the integration */

		void		SetAccuracy(double acc)					{ fAccuracy = acc;			}	/*!< Define the accuracy of the calculations */

		void		SetStopping(StoppingPower s)				{ fStopping = s;			}	/*!< Define the stopping power (dE/dX) for the integration */

		void		SetDetectorEff(TGraph* p)				{ fDetectorEff = p;			}	/*!< Define the detector efficiency histogram (efficiency vs theta) */

		void		UseEfficiency(bool b = true)				{ fUseEfficiency = b;			}	/*!< Use the detector efficiencies for particle detection */

		void		SetProjectileExcitation(bool b = true)			{ fProjectileExcitation = b;		} 	/*!< Determines whether the beam or target is being excited */

	private:

		double		IntegrateThetaEnergy(int s);	

		double		fAccuracy;

		Nucleus		*fNucleus;
		Reaction	*fReaction;

		StoppingPower	fStopping;
	
		double		thetamin;
		double		thetamax;
		
		double		energymin;
		double		energymax;

		int 		nTheta;
		int		nEnergy;

		double		meanThetaLab;
		double		meanThetaCM;
		double		meanEnergy;

		double		thetaminCM;
		double		thetamaxCM;
	
		Integral*	fIntegral;

		bool		targetDetection;
		bool		verbose;

		TVectorD		NewIntegratedCrossSection;
		TVectorD		IntegratedCrossSection_TVec;

		int		nThreads;

		bool		fUseEfficiency;
		bool		fProjectileExcitation;
	
		TGraph		*fDetectorEff;

};
#endif
