#ifndef ParticleDetector_h
#define ParticleDetector_h

#include "TFile.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include <iostream>

///
///	\class ParticleDetector
///
///	\brief Generic class to defined particle detection efficiencies 
///	in theta and phi for corrections to the integrated yields
///

class ParticleDetector {

	public:
		ParticleDetector() 	{ setFlag = false;}
		virtual ~ParticleDetector() 	{;}
		ParticleDetector(const ParticleDetector &d);		/*!< Copy constructor */
		ParticleDetector& operator = (const ParticleDetector &d);	/*!< Assignment operator */

		virtual	double		GetEfficiencyTheta(double);		/*!< Get the particle detection efficiency for a given thete (lab frame)  */

		virtual	double		GetThetaMin();				/*!< Find the minimum theta covered by the detector. For use to define the integration range. */
		virtual	double		GetThetaMax();				/*!< Find the maximum theta covered by the detector. For use to define the integration range. */

		virtual	TGraph*		GetThetaEfficiencyGraph()		const	{ return gThetaEff;	}	/*!< Return the 1D histogram, defining detection efficiency in theta */
		virtual	TGraph2D*	GetThetaPhiMap()			const	{ return gThetaPhi;	}	/*!< Return the 2D histogram, defining detection efficiency in theta and phi */

		virtual void		WriteParticleDetector(const char* filename, const char* opt = "UPDATE");	/*!< Write the particle detection efficiencies to file, filename, with options, opt */


	protected:

		TGraph2D	*gThetaPhi;					/*!< Histogram defining particle detection efficiency in theta and phi */
		TGraph		*gThetaEff;					/*!< Histogram defining particle detection efficiency in theta */
		bool		setFlag;					/*!< Indicates whether the efficiency maps have been defined */

};
#endif
