#ifndef Reaction_h
#define Reaction_h

#include "TGraph.h"
#include "Math/SpecFunc.h"
#include "TMath.h"
#include <iostream>

///
///	\class Reaction
///	
///	\brief Define the reaction kinematics
///
///	This class deals with the reaction kinematics of the Coulomb 
///	excitation process. The kinematics conversion methods borrow 
///	heavily from the GRSISort/GRUTinizer TReaction classes:
///	https://github.com/GRIFFINCollaboration/GRSISort &
///	https://github.com/pcbend/GRUTinizer
///
///	The class is also used to hold a number of useful physical 
///	constants and do define some other common parameters.
///
///	Particle number syntax:
///	1: incoming beam
///	2: target particle
///	3: beam like ejectile
///	4: target like recoil
///
///	Note that center-of-mass angles are always given with respect
///	to the	beam-like ejectile (particle 2). 
///
///

class Reaction
{

	public:
		Reaction();
		
		///
		///	Define the reaction:
		///
		///	aBeam and zBeam are the beam A and Z. aTarget and zTarget are the 
		///	target A and Z. eBeam is the beam energy in MeV (lab frame).	
		///
		Reaction(int aBeam, int zBeam, int aTarget, int zTarget, double eBeam);
		///
		///	Define the reaction:
		///
		///	mBeam and zBeam are the beam mass and Z. mTarget and zTarget are the 
		///	target mass and Z. eBeam is the beam energy in MeV (lab frame).
		///
		///	Masses are defined in atomic mass units	
		///
		Reaction(double mBeam, int zBeam, double mTarget, int zTarget, double eBeam);
		Reaction(const Reaction& r);			/*!< Copy constructor */
		Reaction& operator = (const Reaction& r);	/*!< Assignment operator */
		~Reaction() {;}

		// Basic setters
		void 	SetBeamA(int a)			{ beamA = a;	}	/*!< Set the beam A */
		void 	SetBeamZ(int z)			{ beamZ = z;	}	/*!< Set the beam Z */
		void 	SetTargetA(int a)		{ targetA = a;	}	/*!< Set the target A */
		void 	SetTargetZ(int z)		{ targetZ = z;	}	/*!< Set the target Z */
		void 	SetLabEnergy(double);					/*!< Set the beam energy in the laboratory frame (MeV) */
		void 	SetCMEnergy(double);					/*!< Set the beam energy in the center-of-mass frame (MeV) */

		// Basic getters
		int	GetBeamA()			const		{ return beamA;		}	/*!< Return the beam A */
		int	GetBeamZ()			const		{ return beamZ;		}	/*!< Return the beam Z */
		int	GetTargetA()			const		{ return targetA;	}	/*!< Return the target A */
		int	GetTargetZ()			const		{ return targetZ;	}	/*!< Return the target Z */
		double 	GetLabEnergy()			const		{ return Elab;		}	/*!< Return the beam energy in the laboratory frame (MeV) */
		double 	GetCMEnergy()			const		{ return Ecm;		}	/*!< Return the beam energy in the center-of-mass frame (MeV) */
		double	GetBeta()			const		{ return beta;		}	/*!< Return the beam v/c */

		// Wigner 3j symbol - from GNU
		double 	ThreeJ(int,int,int,int,int,int) const;						/*!< Wigner 3J symbol from GSL */

		// Distance of closest approach
		double 	ClosestApproach()		const		{ return ClosestApproach(beamA,targetA,beamZ,targetZ,Elab)	;}	/*!< Return the distance of closest approach using stored variables */
		double 	ClosestApproach(int aBeam, int aTarget, int zBeam, int zTarget, double eLab) 	const;					/*!< Return the distance of closest approach from aBeam, aTarget, zBeam, zTarget, eLab (MeV) */
	
		// Calculations related to the relative velocity
		///
		///	Return &eta, the wavenumber of a state of energy eState (MeV), define as:
		///	
		///	\f$ \eta=\frac{Z_1Z_2\sqrt{A_1}}{6.34977} * \sqrt(Ebeam - s\cdotEstate) \f$
		///
		///	where:
		///
		///	 \f$ s = (1 + A_1/A2)\f$ 
		///
		double 	EtaCalc(double eState)			const;
		double	RelativeVelocity(double eState) 	const;		/*!< Return the relative velocity for a given state. Presently unused. */

		// Useful constants
		static 	double 	hbarc; 					/*!< MeV fm */
		static 	double 	finestruc; 				/*!< Fine structure constant */
		static 	double 	nuclearmagneton; 			/*!< Nuclear magneton e.fm */
		static 	double	dipole; 				/*!< Normalisation factor for electric dipole polarization - nominally 0.005 */
		static 	double 	electronCharge;				/*!< Electron charge */

		// Rutherford cross section calculations
		double	Rutherford(double theta_lab, int part = 2);	/*!< Return the Rutherford cross section (mb) for theta defined in the lab frame */
		double	RutherfordCM(double theta_cm);			/*!< Return the Rutherford cross section (mb) for theta defined in the CoM frame */

		// Simple printing function - give user information
		void	PrintReaction()			const;		/*!< Print the details of the reaction */

		// Masses (NOT A) for use in reaction calculations. If masses undefined, assumed to be A
		void	SetMass(double b, double t)			{ beamMass = b; targetMass = t;	InitReaction();			}	/*!< Set the beam (b) and target (t) masses in atomic mass units */
		double 	GetBeamMass()			const		{ return beamMass;						}	/*!< Return the beam mass in atomic mass units */
		double 	GetTargetMass()			const		{ return targetMass;						}	/*!< Return the target mass in atomic mass units */
		double 	GetBeamMassMeV()		const		{ return beamMass * 931.494;					}	/*!< Return the beam mass in MeV/c2 */	
		double 	GetTargetMassMeV()		const		{ return targetMass * 931.494;					}	/*!< Return the beam mass in MeV/c2 */	

		// Excitation energy for kinematics
		void	SetExcitationEnergy(double e)			{ exE = e; fQVal = - exE; InitReaction();			}	/*!< Set the excitation energy for the calculation of reaction kinematics */
		double 	GetExcitationEnergy()		const		{ return exE;							}	/*!< Return the excitation energy */

		// Kinematic calculators
		void	InitReaction();			/*!< Set up the reaction kinematics calculations */
		void	SetCmFrame(double);		/*!< Define parameters in the center of mass frame */

		// GOSIA kinematic flags
		/// Sets whether the GOSIA kinematics will be used. GOSIA kinematics differ from (for example) 
		/// Catkin and the normal kinematics models used here. Primarily intended for debugging and
		/// providing a one-to-one comparison with GOSIA.
		void	SetGOSIAKinematics(bool b = true)	{ fGOSIAKin = b;	}	
		bool	GOSIAKinematics()		const	{ return fGOSIAKin;	}	/*!< Return bool indicating whether GOSIA kinematics are being used */


		double	ConvertThetaCmToLab(double, int);		/*!< Convert theta in the center of mass frame to the lab frame (radians) */
		double	ConvertThetaLabToCm(double, int);		/*!< Convert theta in the lab frame to the center of mass frame (radians) */
		TGraph*	ThetaVsTheta(double,double,int);		/*!< Plot theta lab vs theta center of mass */
		TGraph* RutherfordGraph(double,double,int);		/*!< Plot the Rutherford cross section */

	private:
		int beamA;
		int beamZ;
		int targetA;
		int targetZ;

		double beamMass;
		double targetMass;

		double Elab;
		double Ecm;
		double beta;

		double exE;
		double fQVal;

		// Taken from GRUTinizer/GRSISort TReaction code:
		double	fM[4];

		// 	CM frame motion
		double 	fS;
		double 	fInvariantMass;
		double 	fCmTi;
		double 	fCmTf;
		double 	fCmE;
		double	fCmV;
		double	fCmP;
		double	fCmG;

		//	Particles in CM frame
		double	fTCm[4];
		double	fECm[4];
		double	fPCm[4];
		double	fVCm[4];
		double	fGCm[4];

		//	Particles in Lab frame
		//	Note that in the lab frame only the initial beam/target state is fixed
		double	fTLab[2];
		double	fELab[2];
		double	fPLab[2];
		double	fVLab[2];
		double	fGLab[2];
		double	fThetaMax[4];

		bool	reactionSet;

		bool	fGOSIAKin;
		// 	Variables for GOSIA kinematics

		double	fVinf;	// Velocity at infinity
		double	fAred;	// Reduced mass
		double	fEmax;	// Maximum excitation energy
		double	fEPmin;
		double	fTauP;
		double	fTau;

};
#endif
