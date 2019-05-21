#ifndef CoulExSimMinFCN_h
#define CoulExSimMinFCN_h

#include "ScalingParameter.h"
#include "Nucleus.h"
#include "Minuit2/FCNBase.h"
#include "Literature.h"
#include "ExperimentalInput.h"
#include "MatrixElement.h"
#include "TransitionRates.h"
#include "PointCoulEx.h"
#include "Reaction.h"
#include "GammaYield.h"

#include <ctime>
#include <thread>
#include <iomanip>
#include <vector>
#include <chrono>

class MatrixElements;
class ExperimentData;
class LitLifetime;
class LitBranchingRatio;
class LitMixingRatio;

///
///	\class CoulExSimMinFCN
///
///	\brief Contains the definition of the chi-squared function used in simultaneous minimization
///	of beam and target nuclei
///
///	Normalization to target excitation is a powerful method for the study of radioactive ion beams
///	with Coulomb excitation. The normalization for the target and beam are the same (same target
///	thickness, beam intensity, etc.), so they can be fitted simultaneously.
///
///	This class is very similar to the "normal" CoulExMinFCN, except that the inputs and parameters
///	are provided for both the beam and target, with the exception of common scaling parameters.
///
class CoulExSimMinFCN { // : public ROOT::Minuit2::FCNBase{

	public : 

		CoulExSimMinFCN(std::vector<ExperimentData> d_beam, std::vector<ExperimentData> d_target) 
											{ 
												exptData_Beam	= d_beam;
												exptData_Target	= d_target;													
												verbose = false;  
												iter = 0;	
												nThreads = 1;
												fLikelihood 	= false;
											}	/*!< Construct object with vector of experimental data to be fit */
		virtual ~CoulExSimMinFCN()						{;					}

		void	SetupCalculation();	/*!< Prepare the calculation */

		virtual void ClearAll();	/*!< Clear all vectors */

		std::vector<ScalingParameter>	GetScalingParameters()			{ return scalingParameters;		}	/*!< Return the vector of ScalingParameter objects for fitting */
		void	SetScalingParameters(std::vector<ScalingParameter> s)		{ scalingParameters = s;		}	/*!< Define the vector of ScalingParameter objects for fitting */
		void	AddScalingParameter(ScalingParameter s)				{ scalingParameters.push_back(s);	}	/*!< Append a new ScalingParameter to the vector */
		void	ClearScalingParameters()					{ scalingParameters.clear();		}	/*!< Clear the vector of ScalingParameter objects */

		double up() const	 						{ return theErrorDef;			}	/*!< Required by ROOT::Minimizer */
		double operator()(const double*);											/*!< Required by ROOT::Minimizer */
		void	setErrorDef(double def)						{ theErrorDef = def;			}	/*!< Required by ROOT::Minimizer */

		// The below are duplicated for beam and target excitation:
		void	SetBeamMatrixElements(std::vector<MatrixElement> m)		{ ME_Beam = m;				}	/*!< Define the vector of beam MatrixElement objects to be fitted */
		std::vector<MatrixElement>	GetBeamMatrixElements() 	const	{ return ME_Beam;			}	/*!< Return the vector of beam MatrixElement objects to be fitted */
		void	SetTargetMatrixElements(std::vector<MatrixElement> m)		{ ME_Target = m;			}	/*!< Define the vector of target MatrixElement objects to be fitted */
		std::vector<MatrixElement>	GetTargetMatrixElements() 	const	{ return ME_Target;			}	/*!< Return the vector of target MatrixElement objects to be fitted */

		void	SetBeamPointCalcs(std::vector<PointCoulEx> p)			{ pointCalcs_Beam = p;			}	/*!< Define the beam vector of PointCoulEx objects (one for each experiment) */
		std::vector<PointCoulEx>	GetBeamPointCalcs() 		const	{ return pointCalcs_Beam;		}	/*!< Return the beam vector of PointCoulEx objects (one for each experiment) */
		void	SetTargetPointCalcs(std::vector<PointCoulEx> p)			{ pointCalcs_Target = p;		}	/*!< Define the target vector of PointCoulEx objects (one for each experiment) */
		std::vector<PointCoulEx>	GetTargetPointCalcs() 		const	{ return pointCalcs_Target;		}	/*!< Return the target vector of PointCoulEx objects (one for each experiment) */

		void	SetBeamData(std::vector<ExperimentData> d)			{ exptData_Beam = d;			}	/*!< Define the vector of beam ExperimentData objects (one for each experiment) */
		std::vector<ExperimentData>	GetBeamData() 			const	{ return exptData_Beam;			}	/*!< Return the vector of beam ExperimentData objects (one for each experiment) */
		void	SetTargetData(std::vector<ExperimentData> d)			{ exptData_Target = d;			}	/*!< Define the vector of target ExperimentData objects (one for each experiment) */
		std::vector<ExperimentData>	GetTargetData()			const	{ return exptData_Target;		}	/*!< Return the vector of target ExperimentData objects (one for each experiment) */

		void	SetBeamLitLifetimes(std::vector<LitLifetime> l)			{ litLifetimes_Beam = l;		}	/*!< Define the vector of beam LitLifetime objects defining the literature lifetime data for fitting */
		std::vector<LitLifetime>	GetBeamLitLifetimes() 		const	{ return litLifetimes_Beam;		}	/*!< Return the vector of beam LitLifetime objects defining the literature lifetime data for fitting */
		void	SetTargetLitLifetimes(std::vector<LitLifetime> l)		{ litLifetimes_Target = l;		}	/*!< Define the vector of target LitLifetime objects defining the literature lifetime data for fitting */
		std::vector<LitLifetime>	GetTargetLitLifetimes() 	const	{ return litLifetimes_Target;		}	/*!< Return the vector of target LitLifetime objects defining the literature lifetime data for fitting */

		void	SetBeamLitBranching(std::vector<LitBranchingRatio> b) 		{ litBranchingRatios_Beam = b;		}	/*!< Define the vector of beam LitBranchingRatio objects defining the literature branching ratio data for fitting */
		std::vector<LitBranchingRatio>	GetBeamLitBranching() 		const	{ return litBranchingRatios_Beam;	}	/*!< Return the vector of beam LitBranchingRatio objects defining the literature branching ratio data for fitting */
		void	SetTargetLitBranching(std::vector<LitBranchingRatio> b) 	{ litBranchingRatios_Target = b;	}	/*!< Define the vector of target LitBranchingRatio objects defining the literature branching ratio data for fitting */
		std::vector<LitBranchingRatio>	GetTargetLitBranching() 	const	{ return litBranchingRatios_Target;	}	/*!< Return the vector of target LitBranchingRatio objects defining the literature branching ratio data for fitting */
	
		void	SetBeamLitMixing(std::vector<LitMixingRatio> m)			{ litMixingRatios_Beam = m;		}	/*!< Define the vector of beam LitMixingRatio objects defining the literature mixing ratio data for fitting */
		std::vector<LitMixingRatio>	GetBeamLitMixing() 		const	{ return litMixingRatios_Beam;		}	/*!< Return the vector of beam LitMixingRatio objects defining the literature mixing ratio data for fitting */
		void	SetTargetLitMixing(std::vector<LitMixingRatio> m)		{ litMixingRatios_Target = m;		}	/*!< Define the vector of target LitMixingRatio objects defining the literature mixing ratio data for fitting */
		std::vector<LitMixingRatio>	GetTargetLitMixing() 		const	{ return litMixingRatios_Target;	}	/*!< Return the vector of target LitMixingRatio objects defining the literature mixing ratio data for fitting */

		std::vector<TMatrixD>	GetEffectiveCrossSection_Beam()		const	{ return EffectiveCrossSection_Beam;	}	/*!< Return the beam's "effective cross section" = direct population + feeding */
		std::vector<TMatrixD>	GetEffectiveCrossSection_Target() 	const	{ return EffectiveCrossSection_Target;	}	/*!< Return the target's "effective cross section" = direct population + feeding */

		void	SetBaseBeamNucleus(Nucleus* nucl)				{ fNucleus_Target_Base = *nucl;		}	/*!< Define the base beam nucleus (not to be varied in fitting) */
		Nucleus				GetBaseBeamNucleus() 		const	{ return fNucleus_Beam_Base;		}	/*!< Return the base beam nucleus (not to be varied in fitting) */
		void	SetBaseTargetNucleus(Nucleus* nucl)				{ fNucleus_Target_Base = *nucl;		}	/*!< Define the base target nucleus (not to be varied in fitting) */
		Nucleus				GetBaseTargetNucleus() 		const	{ return fNucleus_Target_Base;		}	/*!< Return the base target nucleus (not to be varied in fitting) */
	
		void	SetBeamNucleus(Nucleus *nucl)					{ fNucleus_Beam = *nucl;		}	/*!< Define the fitting beam nucleus (varied in fitting) */
		Nucleus				GetBeamNucleus() 		const	{ return fNucleus_Beam;			}	/*!< Return the fitting beam nucleus (varied in fitting) */
		void	SetTargetNucleus(Nucleus *nucl)					{ fNucleus_Target = *nucl;		}	/*!< Define the fitting target nucleus (varied in fitting) */
		Nucleus				GetTargetNucleus() 		const	{ return fNucleus_Target;		}	/*!< Return the fitting target nucleus (varied in fitting) */

		void	SetBeamCorrectionFactors(std::vector<TVectorD> v)		{ correctionFactors_Beam = v;		}	/*!< Define the vector of correction factors between point and integrated cross sections for the beam */
		std::vector<TVectorD> GetBeamCorrectionFactors()		const	{ return correctionFactors_Beam;	}	/*!< Return the vector of correciton factors between point and integrated cross sections for the beam */
		void	SetTargetCorrectionFactors(std::vector<TVectorD> v)		{ correctionFactors_Target = v;		}	/*!< Define the vector of correction factors between point and integrated cross sections for the target */
		std::vector<TVectorD> GetTargetCorrectionFactors()		const	{ return correctionFactors_Target;	}	/*!< Return the vector of correciton factors between point and integrated cross sections for the target */

		// The following is unchanged from CoulExMinFCN:
		void	SetNpar(int n)							{ nPar = n;				}	/*!< Define the number of fitting parameters (fitting matrix elements + scaling parameters) */
		int	GetNpar()						const	{ return nPar;				}	/*!< Return the number of fitting parameters (fitting matrix elements + scaling parameters) */

		void	SetVerbose(bool b = true)					{ verbose = b;				}	/*!< Define the verbocity of the minimization */	
		bool	GetVerbose() 						const	{ return verbose;			}	/*!< Return the verbocity of the minimization */

		void	SetIter(int i)							{ nIterations = i;			}	/*!< Define the number of iterations (MINUIT2) */
		void 	SetCalls(int i)							{ nCalls = i;				}	/*!< Define the number of function calls (GSL) */

		void	SetNthreads(int n)						{ nThreads = n;				}	/*!< Define the number of cores the function is allowed to use */
		int	GetNthreads()						const	{ return nThreads;			}	/*!< Return the number of cores the function is allowed to use */

		double	GetParameter(int i)					const	{ return parameters.at(i);		}	/*!< Return the fitting parameter indexed by i */
		std::vector<double> GetParameters()				const	{ return parameters;			}	/*!< Return the vector of fitting parameters */

		void	ResetIter()							{ iter = 0;				}	/*!< Reset the iteration number */

		void	SetLikelihoodFit(bool b = true)					{ fLikelihood = b;			}	/*!< Define whether we do a log-likelihood based fit (default: chi-squared) */
		bool	LikelihoodFit()						const	{ return fLikelihood;			}	/*!< Return whether we do a log-likelihood based fit (default: chi-squared) */

	private :

		std::vector<double>		parameters;			/*!< Matrix elements for both beam and target, and common scaling factors */
		std::vector<MatrixElement>	ME_Beam;			/*!< Beam matrix elements - Preset to relate parameters to beam matrix elements */
		std::vector<MatrixElement>	ME_Target;			/*!< Target matrix elements - Preset to relate parameters to target matrix elements */
		std::vector<ScalingParameter>	scalingParameters;		/*!< Scaling parameters - common to both target and beam excitations */

		std::vector<TVectorD>		correctionFactors_Beam;		/*!< Point correction factor for the beam (accounts for the angular distribution of the cross section) */
		std::vector<TVectorD>		correctionFactors_Target;	/*!< Point correction factor for the target (accounts for the angular distribution of the cross section) */

		std::vector<PointCoulEx>	pointCalcs_Beam;		/*!< Point calculations for beam excitation */
		std::vector<PointCoulEx>	pointCalcs_Target;		/*!< Point calculations for target excitation */

		std::vector<ExperimentData>	exptData_Beam;			/*!< Beam excitation experimental data (one vector entry for each data subset) */
		std::vector<ExperimentData>	exptData_Target;		/*!< Target excitation experimental data (one vector entry for each data subset) */

		std::vector<LitLifetime>	litLifetimes_Beam;		/*!< Literature data for the beam, lifetimes */
		std::vector<LitLifetime>	litLifetimes_Target;		/*!< Literature data for the target, lifetimes */

		std::vector<LitBranchingRatio>	litBranchingRatios_Beam;	/*!< Literature data for the beam, branching ratios */
		std::vector<LitBranchingRatio>	litBranchingRatios_Target;	/*!< Literature data for the target, branching ratios */

		std::vector<LitMixingRatio>	litMixingRatios_Beam;		/*!< Literature data for the beam, mixing ratios */
		std::vector<LitMixingRatio>	litMixingRatios_Target;		/*!< Literature data for the target, mixing ratios */

		double				theErrorDef;

		std::vector<TMatrixD>		EffectiveCrossSection_Beam;
		std::vector<TMatrixD>		EffectiveCrossSection_Target;

		Nucleus				fNucleus_Beam;
		Nucleus				fNucleus_Target;
		Nucleus				fNucleus_Beam_Base;
		Nucleus				fNucleus_Target_Base;

		int				nPar;

		TransitionRates			fRates;

		bool				verbose;

		int				nCalls;
		int				nIterations;
		int				iter;

		int				nThreads;

		std::vector<int>		exptIndex;

		bool				fLikelihood;

};

#endif
