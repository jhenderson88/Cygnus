#ifndef CoulExFitter_h
#define CoulExFitter_h

#include "CoulExMinFCN.h"
#include "ScalingParameter.h"
#include <chrono>
#include <iterator>

#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "Math/Minimizer.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnMigrad.h" 
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnContours.h"
#include "Minuit2/MnPlot.h"

#include <iostream>

class CoulExMinFCN;

///
///	\class CoulExFitter
///
///	\brief Calls the CoulExMinFCN class along with the ROOT::Minimizer
///	to perform a chi-squared minimization
///
///	A number of methods and algorithms can be used (as defined in the
///	ROOT::Minimizer documentation in more detail):
///
///	Method		| Algorithm
///	-------------------------------------
///	Minuit2		| Migrad
///			| Simplex
///			| Combined (default)
///			| Scan
///			| Fumili2
///	GSLMultiMin	| ConjugateFR
///			| ConjugatePR
///			| BFGS
///			| BFGS2
///			| SteepestDescent
///	GSLMultiFit	|
///	GSLSimAn	|
///	Genetic		|
///

class CoulExFitter {

	public:
		CoulExFitter();
		virtual	~CoulExFitter()	{;}
		
		void	DoFit(const char* method = "Minuit2", const char* algorithm = "Combined" );	/*!< Perform fitting routine with a user defined method and algorithm (default: Minuit2, Combined) */
		
		virtual	void	ClearAll();			/*!< Completely clear all previous input */
		
		void	DefineExperiment(double t);		/*!< Defines a new experiment with t = theta_cm */
		void	AddData(int i, int init, int fin, double c, double e);	/*!< Add experimental data to experiment defined by index i between initial (init) and final (fin) states with yield (c) and uncertainty (e)  */
		
		void	AddLifetime(int,double,double);		/*!< Add literature lifetime data */
		void	AddBranchingRatio(int,int,int,double,double);	/*!< Add literature branching ratio data */
		void	AddMixingRatio(int,int,double,double);		/*!< Add literature mixing ratio data */

		void	AddFittingMatrixElement(int,int,int,double,double,double);	/*!< Add a fitting matrix element */
		void	CreateScalingParameter(std::vector<int>,double,double,double);	/*!< Add a scaling parameter, with common scaling experiments defined by their indices in a vector of int */


		void	SetMatrixElements(std::vector<MatrixElement> m)		   	{ matrixElements = m;			}	/*!< Define vector of fitting MatrixElements */
		std::vector<MatrixElement>	GetMatrixElements() 			{ return matrixElements;		}	/*!< Return vector of fitting MatrixElements */
		void	ClearMatrixElements()						{ matrixElements.clear();		}	/*!< Clear vector of fitting MatrixElements */

		std::vector<ScalingParameter>	GetScalingParameters()			{ return scalingParameters;		}	/*!< Return vector of ScalingParameter objects */
		void	SetScalingParameters(std::vector<ScalingParameter> s)		{ scalingParameters = s;		}	/*!< Define vector of ScalingParameter objects */
		void	AddScalingParameter(ScalingParameter s)				{ scalingParameters.push_back(s);	}	/*!< Append ScalingParameter object to vector */
		void	ClearScalingParameters()					{ scalingParameters.clear();		}	/*!< Clear scaling parameters */

		void	AddPointCalc(PointCoulEx p)					{ pointCalcs.push_back(p);		}	/*!< Append PointCoulEx calculation object */	
		void	SetPointCalcs(std::vector<PointCoulEx> p)		   	{ pointCalcs = p;			}	/*!< Define vector of PointCoulEx calculation objects */
		std::vector<PointCoulEx>	GetPointCalcs() 			{ return pointCalcs;			}	/*!< Return vector of PointCoulEx calculation objects */

		void	SetData(std::vector<ExperimentData> d)			   	{ exptData = d;				}	/*!< Define vector of ExperimentData */
		std::vector<ExperimentData>	GetData() 				{ return exptData;			}	/*!< Return vector of ExperimentData */

		void	SetLitLifetimes(std::vector<LitLifetime> l)		   	{ litLifetimes = l;			}	/*!< Define vectoe of LitLifetime objects */
		std::vector<LitLifetime>	GetLitLifetimes() 			{ return litLifetimes;			}	/*!< Return vector of LitLifetime objects */

		void	SetLitBranching(std::vector<LitBranchingRatio> b) 	   	{ litBranchingRatios = b;		}	/*!< Define vector of LitBranchingRatio objects */
		std::vector<LitBranchingRatio>	GetLitBranching() 			{ return litBranchingRatios;		}	/*!< Return vector of LitBranchingRatio objects */
	
		void	SetLitMixing(std::vector<LitMixingRatio> m)		   	{ litMixingRatios = m;			}	/*!< Define vector of LitMixingRatio objects */
		std::vector<LitMixingRatio>	GetLitMixing() 				{ return litMixingRatios;		}	/*!< Return vector of LitMixingRatio objects */

		std::vector<TMatrixD>		GetEffectiveCrossSection() 		{ return EffectiveCrossSection;		}	/*!< Return vector of effective cross sections (direct population + feeding) */

		void	SetBaseNucleus(Nucleus* nucl)				   	{ fNucleus_Base = *nucl;		}	/*!< Define base (unmodified) Nucleus object */
		Nucleus				GetBaseNucleus() 			{ return fNucleus_Base;			}	/*!< Return base (unmodified) Nucleus object */
	
		void	SetNucleus(Nucleus *nucl)				   	{ fNucleus = *nucl;			}	/*!< Define Nucleus object */
		Nucleus				GetNucleus() 				{ return fNucleus;			}	/*!< Return Nucleus object */	

		void	AddCorrectionFactor(TVectorD);	/*!< Add point calculation correction factors (append) */
		void	SetCorrectionFactor(int i, TVectorD);	/*!< Define point calculation correction factors for experiment i */
		std::vector<TVectorD> GetCorrectionFactors()				{ return correctionFactors;		}	/*!< Return point calculation correction factors */

		void	Print() const;	/*!< Print fitting details (formatted) */

		void	SetMaxIterations(int i)						{ maxIter = i;				}	/*!< Define number of iterations (Minuit) */
		void	SetMaxFunctionCalls(int i)					{ maxCalls = i;				}	/*!< Define number of calls (GSL) */
		void	SetTolerance(double d)						{ fitTolerance = d;			}	/*!< Define required tolerance */

		int	GetMaxIterations()					const	{ return maxIter;			}	/*!< Return number of iterations (Minuit) */
		int	GetMaxFunctionCalls()					const	{ return maxCalls;			}	/*!< Return number of calls (GSL) */
		double	GetTolerance()						const	{ return fitTolerance;			}	/*!< Return required tolerance */

		void	SetNthreads(int n)						{ nThreads = n;				}	/*!< Define number of allowed cores */
		int	GetNthreads()						const 	{ return nThreads;			}	/*!< Return number of allowed cores */
	
		void	SetVerbose(bool b = true)					{ verbose = b;				}	/*!< Define verbocity */
		bool	GetVerbose()						const	{ return verbose;			}	/*!< Return verbocity */

		void	SetPoissonUncertainties(bool b = true)				{ fUsePoisson = b;			}	/*!< 	*/
		bool	UsePoissonUncertainties()				const	{ return fUsePoisson;			}	/*!<	*/

	private:

		//CoulExMinFCN		theFCN;

		std::vector<int>		index;

		std::vector<double>		parameters;			// Matrix elements + scaling factors
		std::vector<double>		par_LL;				// Matrix elements + scaling factors - LOWER LIMIT
		std::vector<double>		par_UL;				// Matrix elements + scaling factors - UPPER LIMIT
		std::vector<MatrixElement>	matrixElements;			// Preset to related parameters to matrix elements
		std::vector<ScalingParameter>	scalingParameters;

		std::vector<TVectorD>		correctionFactors;

		std::vector<PointCoulEx>	pointCalcs;			// Point calculations
		std::vector<ExperimentData>	exptData;			// Experimental data (one vector entry for each data subset)
		std::vector<LitLifetime>	litLifetimes;			// Literature data, lifetimes
		std::vector<LitBranchingRatio>	litBranchingRatios;		// Literature data, branching ratios
		std::vector<LitMixingRatio>	litMixingRatios;		// Literature data, mixing ratios
		double				theErrorDef;

		std::vector<TMatrixD>		EffectiveCrossSection;

		Nucleus				fNucleus;
		Nucleus				fNucleus_Base;

		int				maxIter;
		int				maxCalls;
		double				fitTolerance;

		int				nThreads;

		bool				first;
		bool				verbose;

		bool				fUsePoisson;

};
#endif
