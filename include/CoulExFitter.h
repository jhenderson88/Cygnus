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
///	By default, the fitter will use the uncertainties provided by the user
///	to perform the minimization. The user can specify, however, to use 
///	a likelihood derived fit based on assumed Poisson statistics by calling
///	the SetPoissonUncertainties() function.
///
///	<b>If using this method, the user should pass non-efficiency corrected yields</b>
///	to the fitter, along with efficiencies (stored in the ExptData class)
///	in order to properly determine the Poisson uncertainties. This option is especially
///	important for low-statistics data, where the assumed symmetric, SQRT(N) behaviour 
///	of the data breaks down.
///
///	By default, the fitter will create a standard covariance matrix using the 
///	minimizer. In reality, most matrix elements will be have asymmetric uncertainties 
///	which must be determined using the MINOS package. To do this, the user should
///	call the DoFullUncertainty() function, prior to performing the minimization.
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
		void	AddMatrixElement(int,int,int,double,double);		/*!< Add literature matrix element data */

		void	AddFittingMatrixElement(int,int,int,double,double,double);	/*!< Add a fitting matrix element */
		void	CreateScalingParameter(std::vector<int>,double,double,double);	/*!< Add a scaling parameter, with common scaling experiments defined by their indices in a vector of int */

		void	SetMatrixElements(std::vector<MatrixElement> m)		   	{ matrixElements = m;			}	/*!< Define vector of fitting MatrixElements */
		std::vector<MatrixElement>	GetMatrixElements() 			{ return matrixElements;		}	/*!< Return vector of fitting MatrixElements */
		void	AddMatrixElement(MatrixElement m)				{ matrixElements.push_back(m);		}
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
	
		void	SetLitMatrixElements(std::vector<LitMatrixElement> m)	   	{ litMatrixElements = m;		}	/*!< Define vector of LitMatrixElement objects */
		std::vector<LitMatrixElement>	GetLitMatrixElements() 			{ return litMatrixElements;		}	/*!< Return vector of LitMatrixElement objects */

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

		void	SetPoissonUncertainties(bool b = true)				{ fUsePoisson = b;			}	/*!< Define whether to use user defined data uncertainties or to use Poisson uncertainties (and likelihood fit) 	*/
		bool	UsePoissonUncertainties()				const	{ return fUsePoisson;			}	/*!< Return whether to use user defined data uncertainties or to use Poisson uncertainties (and likelihood fit) 	*/

		TMatrixD	GetCovarianceMatrix()				const	{ return covMat;			}	/*!< Return covariance matrix from fit 	*/
		TMatrixD	GetCorrelationMatrix()				const	{ return corMat;			}	/*!< Return correlation matrix from fit 	*/

		void	SetDoFullUncertainty(bool b = true)				{ fDoFullUnc = b;			} 	/*!< Define whether to do a complete MINOS uncertainty analysis (slow)	*/
		bool	DoFullUncertainty()					const	{ return fDoFullUnc;			} 	/*!< Return whether to do a complete MINOS uncertainty analysis (slow)	*/

		void	SetLikelihoodFit(bool b = true)					{ fLikelihood = b;			}	/*!< Define whether we do a log-likelihood based fit (default: chi-squared) */
		bool	LikelihoodFit()						const	{ return fLikelihood;			}	/*!< Return whether we do a log-likelihood based fit (default: chi-squared) */

		std::vector<double>	GetFitParameters()				{ return parameters;			}	/*!< Return fit parameters - note that these will be updated with the fit result after the fit - Used in MCMC methods */
		std::vector<double>	GetFitUL()					{ return par_UL;			}	/*!< Return the fit parameter upper limits - Used in MCMC methods */
		std::vector<double>	GetFitLL()					{ return par_LL;			}	/*!< Return the fit parameter lower limits - Used in MCMC methods */

		void	SetFitParameters(std::vector<double> p)				{ parameters = p;			}	/*!< Set the fit parameters - note that these will be updated after the fit has been performed - Used in MCMC methods */
		void	SetFitUL(std::vector<double> p)					{ par_UL = p;				}	/*!< Set the fit parameter upper limits - Used in MCMC methods */
		void	SetFitLL(std::vector<double> p)					{ par_LL = p;				}	/*!< Set the fit parameter lower limits - Used in MCMC methods */

		void	SetFittingParameter(size_t i, double v)				{ parameters.at(i) = v;			}	/*!< Set an individual fitting parameter - Used in MCMC methods */
	
	private:

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
		std::vector<LitMatrixElement>	litMatrixElements;		// Literature data, matrix elements
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

		TMatrixD			covMat;
		TMatrixD			corMat;

		bool				fDoFullUnc;

		bool				fLikelihood;

};
#endif
