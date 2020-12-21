#include "CoulExSimFitter.h"
#include "TRandom3.h"
#include <ctime>
#include <chrono>
#include <thread>

///
///	\Class MCMCSimWalker
///
///	\brief Markov-Chain Monte Carlo iterator
///
///	Using a Metropolis-Hastings algorithm the MCMCSimWalker class iterates
///	through the parameter space, with a weighting towards the most likely
///	solutions, causing eventual convergence while sampling the probability
///	distribution of the parameters.
///
///	Inherits from CoulExFitter
///	


class MCMCSimWalker : public CoulExSimFitter {

	public:
		MCMCSimWalker();									/*!< Constructor */
		virtual ~MCMCSimWalker()	{;}							/*!< Destructor */

		void	AddBeamFittingMatrixElement(int,int,int,double,double,double,double);	/*!< Add a fitting matrix element */
		void	AddTargetFittingMatrixElement(int,int,int,double,double,double,double);	/*!< Add a fitting matrix element */

		void	DoMCMCFit();								/*!< Perform the MCMC walk through parameter space */

		double	SteppingKernel(double,double,double,double);				/*!< Defines the transition kernel between points in the Markov Chain - by default a normal distribution about the current point */

		void	WriteTrackedParameters(const char*);					/*!< Write the chain to a text file */

		void 	RandomizeStarting();							/*!< Randomize the starting points */
		
	
	private:
	
		std::vector<std::vector<double>>	parTracker;
		std::vector<double>	probTracker;	
		std::vector<double>	ScalingSigma;
		std::vector<double>	BeamMatrixElementSigma;
		std::vector<double>	TargetMatrixElementSigma;
		TRandom3		sampler;
		TRandom3		rand;	

};
