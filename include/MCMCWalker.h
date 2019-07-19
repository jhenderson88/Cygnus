#include "CoulExFitter.h"
#include "TRandom3.h"
#include <ctime>
#include <chrono>
#include <thread>

///
///	\Class MCMCWalker
///
///	\brief Markov-Chain Monte Carlo iterator
///
///	Using a Metropolis-Hastings algorithm the MCMCWalker class iterates
///	through the parameter space, with a weighting towards the most likely
///	solutions, causing eventual convergence while sampling the probability
///	distribution of the parameters.
///
///	Inherits from CoulExFitter
///	


class MCMCWalker : public CoulExFitter {

	public:
		MCMCWalker();									/*!< Constructor */
		virtual ~MCMCWalker()	{;}							/*!< Destructor */

		void	AddFittingMatrixElement(int,int,int,double,double,double,double);	/*!< Add a fitting matrix element */
		void	CreateScalingParameter(std::vector<int>,double,double,double,double);	/*!< Add a scaling parameter, with common scaling experiments defined by their indices in a vector of int */

		void	DoMCMCFit();								/*!< Perform the MCMC walk through parameter space */

		double	SteppingKernel(double,double,double,double);				/*!< Defines the transition kernel between points in the Markov Chain - by default a normal distribution about the current point */

		void	WriteTrackedParameters(const char*);					/*!< Write the chain to a text file */

		void 	RandomizeStarting();							/*!< Randomize the starting points */
		
	
	private:
	
		std::vector<std::vector<double>>	parTracker;
		std::vector<double>	probTracker;	
		std::vector<double>	ScalingSigma;
		std::vector<double>	MatrixElementSigma;
		TRandom3		sampler;
		TRandom3		rand;	

};
