#include "CoulExFitter.h"
#include "TRandom3.h"
#include <ctime>
#include <chrono>
#include <thread>

class MCMCWalker : public CoulExFitter {

	public:
		MCMCWalker();
		virtual ~MCMCWalker()	{;}

		void	AddFittingMatrixElement(int,int,int,double,double,double,double);	/*!< Add a fitting matrix element */
		void	CreateScalingParameter(std::vector<int>,double,double,double,double);	/*!< Add a scaling parameter, with common scaling experiments defined by their indices in a vector of int */

		void	DoMCMCFit();

		double	SteppingKernel(double,double,double,double);

		void	WriteTrackedParameters(const char*);

		void 	RandomizeStarting();
		
	
	private:
	
		std::vector<std::vector<double>>	parTracker;
		std::vector<double>	probTracker;	
		std::vector<double>	ScalingSigma;
		std::vector<double>	MatrixElementSigma;
		TRandom3		sampler;
		TRandom3		rand;	

};
