#ifndef StoppingPower_h
#define StoppingPower_h

#include "TGraph.h"
#include "TF1.h"

///
///	\class StoppingPower
///
///	\brief Class for defining stopping powers for use in the integration process
///
///	A simple class which takes user provided stopping powers and fits them to
///	give a stopping power distribution. This is combined with energy integrated 
///	meshpoint calculations to determine accurate experimental yields.
///

class StoppingPower {

	public:
		StoppingPower(); 
		~StoppingPower() {;}
		StoppingPower(const StoppingPower& s);			/*!< Copy constructor */
		StoppingPower& operator = (const StoppingPower& s);	/*!< Assignment operator */

		int	StoppingSize()				{ 	return stoppingE.size();	}	/*!< Return the number of stopping powers defined */

		void	ClearStoppingPowers()			{ 
									stoppingE.clear();		
									stoppingP.clear();  
									fitted = false;	
								}						/*!< Delete the stored stopping power values */
		void	AddStoppingPower(double E, double P)	{ 	
									stoppingE.push_back(E); 	
									stoppingP.push_back(P); 
									fitted = false;	
								}						/*!< Add a new stopping power value */

		void	FitStoppingPowers();	/*!< Perform a fit to the provided stopping power values for easy and accurate interpolation */

		TF1	GetStoppingFit()			{ 
									if(fitted)
										return fit;				
									else{
										FitStoppingPowers();
										return fit;
									}
								}						/*!< Return a fit to the user-provided stopping power values */

	private:
		
		bool			fitted;
		std::vector<double>	stoppingE;
		std::vector<double>	stoppingP;
		TF1			fit;


};

#endif
