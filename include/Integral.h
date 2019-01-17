#ifndef Integral_h
#define Integral_h

#include "StatisticalTensor.h"
#include "Nucleus.h"
#include "Reaction.h"
#include "PointCoulEx.h"

#include <thread>
#include <iomanip>
#include <chrono>
#include "TVectorD.h"

///
///	\class Integral
///
///	\brief Perform the point calculations at the user defined energy
///	and theta meshpoints
///
///	PointCoulEx can only determine the point Coulomb excitation 
///	probabilities. In order to calculate a yield for a real experiment 
///	with a range of detection angles and allowing for energy loss within 
///	the target, the Integral class performs PointCoulEx calculations 
///	at a number of theta and energy meshpoints, as defined by the user.
///
///	The class supports some degree of multithreading, albeit relatively 
///	crude.	
///
///	Note that the class name is something of a misnomer. Integral only
///	calculates the meshpoint values. The integration itself is performed
///	in the ExperimentRange class.
///

class Nucleus;
class Reaction;
class PointCoulEx;

class Integral{

	public:
		Integral();
		Integral(Nucleus*,Reaction*);			/*!< Construct the Integral object with a Nucleus and Reaction */
		Integral(const Integral&);			/*!< Copy constructor */
		Integral& operator = (const Integral&);		/*!< Assigment operator */
		~Integral() { ; }

		void	SetThetaMeshpoint(int i, double t);		/*!< Add a new theta meshpoint. If i is smaller than the current size of the meshpoints vector, t is added at position i. Otherwise, the meshpoint t is appended to the vector.  */
		void	SetEnergyMeshpoint(int i, double e);		/*!< Add a new energy meshpoint. If i is smaller than the current size of the meshpoints vector, e is added at position i. Otherwise, the meshpoint t is appended to the vector. */

		void	AddThetaMeshpoint(double t)	{ theta_meshpoints.push_back(t);	}	/*!< Append meshpoint t to the theta meshpoints vector */
		void	AddEnergyMeshpoint(double e)	{ energy_meshpoints.push_back(e);	}	/*!< Append meshpoint e to the energy meshpoints vector */

		void	ClearThetaMeshpoints()		{ theta_meshpoints.clear();		}	/*!< Clear the theta meshpoint vector */
		void	ClearEnergyMeshpoints()		{ energy_meshpoints.clear();		}	/*!< Clear the energy meshpoint vector */

		void	CalculateIntegral();								/*!< Perform the integration over the theta and energy meshpoints */

		std::vector< TVectorD >				GetProbabilities();			/*!< Returns the vector of TVectorD values containing the probabilities at each meshpoint */
		std::vector< std::vector < TVectorD > >		GetMeshPointCrossSections() 	const	{ return meshpointCrossSections;	}	/*!< Return a vector (energies) of vectors (theta) of TVectorD objects containing the calculated cross sections */
		std::vector< std::vector < TVectorD > >		GetMeshPointProbabilities() 	const	{ return meshpointProbabilities;	}	/*!< Return a vector (energies) of vectors (theta) of TVectorD objects containing the calculated probabilities */
		std::vector< std::vector < double   > >		GetCMThetaPoints()		const	{ return cmTheta;			}	/*!< Return a vector (energies) of vectors (theta) containing the center of mass angles of each meshpoint */

		int	GetNthetaMeshpoints()			const	{ return (int)theta_meshpoints.size();	}	/*!< Return the number of theta meshpoints */
		int	GetNenergyMeshpoints()			const	{ return (int)energy_meshpoints.size();	}	/*!< Return the number of energy meshpoints */

		double	GetTheta(int i)				const	{ return theta_meshpoints.at(i);	}	/*!< Return the theta meshpoint at index i */
		double	GetEnergy(int i)			const	{ return energy_meshpoints.at(i);	}	/*!< Return the energy meshpoint at index i */

		void	SetTargetDetection(bool b)			{ targetDetection = b;			}	/*!< Define whether beam (default) or target detection is being performed */

		void	SetVerbose(bool b = true)			{ verbose = b;				}	/*!< Define verbocity of calculations */
	
		void	SetNthreads(int n)				{ nThreads = n;				}	/*!< Define the number of cores the class can use in calculating the meshpoints */
		int	GetNthreads()				const	{ return nThreads;			}	/*!< Return the number of cores the class can use in calculating the meshpoints */

		void	SetAccuracy(double acc)				{ fAccuracy = acc;			}	/*!< Define the accuracy of the calculation */

		void	SetProjectileExcitation(bool b = true)		{ fProjectileExcitation = b;		}	/*!< Defines whether beam or target excitation is being calculated */

		bool	IntegralComplete()				{ return fComplete;			}	/*!< Returns flag indicating whether the integration has been performed */

	private:

		Nucleus*			fNucleus;
		Reaction*			fReaction;

		std::vector<double>		theta_meshpoints;
		std::vector<double>		energy_meshpoints;

		std::vector<PointCoulEx>	point_calculations;
		std::vector<Reaction>		energymeshpoint_reaction;
		std::vector<double>		track_theta;

		std::vector< std::vector < TVectorD > >		meshpointCrossSections;
		std::vector< std::vector < TVectorD > >		meshpointProbabilities;
		std::vector< std::vector < double > >		cmTheta;

		bool				fProjectileExcitation;
		bool				targetDetection;
		bool				verbose;

		int				nThreads;
		double				fAccuracy;

		bool				fComplete;

		// For every PointCoulEx calculation we will need to determine statistical tensors
		std::vector<StatisticalTensor>		fTensors;

};
#endif
