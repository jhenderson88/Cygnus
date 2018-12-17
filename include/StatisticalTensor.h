#ifndef StatisticalTensor_h
#define StatisticalTensor_h

#include <cstddef>
#include <vector>
#include <iostream>

///
///	\class StateTensor
///
///	\brief Statistical tensor information for a single state
///

class StateTensor {

	public:	
		StateTensor() 	{;}
		~StateTensor()	{;}

		void	SetState(int i)			{ fStateIndex = i;	}			/*!< Define the state */
		void	AddElement(double rho, double k, double kappa);					/*!< Define the density matrix component, rho for given k and kappa */

		int	GetNelements()		const	{ return (int)fTensor.size();		}	/*!< Return the number of components to the tensor */

		int	GetState()		const	{ return fStateIndex;			}	/*!< Return the state index */
		double	GetTensor(int i)	const	{ return fTensor.at((size_t)i);		}	/*!< Return the tensor for element i */
		double	GetK(int i)		const	{ return fK.at((size_t)i);		}	/*!< Return the k value for element i */
		double 	GetKappa(int i)		const 	{ return fKappa.at((size_t)i);		}	/*!< Return the kappa value for element i */

		size_t	IndexFromKkappa(double,double);		/*!< For a given k and kappa combination, find the element index */
	
	private:
		int					fStateIndex;
		std::vector<double>		 	fTensor;
		std::vector<double>			fK;
		std::vector<double>			fKappa;

};

///
///	\class StatisticalTensor
///
///	\brief Holder for StateTensor classes containing statistical
///	tensor information for all populated states
///
class StatisticalTensor{

	public:
		StatisticalTensor()	{;}
		~StatisticalTensor()	{;}

		void SetNstates(int n)					{ fStateTensors.clear(); fStateTensors.resize(n);	}	/*!< Define the number of states */

		int GetNstates()					{ return (int)fStateTensors.size();			}	/*!< Return the number of states */

		StateTensor GetStateTensor(int i)		const	{ return fStateTensors.at(i);				}	/*!< Return StateTensor by state index */
		std::vector<StateTensor> GetStateTensors()	const	{ return fStateTensors;					}	/*!< Return vector of StateTensor objects */

		void AddStateTensor(StateTensor t)			{ fStateTensors.push_back(t);				}	/*!< Add a new state tensor */
		void SetStateTensor(int i, StateTensor t)		{ fStateTensors.at(i) = t;				}	/*!< Set a tensors for state i */

	private:
		std::vector<StateTensor>		fStateTensors;


};

#endif
