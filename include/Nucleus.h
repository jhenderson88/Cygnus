#ifndef Nucleus_h
#define Nucleus_h

#include <TMath.h>
#include <TMatrixD.h>
#include <iostream>
#include <iomanip>
#include <vector>

///
///	\class Nucleus
///
///	\brief Holds information about the nucleus for which the CoulEx 
///	calculation will be performed
///
///	This class stores the information for the nucleus being investigated. 
///
///	In particular, state energies, spins, parities and matrix elements 
///	are stored here. 
///

class Nucleus
{

	public:

		Nucleus();						
		Nucleus(int Z, int A, int nS, int nL = 7);	/*!< Construct nucleus with Z protons, A nucleons, nS total states and nL multipolarities allowed */				
		Nucleus(const Nucleus& n);			/*!< Copy constructor */
		Nucleus& operator = (const Nucleus& n);		/*!< Assignment operator */
		~Nucleus() {;}
	
		void SetNstates(int);				/*!< Define the number of states in the calculation */
		void SetMaxLambda(int);				/*!< Define the maximum multipolarity to be used: Recommended to use default = 7 */

		void SetState(int s, double E, double J, int P);	/*!< Define state of index s, energy E, spin J and parity P */
		void SetStateE(int s, double E);			/*!< Define energy E of state with index s */
		void SetStateJ(int s, double J);			/*!< Define spin J of state with index s */
		void SetStateP(int s, int P);				/*!< Define parity P of state with index s */
		void SetMatrixElement(int l, int s1, int s2, double me);	/*!< Define matrix element me, of multipolarity l, between states s1 and s2 */

		void SetZ(int z)	{ nucleusZ = z;	}	/*!< Define nucleus proton number, Z */
		void SetA(int a)	{ nucleusA = a;	}	/*!< Define nucleus mass number, A */

		void MultipolarityReminder() const;		/*!< Prints a reminder regarding multipolarity indexing */

		int GetZ()	const		{ return nucleusZ;	}	/*!< Return nucleus proton number */
		int GetA()	const		{ return nucleusA;	}	/*!< Return nucleus mass number	*/

		int			GetNstates()		const	{ return nStates;		}	/*!< Return number of state */
		int			GetMaxLambda()		const	{ return maxLambda;		}	/*!< Return maximum lambda (default = 7) */
	
		std::vector<double>	GetLevelEnergies()	const	{ return LevelEnergies;		}	/*!< Return vector of level energies */
		std::vector<double>	GetLevelJ()		const	{ return LevelJ;		}	/*!< Return vector of level spins */
		std::vector<int>	GetLevelP()		const	{ return LevelP;		}	/*!< Return vector of level parities */
		std::vector<TMatrixD>	GetMatrixElements()	const	{ return MatrixElements;	}	/*!< Return vector of transition matrices (one per multipolarity) */
	
		void			PrintNucleus()	const;		/*!< Print Nucleus information */
		void			PrintState(int) const;		/*!< Print State information */
		
	private:

		int			nStates;			/*!< Number of states in the nucleus */
		int			maxLambda;			/*!< Number of multipolarities to be considered */

		int			nucleusZ;			/*!< Proton number of the nucleus */
		int			nucleusA;			/*!< Nucleon number of the nucleus */

		std::vector<double>	LevelEnergies;			/*!< Vector containing the level energies of the states */
		std::vector<double>	LevelJ;				/*!< Vector containing the spins of the states */
		std::vector<int>	LevelP;				/*!< Vector containing the parities of the states */

		std::vector<TMatrixD>	MatrixElements;			/*!< Vector (length nLambda) containg the matrix of matrix elements, dimensions nStates x nStates */
		std::vector<TMatrixD>	MatrixElementsUL;		/*!< Vector (length nLambda) containg the matrix of matrix element upper limits, dimensions nStates x nStates */
		std::vector<TMatrixD>	MatrixElementsLL;		/*!< Vector (length nLambda) containg the matrix of matrix element lower limits, dimensions nStates x nStates */

};
#endif
