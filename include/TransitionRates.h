#ifndef TransitionRates_h
#define TransitionRates_h

#include "MiscFunctions.h"
#include "Nucleus.h"
#include <vector>
#include <iostream>
#include "TMatrixD.h"
#include "TVectorD.h"

class Nucleus;

///
///	\class TransitionRates
///
///	\brief Determine spectroscopic information from matrix elements
///
///	From the matrix elements and energies defined in the Nucleus class
///	this class determines transition rates, lifetimes, mixing ratios
///	and branching ratios. Specifically, this is for comparison with
///	input experimental (literature) values in the minimization process.
///

class TransitionRates
{

	public:
		TransitionRates();
		TransitionRates(Nucleus*);			/*!< Constructor with a Nucleus class */
		~TransitionRates() {;}

		void		SetMatrixElements();		/*!< Set up the experimental data based on the provided Nucleus class */	
		void		Print() const;			/*!< Print the data, as determined from energies and matrix elements */

		TVectorD	GetLifetimes()			const			{ return StateLifetimes;	}	/*!< Return a TVectorD of state lifetimes (ps) */	
		TMatrixD	GetBranchingRatios()		const			{ return BranchingRatios;	}	/*!< Return a TVectorD of branching ratios */
		TMatrixD	GetMixingRatios()		const			{ return MixingRatios;		}	/*!< Return a TMatrixD of mixing ratios (E/M) */

		void		SetNucleus(Nucleus* nucl)	{ fNucleus = nucl; }	/*!< Define the nucleus (energies and matrix elements) */

	private:
		Nucleus*			fNucleus;

		std::vector<TMatrixD>		MatrixElements;
		std::vector<TMatrixD>		TransitionStrengths;
		std::vector<TMatrixD>		TransitionStrengths_Abs;
		TMatrixD			SummedTransitionStrengths;
		TMatrixD			Lifetimes;
		TMatrixD			BranchingRatios;			// Column == initial state, row == final states
		TMatrixD			MixingRatios;
		std::vector<double>		StateJ;
		std::vector<double>		StateE;

		TVectorD			StateLifetimes;
		TVectorD			StateDecayProb;
		
		double				SumColumn(TMatrixD,int) const;
		double				MaxAbsMatrix(TMatrixD)	const;

};

#endif
