#ifndef MiscFunctions_h
#define MiscFunctions_h

#include "TMatrixD.h"
#include "TVectorD.h"
#include "Nucleus.h"
#include <iostream>
#include <fstream>

///
///	\class MiscFunctions
///
///	\brief Useful static functions for use elsewhere in the software package
///

class MiscFunctions {

	public:

		static	double	GetMaxMatrix(TMatrixD);			/*!< Return the maximum value in a TMatrixD */	
		static	double	GetMaxAbsMatrix(TMatrixD);		/*!< Return the maximum absolute value in a TMatrixD */
		static	void	WriteMatrixNucleus(std::ofstream&, TMatrixD, Nucleus);	/*!< Write a well formatted TMatrixD based on states from a Nucleus object */
		static	void	PrintMatrixNucleus(TMatrixD, Nucleus);			/*!< Print a well formatted TMatrixD based on states from a Nucleus object */
		static	void	PrintVectorNucleus(TVectorD, Nucleus, const char*);	/*!< Print a well formatted TVectorD based on states from a Nucleus object */

		static	double	RotationFunction(double,int,int,int);	/*!< Rotation function to convert Reaction frame tensors to laboratory frame */

		static	double	c;	
		static 	double	hbar;
	
};

#endif
