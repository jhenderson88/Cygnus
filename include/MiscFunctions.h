#ifndef MiscFunctions_h
#define MiscFunctions_h

#include "PointCoulEx.h"
#include "Reaction.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TGraph.h"
#include "Nucleus.h"
#include <iostream>
#include <fstream>

///
///	\class MiscFunctions
///
///	\brief Useful static functions for use elsewhere in the software package
///

class PointCoulEx;

class MiscFunctions {

	public:

		static	double	GetMaxMatrix(TMatrixD);			/*!< Return the maximum value in a TMatrixD */	
		static	double	GetMaxAbsMatrix(TMatrixD);		/*!< Return the maximum absolute value in a TMatrixD */
		static	void	WriteMatrixNucleus(std::ofstream&, TMatrixD, Nucleus);	/*!< Write a well formatted TMatrixD based on states from a Nucleus object */
		static	void	PrintMatrixNucleus(TMatrixD, Nucleus);			/*!< Print a well formatted TMatrixD based on states from a Nucleus object */
		static	void	PrintVectorNucleus(TVectorD, Nucleus, const char*);	/*!< Print a well formatted TVectorD based on states from a Nucleus object */

		static	double	RotationFunction(double,int,int,int);	/*!< Rotation function to convert Reaction frame tensors to laboratory frame */

		static	double	c;					/*!< Speed of light for use in calculations */	
		static 	double	hbar;					/*!< Planks constant / 2pi for use in calculations */

		static	unsigned int	doublefactorial(unsigned int);	/*!<	*/

		static 	double	SphericalHarmonics(double,int,int,bool b = true);	/*!<	*/

		static	double	SimpsonsRule(TGraph*,int nSteps, double xMin, double xMax);	/*!<	Perform a numerical integration of a TGraph using Simpson's rule using nSteps between xMin and xMax	*/

		static	TGraph*	PlotCrossSection(PointCoulEx*, double, double, int, int, bool cm = true, int nPart = 2);	/*!<	Plot the excitation cross section to a state over a range of angles in the center of mass (cm = true) or lab (cm = false) frame	*/
		static	TGraph*	PlotProbability(PointCoulEx*, double, double, int, int, bool cm = true, int nPart = 2);		/*!<	Plot the excitation probability to a state over a range of angles in the center of mass (cm = true) or lab (cm = false) frame	*/
	
};

#endif
