#ifndef GammaYield_h
#define GammaYield_h

#include "ExperimentRange.h"
#include "TransitionRates.h"
#include "Nucleus.h"
#include "TMatrixD.h"
#include "TVectorD.h"

#include <iostream>
#include <fstream>

///
///	\class GammaYield
///
///	\brief Converts cross sections to gamma-ray yields
///

class GammaYield {

	public:
		GammaYield()	{;}
		~GammaYield()	{;}

		static TMatrixD	GammaRayYield(TVectorD, TMatrixD);			/*!< Return a matrix of gamma-ray yields from a TVectorD (cross-section) and TMatrixD (branching ratios) */
		static TMatrixD GammaRayYield(TVectorD v, TransitionRates t)				{ return GammaRayYield(v,t.GetBranchingRatios());					}	/*!< Matrix of gamma-ray yields from TVectorD (cross-sections) and TransitionRates object */
		static TMatrixD GammaRayYield(ExperimentRange r, TransitionRates t)			{ return GammaRayYield(r.GetIntegratedCrossSection_TVec(),t.GetBranchingRatios());	}	/*!< Matrix of gamma-ray yields from an ExperimentRange (cross-sections) and TransitionRates object */

		static void	PrintYields(ExperimentRange r, TransitionRates t, Nucleus n, double scale = 1);		/*!< Print well-formatted yields based on ExperimentRange (cross-sections), TransitionRates (branching ratios) and Nucleus (level labelling) objects */
		static void	WriteYields(ExperimentRange r, TransitionRates t, Nucleus n, std::ofstream& outfile, double scale = 1);		/*!< Write well-formatted yields to file based on ExperimentRange (cross-sections), TransitionRates (branching ratios) and Nucleus (level labelling) objects */

		static double	GetYield(ExperimentRange r, TransitionRates t, Nucleus n, int i, int f);	/*!< Return a yield for a single transition between initial state i and final state f */

};
#endif
