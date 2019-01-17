#ifndef ParticleDetectorS3_h
#define ParticleDetectorS3_h

#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include <iostream>
#include "ParticleDetector.h"

///
///	\class ParticleDetectorS3
///
///	\brief Example class inheriting from the ParticleDetector class to define 
///	coverage based on an S3 detector
///

class ParticleDetectorS3 : public ParticleDetector{

	public:
		///
		///	Define an S3 detector (or subset thereof) \n
		///	id - define the experiment index, for use in naming conventions of histograms \n
		///	z - the offset in z of the detector from the target position \n
		///	x_off - the offset (from the beam axis) of the S3 detector in the x-direction \n
		///	y_off - the offset (from the beam axis) of the S3 detector in the y-direction \n
		///	ring_in - the innermost ring of the S3 in the experimental range \n
		///	ring_out - the outermost ring of the S3 in this experimental range \n
		///	tmin & tmax - define the maximum and minimum theta values in the efficiency
		///	histograms, this must cover a larger range than the detector 
		///
		ParticleDetectorS3(TGraph2D*, TGraph*);
		ParticleDetectorS3(int id, double z, double x_off, double y_off, int ring_in, int ring_out, double tmin=0, double tmax=180);

};

#endif
