#include "MiscFunctions.h"

double MiscFunctions::c 	= 299792458; 
double MiscFunctions::hbar	= 1.05457e-34;

double MiscFunctions::GetMaxMatrix(TMatrixD mat){

	double max = 0;
	for(int i=0;i<mat.GetNcols();i++){
		for(int j=0;j<mat.GetNrows();j++){
			if(mat[j][i] > max) max = mat[j][i];
		}
	}
	return max;

}

double MiscFunctions::GetMaxAbsMatrix(TMatrixD mat){

	double max = 0;
	for(int i=0;i<mat.GetNcols();i++){
		for(int j=0;j<mat.GetNrows();j++){
			if(mat[j][i] > TMath::Abs(max)) max = TMath::Abs(mat[j][i]);
		}
	}
	return max;

}

void MiscFunctions::WriteMatrixNucleus(std::ofstream& outfile, TMatrixD tmpMat, Nucleus nucl){

	outfile		<< std::setw(12) << std::left << "Index"
			<< std::setw(12) << " ";
	for(int s = 0; s < nucl.GetNstates(); s++)
		outfile 	<< std::setw(12) << std::left << s;
	outfile		<< std::endl;	

	outfile 	<< std::setw(12) << std::left << " "
			<< std::setw(12) << std::left << "State";
	for(int s = 0; s < nucl.GetNstates(); s++)
		outfile		<< std::setw(12) << std::left << nucl.GetLevelJ().at(s);
	outfile		<< std::endl;

	for(int y = 0; y < tmpMat.GetNrows(); y++){
		outfile 	<< std::setw(12) << std::left << y
				<< std::setw(12) << std::left << nucl.GetLevelJ().at(y);
		for(int x = 0; x < tmpMat.GetNcols(); x++){
			outfile		<< std::setw(12) << std::left << tmpMat[y][x];
		}
		outfile		<< std::endl;
	}

}

void MiscFunctions::PrintMatrixNucleus(TMatrixD tmpMat, Nucleus nucl){

	std::cout	<< std::setw(12) << std::left << "Index"
			<< std::setw(12) << " ";
	for(int s = 0; s < nucl.GetNstates(); s++)
		std::cout 	<< std::setw(12) << std::left << s;
	std::cout	<< std::endl;	

	std::cout 	<< std::setw(12) << std::left << " "
			<< std::setw(12) << std::left << "State";
	for(int s = 0; s < nucl.GetNstates(); s++)
		std::cout	<< std::setw(12) << std::left << nucl.GetLevelJ().at(s);
	std::cout	<< std::endl;

	for(int y = 0; y < tmpMat.GetNrows(); y++){
		std::cout 	<< std::setw(12) << std::left << y
				<< std::setw(12) << std::left << nucl.GetLevelJ().at(y);
		for(int x = 0; x < tmpMat.GetNcols(); x++){
			std::cout	<< std::setw(12) << std::left << tmpMat[y][x];
		}
		std::cout	<< std::endl;
	}

}

void MiscFunctions::PrintVectorNucleus(TVectorD tmpVec, Nucleus nucl, const char *var){

	std::cout	<< std::setw(12) << std::left << "Index" 
			<< std::setw(12) << std::left << "State"
			<< std::setw(20) << std::left << var
			<< std::endl;
	for(int y = 0; y < tmpVec.GetNrows(); y++){
		std::cout 	<< std::setw(12) << std::left << y
				<< std::setw(12) << std::left << nucl.GetLevelJ().at(y)
				<< std::setw(20) << std::left << tmpVec[y]
				<< std::endl;
	}

}

double MiscFunctions::RotationFunction(double beta, int k, int kpp, int kp){ // Beta, K, X, X'

	double cosB 	= TMath::Cos(beta/2.);
	double sinB	= TMath::Sin(beta/2.);
	double ctb	= TMath::Power(cosB/sinB,2);

	int ja	= k + kp;
	int jb	= k - kp;
	int jc 	= k + kpp;
	int jd  = k - kpp;
	double b1	= TMath::Factorial(ja) * TMath::Factorial(jb) * TMath::Factorial(jc) * TMath::Factorial(jd);


	ja	= kp + kpp;
	jb	= 2*k - kp - kpp;
	double f = TMath::Power(-1,k-kp) * TMath::Power(cosB,ja) * TMath::Power(sinB,jb) * TMath::Sqrt(b1);
	
	int mis = 0;
	if(ja < 0)
		mis = -ja;
	int mas = k - kpp;
	if(kpp < kp)
		mas = k - kp;
	ja 	= kp + kpp + mis;
	jb 	= k - kpp - mis;
	jc	= k - kp - mis;
	jd	= mis;
	double b2	= TMath::Factorial(ja) * TMath::Factorial(jb) * TMath::Factorial(jc) * TMath::Factorial(jd);

	double g	= TMath::Power(-ctb,mis)/b2;
	double DJMM	= g;
	ja = mis + 1;
	if(mas >= ja){
		for(int j = ja; j<=mas; j++){
			g *= -ctb*(k-kpp-j+1)*(k-kp-j+1)/((kp+kpp+j)*j);
			DJMM += g;
		}
	}
	DJMM *= f;

	double sum = DJMM;	

	return sum;

}

unsigned int MiscFunctions::doublefactorial(unsigned int n){
	
	if( n ==0 || n == 1)
		return 1;
	return	n * doublefactorial(n-2);

}

//
//	Spherical harmonics tabulated at: https://en.wikipedia.org/wiki/Table_of_spherical_harmonics
//
double MiscFunctions::SphericalHarmonics(double theta, int l, int m, bool sym){

	if(sym){
		if(m != 0)
			return 0;
		switch (l){
			case 0:
				return 	0.50 * TMath::Sqrt(1/TMath::Pi());
			case 2:
				return	0.25 * TMath::Sqrt(5 / TMath::Pi()) * (3 * TMath::Power(TMath::Cos(theta),2) - 1);
			case 4:
				return	(3./16.) * TMath::Sqrt(1 / TMath::Pi()) * (35 * TMath::Power(TMath::Cos(theta),4) - 30 * TMath::Power(TMath::Cos(theta),2) + 3);
			case 6:
				return 	(1. / 32.) * TMath::Sqrt(13. / TMath::Pi()) * (231 * TMath::Power(TMath::Cos(theta),6) - 315 * TMath::Power(TMath::Cos(theta),4) + 105 * TMath::Power(TMath::Cos(theta),2) - 5);
			return 1;
		}
	}
	else{
		switch (l){
			case 0:
				return 	0.50 * TMath::Sqrt(1/TMath::Pi());
			case 2:
				switch (m){
					case -2:
						return	0.25 * TMath::Sqrt(15 / (2 * TMath::Pi())) * (TMath::Power(TMath::Sin(theta),2));
					case -1:
						return	0.50 * TMath::Sqrt(15 / (2 * TMath::Pi())) * (TMath::Sin(theta) * TMath::Cos(theta));
					case 0 :
						return	0.25 * TMath::Sqrt(5 / TMath::Pi()) * (3 * TMath::Power(TMath::Cos(theta),2) - 1);
					case 1 :
						return	-0.50 * TMath::Sqrt(15 / (2 * TMath::Pi())) * (TMath::Sin(theta) * TMath::Cos(theta));
					case 2 :
						return	0.25 * TMath::Sqrt(15 / (2 * TMath::Pi())) * (TMath::Power(TMath::Sin(theta),2));
					return 0;
				}
			case 4:
				switch (m){
					case -4:
						return	(3./16.) * TMath::Sqrt(35. / (2 * TMath::Pi())) * TMath::Power(TMath::Sin(theta),4);
					case -3:
						return	(3./8.) * TMath::Sqrt(35. / (TMath::Pi())) * TMath::Cos(theta) * TMath::Power(TMath::Sin(theta),3);
					case -2:
						return	(3./8.) * TMath::Sqrt(5. / (2 * TMath::Pi())) * TMath::Power(TMath::Sin(theta),2) * (7 * TMath::Power(TMath::Cos(theta),2) - 1);
					case -1:
						return	(3./8.) * TMath::Sqrt(5. / (TMath::Pi())) * TMath::Sin(theta) * (7 * TMath::Power(TMath::Cos(theta),3) - 3 * TMath::Cos(theta));
					case 0 :
						return	(3./16.) * TMath::Sqrt(1 / TMath::Pi()) * (35 * TMath::Power(TMath::Cos(theta),4) - 30 * TMath::Power(TMath::Cos(theta),2) + 3);
					case 1 :
						return	-(3./8.) * TMath::Sqrt(5. / (TMath::Pi())) * TMath::Sin(theta) * (7 * TMath::Power(TMath::Cos(theta),3) - 3 * TMath::Cos(theta));
					case 2 :
						return	(3./8.) * TMath::Sqrt(5. / (2 * TMath::Pi())) * TMath::Power(TMath::Sin(theta),2) * (7 * TMath::Power(TMath::Cos(theta),2) - 1);
					case 3 :
						return	-(3./8.) * TMath::Sqrt(35. / (TMath::Pi())) * TMath::Cos(theta) * TMath::Power(TMath::Sin(theta),3);
					case 4 :
						return	(3./16.) * TMath::Sqrt(35. / (2 * TMath::Pi())) * TMath::Power(TMath::Sin(theta),4);
					return 0;
				}
			case 6:
				switch (m){
					case -6:
						return	(1./64.) * TMath::Sqrt(3003/TMath::Pi()) * TMath::Power(TMath::Sin(theta),6);
					case -5:
						return	(3./32.) * TMath::Sqrt(1001/TMath::Pi()) * TMath::Power(TMath::Sin(theta),5) * TMath::Cos(theta);
					case -4:
						return	(3./32.) * TMath::Sqrt(91/(2*TMath::Pi())) * TMath::Power(TMath::Sin(theta),4) * (11. * TMath::Power(TMath::Cos(theta),2) - 1);
					case -3:
						return	(1./32.) * TMath::Sqrt(1365/TMath::Pi()) * TMath::Power(TMath::Sin(theta),3) * (11 * TMath::Power(TMath::Cos(theta),3) - 3 * TMath::Cos(theta));
					case -2:
						return	(1./64.) * TMath::Sqrt(1365/TMath::Pi()) * TMath::Power(TMath::Sin(theta),2) * (33 * TMath::Power(TMath::Cos(theta),4) - 18 * TMath::Power(TMath::Cos(theta),2) + 1);
					case -1:
						return	(1./16.) * TMath::Sqrt(273/(2*TMath::Pi())) * TMath::Sin(theta) * (33 * TMath::Power(TMath::Cos(theta),5) - 30 * TMath::Power(TMath::Cos(theta),3) + 5 * TMath::Cos(theta));
					case 0 :
						return 	(1. / 32.) * TMath::Sqrt(13. / TMath::Pi()) * (231 * TMath::Power(TMath::Cos(theta),6) - 315 * TMath::Power(TMath::Cos(theta),4) + 105 * TMath::Power(TMath::Cos(theta),2) - 5);
					case 1 :
						return	-(1./16.) * TMath::Sqrt(273/(2*TMath::Pi())) * TMath::Sin(theta) * (33 * TMath::Power(TMath::Cos(theta),5) - 30 * TMath::Power(TMath::Cos(theta),3) + 5 * TMath::Cos(theta));
					case 2 :
						return	(1./64.) * TMath::Sqrt(1365/TMath::Pi()) * TMath::Power(TMath::Sin(theta),2) * (33 * TMath::Power(TMath::Cos(theta),4) - 18 * TMath::Power(TMath::Cos(theta),2) + 1);
					case 3 :
						return	-(1./32.) * TMath::Sqrt(1365/TMath::Pi()) * TMath::Power(TMath::Sin(theta),3) * (11 * TMath::Power(TMath::Cos(theta),3) - 3 * TMath::Cos(theta));
					case 4 :
						return	(3./32.) * TMath::Sqrt(91/(2*TMath::Pi())) * TMath::Power(TMath::Sin(theta),4) * (11. * TMath::Power(TMath::Cos(theta),2) - 1);
					case 5 :
						return	-(3./32.) * TMath::Sqrt(1001/TMath::Pi()) * TMath::Power(TMath::Sin(theta),5) * TMath::Cos(theta);
					case 6 :
						return	(1./64.) * TMath::Sqrt(3003/TMath::Pi()) * TMath::Power(TMath::Sin(theta),6);
					return 0;
				}
			return 0;
		}

	}

	std::cout	<< l << "\t"
			<< m << std::endl;

	return 0;

}


double MiscFunctions::SimpsonsRule(TGraph* g, int nSteps, double xMin, double xMax){

	double xStep = (xMax - xMin)/((double)nSteps - 1.);

	double integral = 0;

	for(int x = 0; x < nSteps; x++){
		double xVal = xMin + x*xStep;

		if(x == 0 || x == (nSteps-1))
			integral += (xStep/3.)*g->Eval(xVal);
		else if((x % 2) == 0)
			integral += (xStep/3.)*(g->Eval(xVal) * 2);
		else 
			integral += (xStep/3.)*(g->Eval(xVal) * 4);
	} 

	return integral;
	
}

TGraph* MiscFunctions::PlotCrossSection(PointCoulEx* poin, double tMin, double tMax, int nTheta, int nState, bool cm, int nPart){

	//double tMinPlot = tMin;
	//double tMaxPlot = tMax;
	double tStep = (tMax - tMin)/(nTheta - 1);
	if(!cm){	// If theta is in lab units
		tMin = poin->GetReaction()->ConvertThetaLabToCm(tMin * TMath::DegToRad(), nPart) * TMath::RadToDeg(); 
		tMax = poin->GetReaction()->ConvertThetaLabToCm(tMax * TMath::DegToRad(), nPart) * TMath::RadToDeg();
		tStep = (tMax - tMin)/(nTheta - 1);
	}

	TGraph *g = new TGraph();
	for(int t = 0; t < nTheta; t++){
		poin->CalculatePointProbabilities(tMin + t*tStep);
		double theta = tMin + t*tStep;
		if(!cm)
			theta = poin->GetReaction()->ConvertThetaLabToCm(theta * TMath::DegToRad(), nPart) * TMath::RadToDeg();
		g->SetPoint(t, theta, poin->GetProbabilitiesVector()[nState]*poin->GetReaction()->RutherfordCM(tMin + t*tStep));
	}

	return g;

}

TGraph* MiscFunctions::PlotProbability(PointCoulEx* poin, double tMin, double tMax, int nTheta, int nState, bool cm, int nPart){

	//double tMinPlot = tMin;
	//double tMaxPlot = tMax;
	double tStep = (tMax - tMin)/(nTheta - 1);
	if(!cm){	// If theta is in lab units
		tMin = poin->GetReaction()->ConvertThetaLabToCm(tMin * TMath::DegToRad(), nPart) * TMath::RadToDeg(); 
		tMax = poin->GetReaction()->ConvertThetaLabToCm(tMax * TMath::DegToRad(), nPart) * TMath::RadToDeg();
		tStep = (tMax - tMin)/(nTheta - 1);
	}

	TGraph *g = new TGraph();
	for(int t = 0; t < nTheta; t++){
		poin->CalculatePointProbabilities(tMin + t*tStep);
		double theta = tMin + t*tStep;
		if(!cm)
			theta = poin->GetReaction()->ConvertThetaCmToLab(theta * TMath::DegToRad(), nPart) * TMath::RadToDeg();
		g->SetPoint(t, theta, poin->GetProbabilitiesVector()[nState]);
	}

	return g;

}
