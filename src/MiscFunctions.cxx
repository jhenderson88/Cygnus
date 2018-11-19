#include "MiscFunctions.h"

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
