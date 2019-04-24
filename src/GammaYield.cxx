#include "GammaYield.h"

TMatrixD GammaYield::GammaRayYield(TVectorD vec, TMatrixD mat){

	TVectorD tmpVec;
	TMatrixD tmpMat;
	tmpVec.ResizeTo(vec.GetNrows());
	tmpMat.ResizeTo(mat.GetNrows(),mat.GetNcols());
	tmpVec = vec;

	// Create effective cross section matrix: cross section for a transition 
	for(int s = tmpVec.GetNrows()-1;s>=0;s--){
		double sum=0;
		for(int sss = 0; sss<tmpVec.GetNrows(); sss++)
			sum += tmpMat[s][sss];
		for(int ss = 0; ss < tmpMat.GetNrows(); ss++){
			tmpMat[ss][s] 	= (vec[s] + sum) * mat[ss][s];
		}
	}

	return tmpMat;
	
}

void GammaYield::PrintYields(ExperimentRange r, TransitionRates t, Nucleus nucl){

	TMatrixD mat;
	mat.ResizeTo(GammaRayYield(r,t).GetNrows(),GammaRayYield(r,t).GetNcols());
	mat 	= GammaRayYield(r,t);

	std::cout	<< std::setw(15) << std::left << "Init. Index:" 
			<< std::setw(15) << std::left << "Final Index:"
			<< std::setw(10) << std::left << "Init. J:"
			<< std::setw(10) << std::left << "Final J:"
			<< std::setw(10) << std::left << "Yield:"
			<< std::endl;

	for(int x=0;x<mat.GetNcols();x++){
		for(int y=0;y<mat.GetNrows();y++){
			if(mat[y][x] > 0){
				std::cout	<< std::setw(15) << std::left << x 
						<< std::setw(15) << std::left << y
						<< std::setw(10) << std::left << nucl.GetLevelJ().at(x)
						<< std::setw(10) << std::left << nucl.GetLevelJ().at(y)
						<< std::setw(10) << std::left << mat[y][x]
						<< std::endl;
			}
		}
	}

}


double GammaYield::GetYield(ExperimentRange r, TransitionRates t, Nucleus nucl, int i, int f){

	TMatrixD mat;
	mat.ResizeTo(GammaRayYield(r,t).GetNrows(),GammaRayYield(r,t).GetNcols());
	mat 	= GammaRayYield(r,t);

	return mat[f][i];

}

