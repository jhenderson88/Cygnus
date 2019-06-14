#include "Integral.h"

void ThreadTaskIntegral(PointCoulEx &p, double theta){

	p.CalculatePointProbabilities(theta);	

	return;
	
}

Integral::Integral() : fNucleus(NULL), fReaction(NULL){

	fAccuracy	= 1e-5;

	fProjectileExcitation	= true;
	targetDetection = false;
	verbose		= false;

	theta_meshpoints.clear();
	energy_meshpoints.clear();
		
	point_calculations.clear();	

	nThreads 	= 1;

	fComplete	= false;

}

Integral::Integral(Nucleus *nucl, Reaction *reac) : fNucleus(nucl), fReaction(reac) {

	fProjectileExcitation	= true;
	targetDetection = false;
	verbose		= false;

	theta_meshpoints.clear();
	energy_meshpoints.clear();
		
	point_calculations.clear();

	nThreads	= 1;
	fAccuracy	= 1e-5;

	fComplete	= false;

}

Integral::Integral(const Integral& x) : 
	fNucleus(x.fNucleus), fReaction(x.fReaction), 
	point_calculations(x.point_calculations.size()), energymeshpoint_reaction(x.energymeshpoint_reaction.size()), 
	meshpointCrossSections(x.meshpointCrossSections.size()), meshpointProbabilities(x.meshpointProbabilities.size())	
{

	fProjectileExcitation	= x.fProjectileExcitation;

	nThreads	= x.nThreads;
	fAccuracy	= x.fAccuracy;

	verbose		= x.verbose;
	targetDetection	= x.targetDetection;

	fComplete	= x.fComplete;

	for(size_t i = 0; i<x.point_calculations.size(); i++)
		point_calculations.at(i) = x.point_calculations.at(i);
	for(size_t i = 0; i<x.energymeshpoint_reaction.size(); i++)
		energymeshpoint_reaction.at(i) = x.energymeshpoint_reaction.at(i);
	for(size_t i = 0; i<x.meshpointCrossSections.size(); i++){
		meshpointCrossSections.at(i).resize(x.meshpointCrossSections.at(i).size());
		for(size_t j=0;j<x.meshpointCrossSections.at(i).size();j++){
			meshpointCrossSections.at(i).at(j).ResizeTo(x.meshpointCrossSections.at(i).at(j));
			meshpointCrossSections.at(i).at(j) = x.meshpointCrossSections.at(i).at(j);
		}
	}
	for(size_t i = 0; i<x.meshpointProbabilities.size(); i++){
		meshpointProbabilities.at(i).resize(x.meshpointProbabilities.at(i).size());
		for(size_t j=0;j<x.meshpointProbabilities.at(i).size();j++){
			meshpointProbabilities.at(i).at(j).ResizeTo(x.meshpointProbabilities.at(i).at(j));
			meshpointProbabilities.at(i).at(j) = x.meshpointProbabilities.at(i).at(j);
		}
	}

}
Integral& Integral::operator = (const Integral& x){

	fProjectileExcitation	= x.fProjectileExcitation;

	nThreads	= x.nThreads;
	fAccuracy	= x.fAccuracy;

	verbose		= x.verbose;
	targetDetection	= x.targetDetection;

	fNucleus 	= x.fNucleus;
	fReaction 	= x.fReaction;

	fComplete	= x.fComplete;
	
	point_calculations.resize(x.point_calculations.size());
	energymeshpoint_reaction.resize(x.energymeshpoint_reaction.size());
	meshpointCrossSections.resize(x.meshpointCrossSections.size());
	meshpointProbabilities.resize(x.meshpointProbabilities.size());

	for(size_t i = 0; i<x.point_calculations.size(); i++)
		point_calculations.at(i) = x.point_calculations.at(i);
	for(size_t i = 0; i<x.energymeshpoint_reaction.size(); i++)
		energymeshpoint_reaction.at(i) = x.energymeshpoint_reaction.at(i);
	for(size_t i = 0; i<x.meshpointCrossSections.size(); i++){
		meshpointCrossSections.at(i).resize(x.meshpointCrossSections.at(i).size());
		for(size_t j=0;j<x.meshpointCrossSections.at(i).size();j++){
			meshpointCrossSections.at(i).at(j).ResizeTo(x.meshpointCrossSections.at(i).at(j));
			meshpointCrossSections.at(i).at(j) = x.meshpointCrossSections.at(i).at(j);
		}
	}
	for(size_t i = 0; i<x.meshpointProbabilities.size(); i++){
		meshpointProbabilities.at(i).resize(x.meshpointProbabilities.at(i).size());
		for(size_t j=0;j<x.meshpointProbabilities.at(i).size();j++){
			meshpointProbabilities.at(i).at(j).ResizeTo(x.meshpointProbabilities.at(i).at(j));
			meshpointProbabilities.at(i).at(j) = x.meshpointProbabilities.at(i).at(j);
		}
	}

	return *this;

}

void Integral::SetThetaMeshpoint(int i, double t){

	if(i<(int)theta_meshpoints.size()){
		theta_meshpoints[i] = t;
	}
	else{
		std::cout << "Energy meshpoint index " << i << " too large, appending to list" << std::endl;
		AddThetaMeshpoint(t);
	}

}
void Integral::SetEnergyMeshpoint(int i, double e){

	if(i<(int)energy_meshpoints.size()){
		energy_meshpoints[i] = e;
	}
	else{
		std::cout << "Theta meshpoint index " << i << " too large, appending to list" << std::endl;
		AddEnergyMeshpoint(e);
	}

}

void Integral::CalculateIntegral(){

	cmTheta.clear();	
	meshpointCrossSections.clear();
	meshpointProbabilities.clear();

	TVectorD out;
	out.ResizeTo(fNucleus->GetNstates());

	int part = 2;
	if(targetDetection)
		part = 3;

	std::vector<int> index;

	//	To simplify things for multithreading, here we pre-create
	//	all of the vectors and matrices containing the point-
	//	calculations.

	energymeshpoint_reaction.clear();
	std::vector< std::vector <PointCoulEx> >	PoinMat;
	std::vector< std::vector <double> >		ThetaMat;

	for(unsigned int mE = 0; mE < energy_meshpoints.size(); mE++){

		Reaction tmpReac = *fReaction;
		energymeshpoint_reaction.push_back(tmpReac);
		energymeshpoint_reaction[mE].SetLabEnergy(energy_meshpoints.at(mE));

		std::vector <PointCoulEx>	tmpPoinVec;
		std::vector <double>		tmpThetaVec;

		for(unsigned int mT = 0; mT < theta_meshpoints.size(); mT++){
			Nucleus nucl = *fNucleus;
			PointCoulEx tmpPoint(&nucl,&energymeshpoint_reaction[mE]);
			tmpPoint.SetAccuracy(fAccuracy);
			tmpPoint.SetProjectileExcitation(fProjectileExcitation);
			tmpPoinVec.push_back(tmpPoint);
			double tmpTheta = energymeshpoint_reaction[mE].ConvertThetaLabToCm(TMath::DegToRad() * theta_meshpoints.at(mT),part) * TMath::RadToDeg();
			tmpThetaVec.push_back(tmpTheta);
		}

		PoinMat.push_back(tmpPoinVec);
		ThetaMat.push_back(tmpThetaVec);

	}

	std::vector<TVectorD>	tmpVector_P;
	std::vector<TVectorD>	tmpVector_CS;	
	std::vector<double>	tmpVector_cmTheta;
	

	//	If we've specified that multiple threads will be used:
	if(nThreads > 1){

		std::vector<std::thread> Threads;
		Threads.resize(nThreads-1);

		size_t energycounter = 0;
		size_t thetacounter = 0;
		while(energycounter < theta_meshpoints.size()){
			size_t activethreads = 0;
			for(size_t t=0; t < (size_t)(nThreads-1); t++){
				Threads[t] = std::thread(ThreadTaskIntegral, std::ref(PoinMat.at(energycounter).at(thetacounter)), ThetaMat.at(energycounter).at(thetacounter));
				activethreads++;
				thetacounter++;
				if(thetacounter == theta_meshpoints.size()){
					thetacounter = 0;
					energycounter++;
					if(energycounter == energy_meshpoints.size())
						break;
				}
			}
			if(energycounter < energy_meshpoints.size()){
				ThreadTaskIntegral(std::ref(PoinMat.at(energycounter).at(thetacounter)), ThetaMat.at(energycounter).at(thetacounter));
				thetacounter++;
				if(thetacounter == theta_meshpoints.size()){
					thetacounter = 0;
					energycounter++;
				}
			}
			for(size_t t=0; t < activethreads; t++)
				Threads[t].join();
			if(energycounter == energy_meshpoints.size())
				break;
		}
		energycounter++;
	
		for(unsigned int mE = 0; mE < energy_meshpoints.size(); mE++){

			tmpVector_cmTheta.clear();
			tmpVector_CS.clear();
			tmpVector_P.clear();

			for(unsigned int mT = 0; mT < theta_meshpoints.size(); mT++){

				TVectorD tmpVec_CS;
				TVectorD tmpVec_P;
				tmpVec_P.ResizeTo(PoinMat.at(mE).at(mT).GetProbabilitiesVector().GetNrows());
				tmpVec_CS.ResizeTo(PoinMat.at(mE).at(mT).GetProbabilitiesVector().GetNrows());

				tmpVec_P = PoinMat.at(mE).at(mT).GetProbabilitiesVector();

				tmpVec_CS = PoinMat.at(mE).at(mT).GetProbabilitiesVector();
				//tmpVec_CS *= energymeshpoint_reaction.at(mE).RutherfordCM(ThetaMat.at(mE).at(mT));
		
				tmpVector_P.push_back(tmpVec_P);
				tmpVector_CS.push_back(tmpVec_CS);
				tmpVector_cmTheta.push_back(ThetaMat.at(mE).at(mT));

				index.push_back(mE);
				track_theta.push_back(theta_meshpoints.at(mT));
				point_calculations.push_back(PoinMat.at(mE).at(mT));

			}

			cmTheta.push_back(tmpVector_cmTheta);	
			meshpointCrossSections.push_back(tmpVector_CS);	
			meshpointProbabilities.push_back(tmpVector_P);
	

		}

	}
	else{
		for(unsigned int mE = 0; mE < energy_meshpoints.size(); mE++){

			tmpVector_cmTheta.clear();
			tmpVector_CS.clear();
			tmpVector_P.clear();

			for(unsigned int mT = 0; mT < theta_meshpoints.size(); mT++){
	
				PoinMat.at(mE).at(mT).CalculatePointProbabilities(ThetaMat.at(mE).at(mT));

				TVectorD tmpVec_CS;
				TVectorD tmpVec_P;
				tmpVec_P.ResizeTo(PoinMat.at(mE).at(mT).GetProbabilitiesVector().GetNrows());
				tmpVec_CS.ResizeTo(PoinMat.at(mE).at(mT).GetProbabilitiesVector().GetNrows());

				tmpVec_P = PoinMat.at(mE).at(mT).GetProbabilitiesVector();

				tmpVec_CS = PoinMat.at(mE).at(mT).GetProbabilitiesVector();
				tmpVec_CS *= energymeshpoint_reaction.at(mE).RutherfordCM(ThetaMat.at(mE).at(mT));
		
				tmpVector_P.push_back(tmpVec_P);
				tmpVector_CS.push_back(tmpVec_CS);
				tmpVector_cmTheta.push_back(ThetaMat.at(mE).at(mT));

				index.push_back(mE);
				track_theta.push_back(theta_meshpoints.at(mT));
				point_calculations.push_back(PoinMat.at(mE).at(mT));

			}

			cmTheta.push_back(tmpVector_cmTheta);	
			meshpointCrossSections.push_back(tmpVector_CS);	
			meshpointProbabilities.push_back(tmpVector_P);
	

		}
	}

	fTensors.resize(point_calculations.size());		

	fComplete = true;

}

std::vector<TVectorD> Integral::GetProbabilities(){
	
	std::vector<TVectorD>	tmpProb;
	for(unsigned int i=0;i < point_calculations.size(); i++)
		tmpProb.push_back(point_calculations.at(i).GetProbabilitiesVector());

	return tmpProb;

}

