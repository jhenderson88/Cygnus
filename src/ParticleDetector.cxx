#include "ParticleDetector.h"

ParticleDetector::ParticleDetector(const ParticleDetector& d){

	gThetaEff	= d.gThetaEff;
	gThetaPhi	= d.gThetaPhi;
	setFlag		= d.setFlag;
	
}
ParticleDetector& ParticleDetector::operator = (const ParticleDetector& d){

	gThetaEff	= d.gThetaEff;
	gThetaPhi	= d.gThetaPhi;
	setFlag		= d.setFlag;

	return *this;

}

ParticleDetector::ParticleDetector(TGraph2D* g2, TGraph* g) : gThetaPhi(g2), gThetaEff(g){

	setFlag = true;	

} 

double ParticleDetector::GetEfficiencyTheta(double t){

	double eff = 1;

	if(setFlag)
		eff = gThetaEff->Eval(t);
	else
		std::cout << "Warning: Theta phi map not defined. Returning efficiency = 1" << std::endl;

	return eff;

}

double ParticleDetector::GetThetaMin(){

	double theta = 180;
	if(setFlag){
		double x,y;
		double step = 0.001;
		//for(int i=0;i<gThetaEff->GetN();i++){
		//	gThetaEff->GetPoint(i,x,y);
		for(int i=0;i<180000;i++){
			x = i * step;
			y = gThetaEff->Eval(x);
			if(y > 0 && x < theta)
				theta = x;
		}
		return theta;
	}
	theta = 0;
	std::cout << "Warning: Theta phi map not well defined. Theta min set to 0 degrees" << std::endl;

	return theta;

}
double ParticleDetector::GetThetaMax(){
	
	double theta = 0;
	if(setFlag){
		double x,y;
		double step = 0.001;
		//for(int i=0;i<gThetaEff->GetN();i++){
		//	gThetaEff->GetPoint(i,x,y);
		for(int i=0;i<180000;i++){
			x = i * step;
			y = gThetaEff->Eval(x);
			if(y > 0 && x > theta)
				theta = x;
		}
		return theta;
	}
	theta = 180;
	std::cout << "Warning: Theta phi map not well defined. Theta max set to 180 degrees" << std::endl;

	return theta;

}

void ParticleDetector::WriteParticleDetector(const char* outfilename, const char* opt){

	TFile *outfile = new TFile(outfilename,opt);
	TDirectory *dir;
	if(outfile->GetDirectory("ParticleDetectors"))
		dir = outfile->GetDirectory("ParticleDetectors");
	else
		dir = outfile->mkdir("ParticleDetectors");
	dir->cd();

	gThetaPhi->Write();
	gThetaEff->Write();
	outfile->Close();
	
}
