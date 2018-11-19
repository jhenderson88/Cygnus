#include "ParticleDetector.h"

ParticleDetector::ParticleDetector(const ParticleDetector& d){

	hThetaEff	= d.hThetaEff;
	hThetaPhi	= d.hThetaPhi;
	setFlag		= d.setFlag;
	
}
ParticleDetector& ParticleDetector::operator = (const ParticleDetector& d){

	hThetaEff	= d.hThetaEff;
	hThetaPhi	= d.hThetaPhi;
	setFlag		= d.setFlag;

	return *this;

}

double ParticleDetector::GetEfficiencyTheta(double t){

	double eff = 1;

	if(setFlag)
		eff = hThetaEff->GetBinContent(hThetaEff->GetXaxis()->FindBin(t));
	else
		std::cout << "Warning: Theta phi map not defined. Returning efficiency = 1" << std::endl;

	return eff;

}

double ParticleDetector::GetThetaMin(){

	double theta = 0;
	if(setFlag){
		for(int i=1;i<hThetaEff->GetNbinsX();i++)
			if(hThetaEff->GetBinContent(i) > 0)
				return hThetaEff->GetBinCenter(i);
	}
	std::cout << "Warning: Theta phi map not well defined. Theta min set to 0 degrees" << std::endl;

	return theta;

}
double ParticleDetector::GetThetaMax(){
	
	double theta = 180;
	if(setFlag){
		for(int i=hThetaEff->GetNbinsX();i>0;i--)
			if(hThetaEff->GetBinContent(i) > 0)
				return hThetaEff->GetBinCenter(i);
	}
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

	hThetaPhi->Write();
	hThetaEff->Write();
	outfile->Close();
	
}
