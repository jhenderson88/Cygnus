#include "Reaction.h"

// Set up some useful nuclear constants
double Reaction::hbarc = 197.326;
double Reaction::finestruc = 0.007297352;
double Reaction::nuclearmagneton = 0.105155;
double Reaction::dipole=0.005;
double Reaction::electronCharge=1.602e-19;

Reaction::Reaction(){

	reactionSet = false;
	SetBeamA(0);
	SetTargetA(0);
	SetMass(0,0);
	SetExcitationEnergy(0);
}

Reaction::Reaction(int aB, int zB, int aT, int zT, double El){

	reactionSet = false;

	SetMass((double)aB,(double)aT);
	SetBeamA(aB);
	SetBeamZ(zB);
	SetTargetA(aT);
	SetTargetZ(zT);
	SetLabEnergy(El);
	SetExcitationEnergy(0);

	InitReaction();

}
Reaction::Reaction(double mB, int zB, double mT, int zT, double El){

	reactionSet =		false;

	SetMass(mB,mT);
	SetBeamA((int)mB);
	SetBeamZ(zB);
	SetTargetA((int)mT);
	SetTargetZ(zT);
	SetLabEnergy(El);
	SetExcitationEnergy(0);

	InitReaction();

}
Reaction::Reaction(const Reaction& r){


	beamA		= r.beamA;
	beamZ		= r.beamZ;
	targetA		= r.targetA;
	targetZ		= r.targetZ;

	beamMass	= r.beamMass;
	targetMass	= r.targetMass;
	
	Elab		= r.Elab;
	Ecm		= r.Ecm;
	beta		= r.beta;

	exE		= r.exE;

	InitReaction();

}
Reaction& Reaction::operator = (const Reaction& r){

	beamA		= r.beamA;
	beamZ		= r.beamZ;
	targetA		= r.targetA;
	targetZ		= r.targetZ;

	beamMass	= r.beamMass;
	targetMass	= r.targetMass;
	
	Elab		= r.Elab;
	Ecm		= r.Ecm;
	beta		= r.beta;

	exE		= r.exE;

	InitReaction();

	return *this;	

}


void Reaction::SetLabEnergy(double E){
	reactionSet =		false;
	Elab = E;
	if(GetBeamMass() > 0 && GetTargetMass() > 0)
	{
		Ecm = GetLabEnergy()/(1+(GetBeamMass()/GetTargetMass()));
		beta = 0.046337*TMath::Sqrt(GetLabEnergy() / GetBeamMass());
	}
	InitReaction();
}
	
void Reaction::SetCMEnergy(double E){
	reactionSet =		false;
	Ecm = E; 
	if(GetBeamMass() > 0 && GetTargetMass() > 0)
	{
		Elab = GetCMEnergy()*(1+(GetBeamMass()/GetTargetMass()));
		beta = 0.046337*TMath::Sqrt(GetLabEnergy() / GetBeamMass());
	}
	InitReaction();
}

double Reaction::ThreeJ(int I1, int I2, int I3, int M1, int M2, int M3) const{

	double threej = (double)ROOT::Math::wigner_3j(I1,I2,I3,M1,M2,M3); // Uses the GNU science libraries

	return threej;

}

double Reaction::ClosestApproach(int PA, int TA, int PZ, int TZ, double E) const{

	double A = (2. * 0.71999 * (1. + ((double)PA)/((double)TA)) * ((double)PZ) * ((double)TZ)) / E;

	return A;

}

double Reaction::EtaCalc(double E) const{

	double b = (1 + (double)GetBeamA()/(double)GetTargetA());

	double s = ((double)GetBeamZ() * (double)GetTargetZ() * TMath::Sqrt((double)GetBeamA())/6.34977) * (1/TMath::Sqrt(GetLabEnergy() - b*E));

	return s;

}

double Reaction::RelativeVelocity(double State_E) const{

	double v = 0.046337*TMath::Sqrt((GetLabEnergy() - State_E*(1 + ((double)GetBeamA()/(double)GetTargetA()) ) )/ ((double)GetBeamMass()));

	return v;

}

double Reaction::RutherfordCM(double theta){
	
	if(!reactionSet)
		InitReaction();
	
	theta *= TMath::DegToRad();

	return TMath::Power( TMath::Power(electronCharge,2) * GetBeamZ() * GetTargetZ() / (16 * TMath::Pi() * 8.85e-12 * Ecm * 1.60e-13) ,2) * 1 / (TMath::Power( TMath::Sin( theta / 2) , 4)) * TMath::Sin(theta) * 1E24;

}

double Reaction::Rutherford(double theta_lab, int part){
	
	if(!reactionSet)
		InitReaction();

	double theta = ConvertThetaLabToCm(theta_lab * TMath::DegToRad(), part);
	
	return TMath::Power( TMath::Power(electronCharge,2) * GetBeamZ() * GetTargetZ() / (16 * TMath::Pi() * 8.85e-12 * Ecm * 1.60e-13) ,2) * 1 / (TMath::Power( TMath::Sin( theta / 2) , 4)) * TMath::Sin(theta) * 1E24;

}

void Reaction::InitReaction(){

	fM[0]	=		GetBeamMassMeV();
	fM[1]	=		GetTargetMassMeV();
	fM[2]	=		fM[0];
	fM[3]	=		fM[1];

	fTLab[0] =	 	GetLabEnergy();									// Kinetic energy of the beam
	fELab[0] =	 	fTLab[0] + fM[0];								// Total energy of the beam (lab)
	fPLab[0] =	 	TMath::Sqrt(TMath::Power(fTLab[0],2) + 2 * fTLab[0] * fM[0]);			// Momentum of the beam (lab)
	fVLab[0] =		fPLab[0] / fELab[0]; 								// Velocity of the beam (lab)
	fGLab[0] =	 	1 / TMath::Sqrt(1 - TMath::Power(fVLab[0],2));					// Gamma factor of the beam

	fTLab[1] =		0;										// Kinetic energy of the target (0)
	fELab[1] =		fM[1];										// Total energy of the target (lab)
	fPLab[1] = 		0;										// Momentum of the target (lab)
	fVLab[1] = 		0;										// Velocity of the target (lab)
	fGLab[1] =		1;										// Gamma factor of the target	

	fS = 			TMath::Power(fM[0],2) + TMath::Power(fM[1],2) + 2 * fELab[0] * fM[1];
	fInvariantMass =	TMath::Sqrt(fS);	

	// CM motion
	
	fCmE 	=		fInvariantMass;
	fCmTi 	= 		fCmE - fM[0] - fM[1];
	fCmTf	=		fCmTi + fQVal;
	fCmV 	=		fPLab[0] / (fELab[0] + fM[1]);
	fCmP 	=		fCmV * fCmE;
	fCmG 	=		1 / TMath::Sqrt(1 - TMath::Power(fCmV,2));

	SetCmFrame(exE);

	reactionSet =		true;
	
}

void Reaction::SetCmFrame(double exc){

	// Particles in CM frame:
	fPCm[0]	=		TMath::Sqrt((fS - TMath::Power(fM[0] - fM[1],2)) * (fS - TMath::Power(fM[0] + fM[1],2)))/(2 * TMath::Sqrt(fS));
	fPCm[1]	=		fPCm[0];
	fPCm[2]	=		TMath::Sqrt((fS - TMath::Power(exc + fM[2] - fM[3],2)) * (fS - TMath::Power(exc + fM[2] + fM[3],2)))/(2 * TMath::Sqrt(fS));	
	fPCm[3]	=		fPCm[2];

	for(int i=0;i<4;i++){

		fECm[i]	=	TMath::Sqrt(TMath::Power(fM[i],2) + TMath::Power(fPCm[i],2));
		fTCm[i]	=	fECm[i] - fM[i];
		fVCm[i]	=	fPCm[i] / fECm[i];
		fGCm[i]	=	1 / TMath::Sqrt(1-TMath::Power(fVCm[i],2));

		if( i < 2)
			fThetaMax[i] =	0;
		else{
			double val = fPCm[i] / (fM[i]*fCmV*fCmG);
			if(val < 1)
				fThetaMax[i] = TMath::ASin(val);
			else if(val < 1.001)
				fThetaMax[i] = TMath::Pi()/2;
			else
				fThetaMax[i] = TMath::Pi();
		}

	}

}

double Reaction::ConvertThetaLabToCm(double theta_lab, int part) const{

	if(theta_lab > fThetaMax[part])
		theta_lab = fThetaMax[part];

	double gtan2 	=	TMath::Power(TMath::Tan(theta_lab)*fCmG,2);
	double x 	=	fCmV / fVCm[part];
	double expr	=	TMath::Sqrt(1 + gtan2 * (1 - TMath::Power(x,2)));
	double theta_cm;

	if(TMath::Tan(theta_lab) >= 0)
		theta_cm = TMath::ACos((-x * gtan2 + expr)/(1 + gtan2));
	else
		theta_cm = TMath::ACos((-x * gtan2 - expr)/(1 + gtan2));
		

	if(part == 3)
		theta_cm = TMath::Pi() - theta_cm;

	return theta_cm;

}
double Reaction::ConvertThetaCmToLab(double theta_cm, int part) const{

	if(part == 3)
		theta_cm = TMath::Pi() - theta_cm;
	
	double theta_lab = TMath::ATan2(TMath::Sin(theta_cm),fCmG * ( TMath::Cos(theta_cm) + fCmV/fVCm[part]));

	if(theta_lab > fThetaMax[part])
		return fThetaMax[part];

	return theta_lab;

}

TGraph* Reaction::ThetaVsTheta(double thmin, double thmax, int part) const{

	TGraph *g = new TGraph();
	
	g->SetName("ThetaVsTheta_Lab_CM");
	g->SetTitle("Angular Conversion Lab to CM");

	double theta_cm, theta_lab;

	int counter = 0;

	for(int i=0;i<180;i++){

		theta_lab = (double)i;
		theta_cm = ConvertThetaLabToCm(theta_lab * TMath::DegToRad(),part) * TMath::RadToDeg();

		if(theta_lab < thmin || theta_lab > thmax)
			continue;
	
		g->SetPoint(counter,theta_lab,theta_cm);

		counter++;

	}

	return g;

}

TGraph*	Reaction::RutherfordGraph(double thmin, double thmax, int part){

	TGraph *g = new TGraph();
		
	g->SetName("Rutherford_CS");
	g->SetTitle("Rutherford cross section in lab");
	
	double theta_lab, rCS;

	int counter = 0;

	for(int i=0;i<180;i++){
		
		theta_lab = (double)i;
		
		if(theta_lab < thmin || theta_lab > thmax)
			continue;
	
		rCS = Rutherford(theta_lab,part);
	
		g->SetPoint(counter,theta_lab,rCS);

		counter++;

	}

	return g;

} 

void Reaction::PrintReaction() const{

	std::cout << "Beam:\tA = " << GetBeamA() << "\t Z = " << GetBeamZ() << std::endl;
	std::cout << "Target:\tA = " << GetTargetA() << "\t Z = " << GetTargetZ() << std::endl;
	std::cout << "Energy: " << GetLabEnergy() << std::endl;

	std::cout << "Beam mass: " << fM[0] << " target mass: " << fM[1] << " [MeV]" << std::endl;
	std::cout << "Beam mass: " << GetBeamMass() << " target mass: " << GetTargetMass() << " [amu]" << std::endl;

	std::cout << "Beam velocity: " << fVLab[0] << std::endl;

}
