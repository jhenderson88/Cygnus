#include "ExperimentRange.h"

ExperimentRange::ExperimentRange() : fNucleus(NULL), fReaction(NULL), fIntegral(NULL), fDetectorEff(NULL)
{

	fAccuracy	= 1e-5;

	fProjectileExcitation	= true;
	fUseEfficiency	= false;
	targetDetection = false;
	verbose		= false;
	nThreads	= 1;

	fUseFit		= false;

	integratedRutherford 	= 0;

}

ExperimentRange::ExperimentRange(Nucleus* nucl, Reaction* reac) : fNucleus(nucl), fReaction(reac), fDetectorEff(NULL)
{

	fAccuracy	= 1e-5;	

	fProjectileExcitation	= true;
	fUseEfficiency	= false;
	targetDetection = false;
	verbose		= false;

	fIntegral = new Integral(nucl,reac);

	nThreads	= 1;

	fUseFit		= false;

	integratedRutherford 	= 0;

}

ExperimentRange::ExperimentRange(const ExperimentRange& e) : fNucleus(e.fNucleus), fReaction(e.fReaction), fIntegral(e.fIntegral)
{

	IntegratedCrossSection_TVec.ResizeTo(e.IntegratedCrossSection_TVec.GetNrows());
	IntegratedCrossSection_TVec = e.IntegratedCrossSection_TVec;

	fProjectileExcitation	= e.fProjectileExcitation;

	fUseEfficiency	= e.fUseEfficiency;

	fDetectorEff	= e.fDetectorEff;

	fStopping	= e.fStopping;
	
	fAccuracy	= e.fAccuracy;

	thetamin 	= e.thetamin;
	thetamax	= e.thetamax;
	energymin	= e.energymin;
	energymax	= e.energymax;
	nTheta		= e.nTheta;
	nEnergy		= e.nEnergy;
	meanThetaCM	= e.meanThetaCM;
	meanThetaLab	= e.meanThetaLab;
	meanEnergy	= e.meanEnergy;
	thetaminCM 	= e.thetaminCM;
	thetamaxCM	= e.thetamaxCM;
	targetDetection = e.targetDetection;
	verbose		= e.verbose;

	nThreads	= e.nThreads;

	fUseFit		= e.fUseFit;

	integratedRutherford	= e.integratedRutherford;

}
ExperimentRange& ExperimentRange::operator = (const ExperimentRange &e){

	fProjectileExcitation	= e.fProjectileExcitation;

	fUseEfficiency	= e.fUseEfficiency;

	fNucleus 	= e.fNucleus;
	fReaction 	= e.fReaction;
	fIntegral	= e.fIntegral;

	IntegratedCrossSection_TVec.ResizeTo(e.IntegratedCrossSection_TVec.GetNrows());
	IntegratedCrossSection_TVec = e.IntegratedCrossSection_TVec;

	fDetectorEff	= e.fDetectorEff;

	fStopping	= e.fStopping;
	
	fAccuracy	= e.fAccuracy;

	thetamin 	= e.thetamin;
	thetamax	= e.thetamax;
	energymin	= e.energymin;
	energymax	= e.energymax;
	nTheta		= e.nTheta;
	nEnergy		= e.nEnergy;
	meanThetaCM	= e.meanThetaCM;
	meanThetaLab	= e.meanThetaLab;
	meanEnergy	= e.meanEnergy;
	thetaminCM 	= e.thetaminCM;
	thetamaxCM	= e.thetamaxCM;
	targetDetection = e.targetDetection;
	verbose		= e.verbose;

	nThreads	= e.nThreads;

	fUseFit		= e.fUseFit;

	integratedRutherford	= e.integratedRutherford;

	return *this;
	
}

void ExperimentRange::IntegrateRange()
{

	fIntegral->SetAccuracy(fAccuracy);

	fIntegral->SetNthreads(nThreads);
	
	fIntegral->SetProjectileExcitation(fProjectileExcitation);
	fIntegral->SetTargetDetection(targetDetection);

	fIntegral->ClearThetaMeshpoints();
	fIntegral->ClearEnergyMeshpoints();

	double thetastep; 
	if(nTheta == 1)
		thetastep = 0;
	else
		thetastep = (thetamax - thetamin) / (nTheta - 1);
	double energystep;
	if(nEnergy == 1)
		energystep = 0;
	else
		energystep = (energymax - energymin) / (nEnergy - 1);

	double meantheta = (thetamax - thetamin)/2. + thetamin;
	double meanenergy = (energymax - energymin)/2. + energymin;

	
	int part = 2;
	if(targetDetection)
		part = 3;
	double meantheta_cm = fReaction->ConvertThetaLabToCm(meantheta * TMath::DegToRad(),part) * TMath::RadToDeg();
	meanThetaCM = meantheta_cm;
	meanThetaLab = meantheta;
	meanEnergy  = meanenergy;

	for(int n=0; n<nTheta; n++)
		fIntegral->AddThetaMeshpoint(thetamin + thetastep*n);
	for(int n=0; n<nEnergy; n++)
		fIntegral->AddEnergyMeshpoint(energymin + energystep*n);

	typedef std::chrono::high_resolution_clock Clock;
	typedef std::chrono::milliseconds milliseconds;
	Clock::time_point t0 = Clock::now();

	fIntegral->CalculateIntegral();			

	Clock::time_point t1 = Clock::now();
	milliseconds ms = std::chrono::duration_cast<milliseconds>(t1-t0);

	double thetacm_min = fIntegral->GetCMThetaPoints().at(0).at(0);
	double thetacm_max = fIntegral->GetCMThetaPoints().at(0).at(fIntegral->GetCMThetaPoints().at(0).size()-1);

	if(thetacm_min > thetacm_max){
		double tmp = thetacm_min;
		thetacm_min = thetacm_max;
		thetacm_max = tmp;
	}

	int nPart = 2;
	if(targetDetection)
		nPart = 3;
	std::cout 	<< std::setw(18) << std::left << thetacm_min 
			<< std::setw(18) << std::left << thetacm_max
		 	<< std::setw(18) << std::left << fReaction->ConvertThetaCmToLab(thetacm_min*TMath::DegToRad(),nPart)*TMath::RadToDeg()
			<< std::setw(18) << std::left << fReaction->ConvertThetaCmToLab(thetacm_max*TMath::DegToRad(),nPart)*TMath::RadToDeg()
			<< std::setw(14) << std::left << nThreads
			<< std::setw(20) << std::left << ms.count()
			<< std::endl;
	
	IntegratedCrossSection_TVec.ResizeTo(fNucleus->GetNstates());	

	for(int s=0;s<fNucleus->GetNstates();s++)
		IntegratedCrossSection_TVec[s] = IntegrateThetaEnergy(s); 

}

double ExperimentRange::IntegrateThetaEnergy(int state){

	double CS = 0;

	TGraph2D *gEvsTheta = InterpolatedEnergyTheta(state,fUseEfficiency);

	int eSteps = 101;
	int tSteps = 201;
	double eMin = fIntegral->GetEnergy(0);
	double eMax = fIntegral->GetEnergy(nEnergy-1);
	if(eMin > eMax){
		double tmp = eMin;
		eMin = eMax;
		eMax = tmp;
	}
	int nPart = 2;
	if(targetDetection)
		nPart = 3;
	double tMin = fReaction->ConvertThetaCmToLab(fIntegral->GetCMThetaPoints().at(0).at(0) * TMath::DegToRad(), nPart) * TMath::RadToDeg();
	double tMax = fReaction->ConvertThetaCmToLab(fIntegral->GetCMThetaPoints().at(0).at(fIntegral->GetCMThetaPoints().at(0).size()-1) * TMath::DegToRad(), nPart) * TMath::RadToDeg();


	if(fDetectorEff && fUseEfficiency){
		fDetectorEff_CM = new TGraph();
		double x,y;
		double tMin_tmp = 180;
		double tMax_tmp = 0;
		for(int n=0;n<fDetectorEff->GetN();n++){
			fDetectorEff->GetPoint(n,x,y);
			fDetectorEff_CM->SetPoint(n,fReaction->ConvertThetaLabToCm(x*TMath::DegToRad(),nPart)*TMath::RadToDeg(),y);
			if(y > 0 && fReaction->ConvertThetaLabToCm(x*TMath::DegToRad(),nPart)*TMath::RadToDeg() < tMin_tmp)
				tMin_tmp = fReaction->ConvertThetaLabToCm(x*TMath::DegToRad(),nPart)*TMath::RadToDeg();
			if(y > 0 && fReaction->ConvertThetaLabToCm(x*TMath::DegToRad(),nPart)*TMath::RadToDeg() > tMax_tmp)
				tMax_tmp = fReaction->ConvertThetaLabToCm(x*TMath::DegToRad(),nPart)*TMath::RadToDeg();
		}

		tMin = tMin_tmp;
		tMax = tMax_tmp;

	}

	if(tMin > tMax){
		double tmp = tMin;
		tMin = tMax;
		tMax = tmp;
	}
	double tStep = (tMax - tMin) / ((double)tSteps - 1.); // Theta stepsize
	double eStep = (eMax - eMin) / ((double)eSteps - 1.); // Energy stepsize
	
	// We can use Simpson's rule to integrate theta and energy and
	// for every energy, divide out dE/dX
	for(int e = 0; e < eSteps; e++){
		double tmpCS = 0;
		double energy = eMin + e*eStep;
		for(int t = 0; t < tSteps; t++){
			if(t==0 || t==(tSteps-1))
				tmpCS += gEvsTheta->Interpolate(tMin + t*tStep,energy);
			else if((t % 2) == 0)
				tmpCS += gEvsTheta->Interpolate(tMin + t*tStep,energy) * 2;
			else
				tmpCS += gEvsTheta->Interpolate(tMin + t*tStep,energy) * 4;
		}
		tmpCS *= (tStep / 3.);
		tmpCS /= fStopping.GetStoppingFit().Eval(energy);

		if(e == 0 || e == (eSteps-1))
			CS += (eStep/3.)*tmpCS;
		else if((e % 2) == 0)
			CS += (eStep/3.)*(tmpCS * 2);
		else 
			CS += (eStep/3.)*(tmpCS * 4);
	} 

	return CS;

}

TGraph2D* ExperimentRange::InterpolatedEnergyTheta(int state, bool useDetector){

	int eSteps = 101;
	int tSteps = 201;
	double eMin = fIntegral->GetEnergy(0);
	double eMax = fIntegral->GetEnergy(nEnergy-1);
	double tMin = fIntegral->GetCMThetaPoints().at(0).at(0);
	double tMax = fIntegral->GetCMThetaPoints().at(0).at(fIntegral->GetCMThetaPoints().at(0).size()-1);

	int nPart = 2;
	if(targetDetection)
		nPart = 3;

	tMin = fReaction->ConvertThetaCmToLab(tMin * TMath::DegToRad(),nPart) * TMath::RadToDeg();
	tMax = fReaction->ConvertThetaCmToLab(tMax * TMath::DegToRad(),nPart) * TMath::RadToDeg();

	if(tMin > tMax)
		std::swap(tMin,tMax);

	TF1 *tFits[fIntegral->GetCMThetaPoints().size()];
	TSpline3 *tSplines[fIntegral->GetCMThetaPoints().size()];
	for(int e = 0; e < nEnergy; e++){ 	// For each energy meshpoint, do a fit to the theta distribution
		TGraph gTmp;
		double *tmp_x = new double[fIntegral->GetCMThetaPoints().at(e).size()];
		double *tmp_y = new double[fIntegral->GetCMThetaPoints().at(e).size()];
		for(size_t t = 0; t < fIntegral->GetCMThetaPoints().at(e).size(); t++){
			double tmp_theta 	= fReaction->ConvertThetaCmToLab(fIntegral->GetCMThetaPoints().at(e).at(t) * TMath::DegToRad(), nPart) * TMath::RadToDeg();
			double tmp_cs 		= fIntegral->GetMeshPointProbabilities().at(e).at(t)[state] * TMath::Sin(tmp_theta * TMath::DegToRad()) * fReaction->Rutherford(tmp_theta,nPart);
			tmp_cs *= TMath::DegToRad() * 2 * TMath::Pi() / (4 * TMath::Pi());
			gTmp.SetPoint((int)t,tmp_theta,tmp_cs);
			tmp_x[t] = tmp_theta;
			tmp_y[t] = tmp_cs;
		}
		char fname[32];
		sprintf(fname,"Energy_%i",(int)e+1);
		tSplines[e] = new TSpline3(fname,tmp_x,tmp_y,fIntegral->GetCMThetaPoints().at(e).size());
		if(UseFit()){
			char fittype[16];
			if(thetamaxCM - thetaminCM > 90.)
				sprintf(fittype,"pol6");
			else
				sprintf(fittype,"pol5");
			tFits[e] = new TF1(fname,fittype,tMin,tMax);
			gTmp.Fit(tFits[e],"RQ0");
		}
	}

	if(false){
	//if(fDetectorEff && fUseEfficiency){
		fDetectorEff_CM = new TGraph();
		double x,y;
		double tMin_tmp = 180;
		double tMax_tmp = 0;
		for(int n=0;n<fDetectorEff->GetN();n++){
			fDetectorEff->GetPoint(n,x,y);
			fDetectorEff_CM->SetPoint(n,fReaction->ConvertThetaLabToCm(x*TMath::DegToRad(),nPart)*TMath::RadToDeg(),y);
			if(y > 0 && fReaction->ConvertThetaLabToCm(x*TMath::DegToRad(),nPart)*TMath::RadToDeg() < tMin_tmp)
				tMin_tmp = fReaction->ConvertThetaLabToCm(x*TMath::DegToRad(),nPart)*TMath::RadToDeg();
			if(y > 0 && fReaction->ConvertThetaLabToCm(x*TMath::DegToRad(),nPart)*TMath::RadToDeg() > tMax_tmp)
				tMax_tmp = fReaction->ConvertThetaLabToCm(x*TMath::DegToRad(),nPart)*TMath::RadToDeg();
		}

		tMin = tMin_tmp;
		tMax = tMax_tmp;

	}

	// Now we can use the fitted curve as an interpolation to create a much finer theta-energy surface before using Simpson's rule to integrate
	TF1 *eFits[tSteps];
	TSpline3 *eSplines[tSteps];
	double tStep = (tMax - tMin) / ((double)tSteps - 1.); // Theta stepsize
	double eStep = (eMax - eMin) / ((double)eSteps - 1.); // Energy stepsize
	for(int t = 0; t < tSteps; t++){
		TGraph gTmp;
		double *tmp_x = new double[nEnergy];
		double *tmp_y = new double[nEnergy];
		for(int e = 0; e < nEnergy; e++){
			if(UseFit()){
				gTmp.SetPoint(e,fIntegral->GetEnergy(e),tFits[e]->Eval(tMin + t*tStep));
			}
			else{
				gTmp.SetPoint(e,fIntegral->GetEnergy(e),tSplines[e]->Eval(tMin + t*tStep));
				tmp_x[e] = fIntegral->GetEnergy(e);
				tmp_y[e] = tSplines[e]->Eval(tMin + t*tStep);
			}
		}
		char fname[32];
		sprintf(fname,"Theta_%i",t);
		eSplines[t] = new TSpline3(fname,tmp_x,tmp_y,nEnergy);
		if(UseFit()){
			eFits[t] = new TF1(fname,"pol3",eMin,eMax);
			gTmp.Fit(eFits[t],"RQ0");
		}	
	}

	TGraph2D *g = new TGraph2D;
	int pointCounter = 0;
	for(int t = 0; t < tSteps; t++){
		double eff = 1;
		if(fDetectorEff && fUseEfficiency){
			eff = fDetectorEff_CM->Eval(tMin + t*tStep);
		}
		for(int e = 0; e < eSteps; e++){
			if(UseFit())
				g->SetPoint(pointCounter,tMin + t*tStep, eMin + e*eStep, eFits[t]->Eval(eMin + e*eStep) * eff);
			else
				g->SetPoint(pointCounter,tMin + t*tStep, eMin + e*eStep, eSplines[t]->Eval(eMin + e*eStep) * eff);
			pointCounter++;
		}
	}

	return g;

}

double ExperimentRange::IntegrateRutherford(){

	double CS = 0;

	//TGraph2D *g = GetRutherfordThetaEnergy();
	Reaction tmpReac = *fReaction;

	int eSteps = 101;
	int tSteps = 101;
	double eMin;
	double eMax;
	double tMin;
	double tMax;

	int nPart = 2;
	if(targetDetection)
		nPart = 3;
	// If the integration hasn't been performed, we can still integrate the
	// Rutherford cross section using the user-provided ranges and efficiency.
	// If the integration has been performed, we use those angles/energies 
	// to ensure consistency.
	if(fIntegral->IntegralComplete()){
		eMin = fIntegral->GetEnergy(0);
		eMax = fIntegral->GetEnergy(nEnergy-1);
		if(eMin > eMax){
			double tmp = eMin;
			eMin = eMax;
			eMax = tmp;
		}
		tMin = tmpReac.ConvertThetaCmToLab(fIntegral->GetCMThetaPoints().at(0).at(0) * TMath::DegToRad(),nPart) * TMath::RadToDeg();
		tMax = tmpReac.ConvertThetaCmToLab(fIntegral->GetCMThetaPoints().at(0).at(fIntegral->GetCMThetaPoints().at(0).size()-1) * TMath::DegToRad(),nPart) * TMath::RadToDeg();
		if(tMin > tMax){
			double tmp = tMin;
			tMin = tMax;
			tMax = tmp;
		}
	}
	else{
		eMin = energymin;
		eMax = energymax;
		if(eMin > eMax){
			double tmp = eMin;
			eMin = eMax;
			eMax = tmp;
		}
		tMin = thetamin;//fReaction->ConvertThetaLabToCm(thetamin,nPart);
		tMax = thetamax;//fReaction->ConvertThetaLabToCm(thetamax,nPart);
		if(tMin > tMax){
			double tmp = tMin;
			tMin = tMax;
			tMax = tmp;
		}
	}

	if(fDetectorEff && fUseEfficiency){
		fDetectorEff_CM = new TGraph();
		double x,y;
		double tMin_tmp = 180;
		double tMax_tmp = 0;
		for(int n=0;n<fDetectorEff->GetN();n++){
			fDetectorEff->GetPoint(n,x,y);
			fDetectorEff_CM->SetPoint(n,fReaction->ConvertThetaLabToCm(x*TMath::DegToRad(),nPart)*TMath::RadToDeg(),y);
			if(y > 0 && fReaction->ConvertThetaLabToCm(x*TMath::DegToRad(),nPart)*TMath::RadToDeg() < tMin_tmp)
				tMin_tmp = fReaction->ConvertThetaLabToCm(x*TMath::DegToRad(),nPart)*TMath::RadToDeg();
			if(y > 0 && fReaction->ConvertThetaLabToCm(x*TMath::DegToRad(),nPart)*TMath::RadToDeg() > tMax_tmp)
				tMax_tmp = fReaction->ConvertThetaLabToCm(x*TMath::DegToRad(),nPart)*TMath::RadToDeg();
		}

		tMin = tMin_tmp;
		tMax = tMax_tmp;

	}

	if(tMin > tMax){
		double tmp = tMin;
		tMin = tMax;
		tMax = tmp;
	}
	double tStep = (tMax - tMin) / ((double)tSteps - 1.); // Theta stepsize
	double eStep = (eMax - eMin) / ((double)eSteps - 1.); // Energy stepsize

	// We can use Simpson's rule to integrate theta and energy and
	// for every energy, divide out dE/dX
	for(int e = 0; e < eSteps; e++){
		double tmpCS = 0;
		double energy = eMin + e*eStep;
		tmpReac.SetLabEnergy(energy);
		for(int t = 0; t < tSteps; t++){
			double eff = 1;
			if(fDetectorEff && fUseEfficiency){
				eff = fDetectorEff_CM->Eval(tMin + t*tStep);
			}
			if(t==0 || t==(tSteps-1)){
				tmpCS += tmpReac.Rutherford(tMin + t*tStep, nPart) * eff * TMath::Sin((tMin + t*tStep) * TMath::DegToRad()) * TMath::DegToRad() * 2 * TMath::Pi();
			}
			else if((t % 2) == 0){
				tmpCS += tmpReac.Rutherford(tMin + t*tStep, nPart) * eff * 2 * TMath::Sin((tMin + t*tStep) * TMath::DegToRad()) * TMath::DegToRad() * 2 * TMath::Pi();
			}
			else{
				tmpCS += tmpReac.Rutherford(tMin + t*tStep, nPart) * eff * 4 * TMath::Sin((tMin + t*tStep) * TMath::DegToRad()) * TMath::DegToRad() * 2 * TMath::Pi();
			}
		}
		tmpCS *= (tStep / 3.);
		tmpCS /= fStopping.GetStoppingFit().Eval(energy);

		if(e == 0 || e == (eSteps-1))
			CS += (eStep/3.) * tmpCS;
		else if((e % 2) == 0)
			CS += (eStep/3.) * tmpCS * 2;
		else 
			CS += (eStep/3.) * tmpCS * 4;
	} 

	integratedRutherford = CS;

	return CS;

}

TGraph2D* ExperimentRange::GetRutherfordThetaEnergy(){

	TGraph2D *g = new TGraph2D();
	
	Reaction tmpReac = *fReaction;

	int eSteps = 101;
	int tSteps = 101;
	double eMin;
	double eMax;
	double eStep;
	double tMin;
	double tMax;
	double tStep;

	int nPart = 2;
	if(targetDetection)
		nPart = 3;

	// If the integration hasn't been performed, we can still integrate the
	// Rutherford cross section using the user-provided ranges and efficiency.
	// If the integration has been performed, we use those angles/energies 
	// to ensure consistency.
	if(fIntegral->IntegralComplete()){
		eMin = fIntegral->GetEnergy(0);
		eMax = fIntegral->GetEnergy(nEnergy-1);
		if(eMin > eMax){
			double tmp = eMin;
			eMin = eMax;
			eMax = tmp;
		}
		tMin = tmpReac.ConvertThetaCmToLab(fIntegral->GetCMThetaPoints().at(0).at(0) * TMath::DegToRad(), nPart)*TMath::RadToDeg();
		tMax = tmpReac.ConvertThetaCmToLab(fIntegral->GetCMThetaPoints().at(0).at(fIntegral->GetCMThetaPoints().at(0).size()-1) * TMath::DegToRad(), nPart) * TMath::RadToDeg();
		if(tMin > tMax){
			double tmp = tMin;
			tMin = tMax;
			tMax = tmp;
		}
		tStep = (tMax - tMin) / ((double)tSteps - 1.); // Theta stepsize
		eStep = (eMax - eMin) / ((double)eSteps - 1.); // Energy stepsize
	}
	else{
		eMin = energymin;
		eMax = energymax;
		if(eMin > eMax){
			double tmp = eMin;
			eMin = eMax;
			eMax = tmp;
		}
		tMin = thetamin;
		tMax = thetamax;
		if(tMin > tMax){
			double tmp = tMin;
			tMin = tMax;
			tMax = tmp;
		}
		tStep = (tMax - tMin) / ((double)tSteps - 1.); // Theta stepsize
		eStep = (eMax - eMin) / ((double)eSteps - 1.); // Energy stepsize
	}

	if(fDetectorEff && fUseEfficiency){
		fDetectorEff_CM = new TGraph();
		double x,y;
		double tMin_tmp = 180;
		double tMax_tmp = 0;
		for(int n=0;n<fDetectorEff->GetN();n++){
			fDetectorEff->GetPoint(n,x,y);
			fDetectorEff_CM->SetPoint(n,fReaction->ConvertThetaLabToCm(x*TMath::DegToRad(),nPart)*TMath::RadToDeg(),y);
			if(y > 0 && fReaction->ConvertThetaLabToCm(x*TMath::DegToRad(),nPart)*TMath::RadToDeg() < tMin_tmp)
				tMin_tmp = fReaction->ConvertThetaLabToCm(x*TMath::DegToRad(),nPart)*TMath::RadToDeg();
			if(y > 0 && fReaction->ConvertThetaLabToCm(x*TMath::DegToRad(),nPart)*TMath::RadToDeg() > tMax_tmp)
				tMax_tmp = fReaction->ConvertThetaLabToCm(x*TMath::DegToRad(),nPart)*TMath::RadToDeg();
		}

		tMin = tMin_tmp;
		tMax = tMax_tmp;

	}

	int pointCounter = 0;
	for(int e = 0; e < eSteps; e++){
		tmpReac.SetLabEnergy(eMin + e*eStep);
		for(int t = 0; t < tSteps; t++){
			double eff = 1;
			if(fDetectorEff && fUseEfficiency){
				eff = fDetectorEff_CM->Eval(tMin + t*tStep);
			} 
			g->SetPoint(pointCounter,tMin + t*tStep, eMin + e*eStep, TMath::Sin((tMin + t*tStep) * TMath::DegToRad()) * tmpReac.Rutherford(tMin + t*tStep, nPart) * eff * TMath::DegToRad() * 2 * TMath::Pi());
			pointCounter++;
		}
	}
	
	return g;

}

TGraph2D* ExperimentRange::GetThetaEnergyGraph(int state){

	TGraph2D *g = new TGraph2D();

	int nPart = 2;
	if(targetDetection)
		nPart = 3;

	int npoint = 0;
	for(int nT = 0; nT < nTheta; nT++){
		for(int nE = 0; nE < nEnergy; nE++){
			double tmp_theta 	= fReaction->ConvertThetaCmToLab(fIntegral->GetCMThetaPoints().at(nE).at(nT) * TMath::DegToRad(), nPart) * TMath::RadToDeg();
			double tmp_cs 		= fIntegral->GetMeshPointProbabilities().at(nE).at(nT)[state] * TMath::Sin(tmp_theta * TMath::DegToRad()) * fReaction->Rutherford(tmp_theta,nPart) * TMath::DegToRad() * 2 * TMath::Pi();
			g->SetPoint(npoint,tmp_theta,fIntegral->GetEnergy(nE),tmp_cs);
			npoint++;
		}
	}

	char gname[64];
	sprintf(gname,"CrossSection_%f_%f_degrees_%f_%f_MeV",g->GetX()[0],g->GetX()[g->GetN()-1],g->GetY()[0],g->GetY()[g->GetN()-1]);
	g->SetName(gname);
	g->SetTitle(gname);

	return g;
	

}

void ExperimentRange::SetMeanValues(){


	double meantheta = (thetamax - thetamin)/2. + thetamin;
	double meanenergy = (energymax - energymin)/2. + energymin;

	int part = 2;
	if(targetDetection)
		part = 3;

	double meantheta_cm = fReaction->ConvertThetaLabToCm(meantheta * TMath::DegToRad(),part) * TMath::RadToDeg();
	meanThetaCM = meantheta_cm;
	meanEnergy  = meanenergy;

}


void ExperimentRange::PrintDetails() const{

	fNucleus->PrintNucleus();
	fReaction->PrintReaction();

	std::cout << "Mean theta [lab]: " << GetMeanThetaLab() << " degrees" << std::endl;
	std::cout << "Mean theta [CM]: " << GetMeanThetaCM() << " degrees" << std::endl;  
	std::cout << "Angular range [lab]: " << thetamin << "\t" << thetamax << " degrees" << std::endl;
	std::cout << "Energy range: " << energymin << "\t" << energymax << " MeV" << std::endl; 

}
