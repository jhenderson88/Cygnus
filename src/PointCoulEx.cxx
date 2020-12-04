#include "PointCoulEx.h"

// Empty constructor:
PointCoulEx::PointCoulEx(){

	bTargetDetection	= false;

	fTrack			= false;
	
	verbose 		= false;
	debug 			= false;

	projectileExcitation	= true;

	fAccuracy		= 1e-8;

	fUseFixedStep		= false;
	fUseSymmetry		= true;

}

// Construct with a nucleus:
PointCoulEx::PointCoulEx(Nucleus *nuc, Reaction *reac) 
{

	fNucleus = *nuc;
	fReaction = *reac;

	bTargetDetection	= false;

	fTrack			= false;

	verbose 		= false;
	debug 			= false;

	projectileExcitation	= true;

	fAccuracy		= 1e-8;

	fUseFixedStep		= false;
	fUseSymmetry		= true;

	// Because we have a reaction and  nucleus we can immediately define
	// connections between substates from state and reaction information
	PrepareConnections();

}

//	Copy constructor:
PointCoulEx::PointCoulEx(const PointCoulEx& p)
{

	fUseFixedStep		= p.fUseFixedStep;
	fUseSymmetry		= p.fUseSymmetry;

	bTargetDetection	= p.bTargetDetection;

	fTrack			= p.fTrack;
	fEpsilon		= p.fEpsilon;
	fOmegaTracking		= p.fOmegaTracking;
	fStateProbTracking	= p.fStateProbTracking;
 	fStateMultTracking	= p.fStateMultTracking;	
	
	fAccuracy		= p.fAccuracy;

	fStates.resize(p.fStates.size());
	for(size_t s = 0; s < p.fStates.size(); s++)
		fStates[s] = p.fStates[s];
	fSubstates.resize(p.fSubstates.size());
	for(size_t s = 0; s < p.fSubstates.size(); s++)
		fSubstates[s] = p.fSubstates[s];

	projectileExcitation	= p.projectileExcitation;

	fA1		= p.fA1;
	fA2		= p.fA2;
	fZ1		= p.fZ1;
	fZ2		= p.fZ2;

	fNucleus 	= p.fNucleus;
	fReaction 	= p.fReaction;

	verbose		= p.verbose;
	debug		= p.debug;
	
	fTheta		= p.fTheta;
	fSubStateProbabilities.ResizeTo(p.fSubStateProbabilities.GetNrows());
	fSubStateProbabilities	= p.fSubStateProbabilities;

	LMax		= p.LMax;

	XiMax		= p.XiMax;

	FinalRealAmplitude.ResizeTo(p.FinalRealAmplitude.GetNrows(),p.FinalRealAmplitude.GetNcols());
	FinalImagAmplitude.ResizeTo(p.FinalImagAmplitude.GetNrows(),p.FinalImagAmplitude.GetNcols());
	FinalRealAmplitude		= p.FinalRealAmplitude;
	FinalImagAmplitude		= p.FinalImagAmplitude;
	Probabilities.ResizeTo(p.Probabilities.GetNrows());
	Probabilities 			= p.Probabilities;

	fTensors	= p.fTensors;
	fTensorsB 	= p.fTensorsB;

}
//	Assignment operator:
PointCoulEx& PointCoulEx::operator = (const PointCoulEx &p)
{

	fUseFixedStep		= p.fUseFixedStep;
	fUseSymmetry		= p.fUseSymmetry;

	bTargetDetection	= p.bTargetDetection;

	fTrack			= p.fTrack;
	fEpsilon		= p.fEpsilon;
	fOmegaTracking		= p.fOmegaTracking;
	fStateProbTracking	= p.fStateProbTracking;
 	fStateMultTracking	= p.fStateMultTracking;	

	fAccuracy		= p.fAccuracy;

	fStates.resize(p.fStates.size());
	for(size_t s = 0; s < p.fStates.size(); s++)
		fStates[s] = p.fStates[s];
	fSubstates.resize(p.fSubstates.size());
	for(size_t s = 0; s < p.fSubstates.size(); s++)
		fSubstates[s] = p.fSubstates[s];;

	projectileExcitation	= p.projectileExcitation;

	fA1		= p.fA1;
	fA2		= p.fA2;
	fZ1		= p.fZ1;
	fZ2		= p.fZ2;

	fNucleus 	= p.fNucleus;
	fReaction 	= p.fReaction;

	verbose		= p.verbose;
	debug		= p.debug;
	
	fTheta		= p.fTheta;
	fSubStateProbabilities.ResizeTo(p.fSubStateProbabilities.GetNrows());
	fSubStateProbabilities	= p.fSubStateProbabilities;

	LMax		= p.LMax;

	XiMax		= p.XiMax;

	FinalRealAmplitude.ResizeTo(p.FinalRealAmplitude.GetNrows(),p.FinalRealAmplitude.GetNcols());
	FinalImagAmplitude.ResizeTo(p.FinalImagAmplitude.GetNrows(),p.FinalImagAmplitude.GetNcols());
	FinalRealAmplitude		= p.FinalRealAmplitude;
	FinalImagAmplitude		= p.FinalImagAmplitude;
	Probabilities.ResizeTo(p.Probabilities.GetNrows());
	Probabilities 			= p.Probabilities;

	fTensors	= p.fTensors;
	fTensorsB 	= p.fTensorsB;

	return *this;
	
}

//********************************* CALCULATION PREPATION SECTION  ***********************************//

//****************************************************************************************************//
//	Before we perform the calculations we will set up the variables defining the connections
//	between (sub)states. These are constant across the integration path and so only need to 
//	be calculated once.
//
//	In the event that the Nucleus (state energies) or reaction (beam energy) changes, these
//	variables need to be recalculated.
//****************************************************************************************************//

void PointCoulEx::PrepareConnections()
{

	//	For a non-zero ground-state spin, the ground state substates
	//	need to be accounted for:
	LMax = fNucleus.GetLevelJ().at(0) + 1;

	//	Indices for mass and proton number in the upcoming variable
	//	calculations depend on whether we're looking at projectile
	//	or target excitation:
	if(projectileExcitation){
		fZ1 	= (double)fReaction.GetTargetZ();
		fZ2	= (double)fReaction.GetBeamZ();
		fA1	= (double)fReaction.GetBeamA();
		fA2	= (double)fReaction.GetTargetA();
	}
	else{
		fZ1 	= (double)fReaction.GetBeamZ();
		fZ2	= (double)fReaction.GetTargetZ();
		fA1	= (double)fReaction.GetBeamA();
		fA2	= (double)fReaction.GetTargetA();
	}

	if(verbose){
		fNucleus.PrintNucleus();
		fReaction.PrintReaction();
	}

	//	The C coefficients depend on the multipolarity of the transition.
	//	For electric transitions:
	//	C(E,Lambda) = 1.116547 * pow(13.889122,Lambda) *
	//			(Lambda - 1)! / (2*Lambda + 1)!!
	//	For magnetic transitions:
	//	C(M,Lambda) = (v/c)C(E,Lambda)/95.0981942:
	double CPsi_Elec[6] = {5.169286,14.359366,56.982577,263.812653,1332.409500,7117.6915}; 
	double CPsi_Mag = fReaction.GetBeta() * CPsi_Elec[0] / 95.0981942; 

	//	Divding factor for determining the Psi variable
	double AAZZ = ((double)fReaction.GetBeamZ() * (double)fReaction.GetTargetZ()) * ((double)1 + fA1/fA2); 

	//	Before defining variables associated with the connections, define
	//	holders for the states and substates. Note that in the excitation
	//	process, substates are treated individually, with states only ho-
	//	lding information common to all substates.
	//
	//	See the State and Substate classes for information on holders
	fStates.clear();
	fSubstates.clear();
	for(int s = 0; s < fNucleus.GetNstates(); s++){
		State	tmpState(fNucleus.GetLevelEnergies().at(s),fNucleus.GetLevelJ().at(s));
		tmpState.SetEta(fReaction.EtaCalc(fNucleus.GetLevelEnergies().at(s)));
		tmpState.SetPsi(fReaction.GetLabEnergy() - (1 + fA1/fA2) * fNucleus.GetLevelEnergies().at(s));
		int phase = 1;
		if(s != 0)
			phase = TMath::Power(-1, 1 - tmpState.GetJ() - fStates.at(0).GetJ());
		tmpState.SetPhase(phase);
		fStates.push_back(tmpState);
		for(double ss = -fNucleus.GetLevelJ().at(s); ss <= fNucleus.GetLevelJ().at(s); ss++){
			if(ss == -0)
				ss = 0;
			Substate tmpSubstate(ss,fSubstates.size(),s);
			fSubstates.push_back(tmpSubstate);
		}
	}

	//	The maximum difference in wavenumber (XiMax) helps us select an
	//	appropriate integration range later
	XiMax	= 0;

	//	Because we are dealing with excitation, the important state is 
	//	the FINAL state, with the connections into that substate from
	//	other states being the important variables
	for(size_t ss1 = 0; ss1 < fSubstates.size(); ss1++){ // Substate 1 is the FINAL state

		//	Select the state corresponding to the substate ss1
		int	s1 	= fSubstates.at(ss1).GetStateIndex();
		for(size_t ss2 = 0; ss2 < fSubstates.size(); ss2++){ // Substate 2 is the INITIAL state
			//	Select the state corresponding to the substate ss1
			int 	s2 	= fSubstates.at(ss2).GetStateIndex();

			if(fSubstates.at(ss1).GetM() == -fSubstates.at(ss2).GetM() && s1 == s2){
				fSubstates.at(ss1).SetMirrorIndex(ss2);
				fSubstates.at(ss2).SetMirrorIndex(ss1);
			}

			//	Phase convention for matrix elements:
			int	mePhase	= TMath::Power(-1,(fStates.at(s2).GetJ()-fStates.at(s1).GetJ()));
			if(s2 > s1)
				mePhase = 1;

			double  psi1	= fStates.at(s1).GetPsi();
			double 	psi2	= fStates.at(s2).GetPsi();
			//	Phase factor for coupling parameter, depends on INITIAL (sub)state only
			double	phase	= TMath::Power(-1,fStates.at(s2).GetJ()-fSubstates.at(ss2).GetM());

			/*std::cout	<< std::setw(10) << std::left << phase
					<< std::setw(10) << std::left << fStates.at(s2).GetJ()
					<< std::setw(10) << std::left << fSubstates.at(ss2).GetM()
					<< std::endl;*/

			//	Difference in wavenumber between final and initial sttae
			double	tmpXi	= fStates.at(fSubstates.at(ss1).GetStateIndex()).GetEta() - fStates.at(fSubstates.at(ss2).GetStateIndex()).GetEta();

			//	Holder class for connection information, will be filled with all conne-
			//	ctions between substates ss1 and ss2 and then held in Substate ss1
			//
			//	See class Connection for details on connection holder class
			Connection tmpConnection;
			tmpConnection.SetConnectedState(ss2);
			//	For convenience, we allow all possible Lamda values (E1-6, M1)
			tmpConnection.SetMaxLambda(fNucleus.GetMaxLambda());

			//	Does this connection actually contain anything?			
			bool connectionFlag = false;
			for(int l=0;l<fNucleus.GetMaxLambda();l++){ // Loop over multipolarities

				//	If there is no matrix element between these states at this
				//	multipolarity, then we don't need to go any further
				if(fNucleus.GetMatrixElements().at(l)[s1][s2] == 0)
					continue;

				//	Because of the indexing convention for magnetic transitions
				//	(M1 = 7), we need to adjust our actual lambda accordingly
				int	lambda 	= l+1;	
				if(l > 5)
					lambda	= l-5;

				//	If the substate -> substate transition is forbidden on 
				//	angular momentum grounds, we can skip this
				if(TMath::Abs(fSubstates.at(ss1).GetM() - fSubstates.at(ss2).GetM()) > lambda)
					continue;
				if(TMath::Abs(fStates.at(fSubstates.at(ss1).GetStateIndex()).GetJ() - fStates.at(fSubstates.at(ss2).GetStateIndex()).GetJ()) > lambda)
					continue;

				//	Update our maximum wavenumber difference for use in the
				//	integration step later
				if(tmpXi > XiMax)
					XiMax = tmpXi;

				double 	Power	= (2*((double)lambda)-1)/4.;
			
				double	s	= TMath::Power(AAZZ,lambda);
				double 	C;
				if(l > 5)
					C	= CPsi_Mag / s;
				else
					C	= CPsi_Elec[l] / s;

				double	tmpPsi	= C * fZ1 * TMath::Sqrt(fA1) * TMath::Power( psi1 * psi2,Power);
			
				//	ThreeJ symbol uses 2* values (to allow for half-integer 
				//	spin) so we create a few new variables
				int	JInit	= 2 * fStates.at(fSubstates.at(ss2).GetStateIndex()).GetJ();
				int	L2	= 2 * lambda;
				int 	JFinal	= 2 * fStates.at(fSubstates.at(ss1).GetStateIndex()).GetJ();
				int	MInit	= 2 * fSubstates.at(ss2).GetM();
				int	Mu	= 2 * (fSubstates.at(ss2).GetM() - fSubstates.at(ss1).GetM());	
				int	MFinal	= 2 * fSubstates.at(ss1).GetM();

				double	tmpZeta	= TMath::Power((2*lambda + 1),0.5) * phase * fReaction.ThreeJ(JInit,L2,JFinal,-MInit,Mu,MFinal) * tmpPsi * mePhase;

				/*std::cout 	<< std::setw(8)  << std::left << fSubstates.at(ss1).GetM() 
						<< std::setw(8)  << std::left << fSubstates.at(ss2).GetM() 
						<< std::setw(15) << std::left << tmpZeta 
						<< std::setw(8)  << std::left << phase 
						<< std::setw(15) << std::left << tmpPsi 
						<< std::setw(15) << std::left << fReaction.ThreeJ(JInit,L2,JFinal,-MInit,Mu,MFinal)
						<< std::setw(12) << std::left << C
						<< std::setw(12) << std::left << Power
						<< std::endl;*/

				//	We've got a valid connection between substates for this
				//	multipolarity, so we can add it to the connection
				tmpConnection.AddLambda(l,tmpXi,tmpPsi,tmpZeta);	
				// 	We also now have something worth adding to our substate	
				connectionFlag = true;
			}
			//	If we've created a valid connection, add it to the final
			//	substate:
			if(connectionFlag)
				fSubstates.at(ss1).AddConnection(tmpConnection);
		}
	}


}

//*********************************** INTEGRATION SECTION BEGINS *************************************//

//****************************************************************************************************//
//	We have set up all of the necessary variables to perform an integration over the 
//	particle paths. Now we are ready to perform the integration. Two functions are
//	available, depending whether the probabilities vector is required immediately or
//	will be accessed later.	
//****************************************************************************************************//

TVectorD PointCoulEx::PointProbabilities(double theta)
{

	//	Preparing the connecitons takes very little time so do it every time, to
	//	be sure to include all changes
	PrepareConnections();
	TVectorD temp; 
	temp.ResizeTo(fNucleus.GetNstates());
	temp = Integration(theta); 
	Probabilities.ResizeTo(temp.GetNrows());
	Probabilities = temp;
	return temp;
}

void PointCoulEx::CalculatePointProbabilities(double theta)
{

	//	Preparing the connecitons takes very little time so do it every time, to
	//	be sure to include all changes	PrepareConnetions();
	TVectorD temp; 
	temp.ResizeTo(fNucleus.GetNstates());
	temp = Integration(theta); 
	Probabilities.ResizeTo(temp.GetNrows());
	Probabilities = temp;

}

//****************************************************************************************************//
//	This is where we're going to actually do the integration over the path.
//
//	Firstly, we need to redefine our coordinate system into something more sensible than x,y,z
//	
//	Then we can begin our integration, using some numerical differentiation methods. We begin 
//	with the more accurate Runge-Kutta method, this gives us the starting point. Once we have 
//	that, we move to the faster Adams-Moulton method. If we notice that our values lie outside 
//	a defined acceptable accuracy, we go back to the Runge-Kutta and redefine our incrementat-
//	ion to be smaller or larger, as necessary.
//****************************************************************************************************//


TVectorD PointCoulEx::Integration(double theta)
{

	fTheta = theta;

	// Here we define ourselves std::vectors of TMatrices and TVectors, note that we define 
	// vector.at(0) to be the real part and vector.at(1) to be the imaginary part
	std::vector<TMatrixD> Amplitude;
	std::vector<TMatrixD> AmplitudeP; 
	std::vector<TMatrixD> AmplitudeDerivative;
	// These two matrices (RealAmpF and ImagAmpF) give us the previous four steps for our
	//  numerical differentiation - differ from the aforementioned real/imaginary vector
	//  convention
	std::vector<TMatrixD> RealAmpF; 
	std::vector<TMatrixD> ImagAmpF;
	std::vector<TMatrixD> Q1_Matrix;

	// These are the vectors which we will eventually store our probabilities, here we 
	// define and resize them to appropriate values
	TVectorD SubStateProbabilities; 
	TVectorD StateProbabilities; 
	SubStateProbabilities.ResizeTo(fSubstates.size());
	StateProbabilities.ResizeTo(fNucleus.GetNstates());

	// Define an initial TotalProbability - soon to be overwritten
	double TotalProbability = 0;
	
	// Resize the vectors for their purpose (2 for most - real, imaginary, and 4 for 
	// RealAmpF and ImagAmpF, for the numerical differentiation)
	Amplitude.resize(2);
	AmplitudeP.resize(2);
	AmplitudeDerivative.resize(2);
	RealAmpF.resize(4);
	ImagAmpF.resize(4);
	Q1_Matrix.resize(2);

	// Resize the various matrices to the appropriate sizes
	Amplitude.at(0).ResizeTo(fSubstates.size(),LMax);
	Amplitude.at(1).ResizeTo(fSubstates.size(),LMax);
	AmplitudeP.at(0).ResizeTo(fSubstates.size(),LMax);
	AmplitudeP.at(1).ResizeTo(fSubstates.size(),LMax);
	AmplitudeDerivative.at(0).ResizeTo(fSubstates.size(),LMax);
	AmplitudeDerivative.at(1).ResizeTo(fSubstates.size(),LMax);
	RealAmpF.at(0).ResizeTo(fSubstates.size(),LMax);
	RealAmpF.at(1).ResizeTo(fSubstates.size(),LMax);
	RealAmpF.at(2).ResizeTo(fSubstates.size(),LMax);
	RealAmpF.at(3).ResizeTo(fSubstates.size(),LMax);
	ImagAmpF.at(0).ResizeTo(fSubstates.size(),LMax);
	ImagAmpF.at(1).ResizeTo(fSubstates.size(),LMax);
	ImagAmpF.at(2).ResizeTo(fSubstates.size(),LMax);
	ImagAmpF.at(3).ResizeTo(fSubstates.size(),LMax);
	Q1_Matrix.at(0).ResizeTo(fSubstates.size(),LMax);
	Q1_Matrix.at(1).ResizeTo(fSubstates.size(),LMax);

	double Up;
	double dOmega;	
	// Redefine theta into units of eliptical eccentricity, a more convenient system for
	// the upcoming calculation:
	double Epsilon = 1 / (TMath::Sin(theta * TMath::DegToRad() / 2.0));
	// We take our accuracy condition here to determine integration range and step size
	double Accuracy = fAccuracy; 
	double Accuracy_50 = Accuracy / 50.;
	// Defines the closest approach given our value of theta 
	double Distance = fReaction.ClosestApproach() * (1 + Epsilon) / 2.0;
	double ABW = 0.0;
	// Now we're onto the serious stuff, here we define the accuracy and the range and
	// steps of the integration we're going to perform. We also redefine everything into 
	// a more convenient coordinate system.
	//
	// For more details on the coordinate system especially, see GOSIA manual section 3.3 
	// - but note that omega = 0 corresponds to the closest approach
	//
	if(UseFixedStep()){
		// If we're using a fixed step (no adjustment - see INT flag in GOSIA) then we have
		// some preset values to use here. 
		if(fabs(MiscFunctions::GetMaxMatrix(fNucleus.GetMatrixElements().at(0))) > 1e-6)
			Up	= 13.12;
		else if (fabs(MiscFunctions::GetMaxMatrix(fNucleus.GetMatrixElements().at(6))) > 1e-6)
			Up	= 7.11;
		else if (fabs(MiscFunctions::GetMaxMatrix(fNucleus.GetMatrixElements().at(1))) > 1e-6)
			Up 	= 7.11;
		else if (fabs(MiscFunctions::GetMaxMatrix(fNucleus.GetMatrixElements().at(3))) > 1e-6)
			Up	= 5.14;
		else{	// No E1, M1, E2 or E3... this is weird (so default to E2 and print warning)
			std::cout	<< "No E1, M1, E2 or E2 matrix elements found, omega defaulting to E2"
					<< std::endl;
			Up	= 7.11;
		}
		dOmega = 0.03;
	}
	else{
		// The logic for determining integration range (Omega -> Up) is based on that used
		// in the CLX code, determined from the accuracy and maximum difference in wavenumber 
		// between initial and final states (XiMax):
		Up = TMath::Log(1 / (Epsilon * TMath::Sqrt(Accuracy))); 
		dOmega = (double) 40.0 * (pow(Accuracy,0.2)) / (10.0 + 48.0 * XiMax + 16.0 * XiMax * Epsilon); 
		int Steps = Up / (8.0*dOmega) + 1.;
		Steps *= 8;
		dOmega = Up / (Steps - 0.25);
		Up = dOmega * Steps;
	}
	double Omega = -Up; 

	if(verbose){ 
		std::cout	<< std::setw(24) << std::left << "Parabola excentricity: " 
				<< std::setw(8) << std::right << Epsilon 
				<< std::endl;
		std::cout	<< std::setw(24) << std::left << "Closest approach: "
				<< std::setw(8) << std::right << Distance
				<< std::endl;
		std::cout	<< "Integration between: " << Omega << " and " << Up << " d2w: " << 2*dOmega
				<< std::endl;
		std::cout 	<< "Largest Xi: " << XiMax 	
				<< std::endl;
	}

	fEpsilon		= Epsilon;
	fOmegaTracking.clear();
	fStateProbTracking.clear();
 	fStateMultTracking.clear();	

	// To begin with, the only non-zero amplitude should be the ground state, 
	// which should have a real amplitude of 1. Everything else should be zero.
	for(int i=0;i<LMax;i++) 
		Amplitude.at(0)[0][i] = 1;

	// The syntax below is based heavily on that used in the CLX Fortran code	

	// Now we get into the integration process itself, this is all based on 
	// numerical methods, no real physics going on:
	while(Omega <= Up) // Loop over Omega until it hits the end point
	{

		// First we're going to use the accurate Runge-Kutta method to get four 
		// intitial points: https://en.wikipedia.org/wiki/Runge-Kutta_methods 

		// Flag that tells us whether our omega steps are too small or large, 
		// we begin assuming our step-size is "good"	
		bool dOmegaFlag = true; 

		// Determine the amplitude derivatives at Omega
		AmplitudeDerivative = ComputeAmpDerivativeMatrices(Amplitude, Epsilon, Omega); 
		RealAmpF.at(0) = AmplitudeDerivative.at(0); 
		ImagAmpF.at(0) = AmplitudeDerivative.at(1);

		// Adams-Moulton requires four starting points, so we use Runge-Kutta
		// to determine another three:
		for(int n=1;n<4;n++) 
		{

			Q1_Matrix.at(0) = dOmega * AmplitudeDerivative.at(0);
			Q1_Matrix.at(1) = dOmega * AmplitudeDerivative.at(1); 
			Amplitude.at(0) = Amplitude.at(0) + Q1_Matrix.at(0); 
			Amplitude.at(1) = Amplitude.at(1) + Q1_Matrix.at(1);
			if(fNucleus.GetLevelJ()[0] == 0 && UseSymmetry()){
				for(size_t ss = 0; ss < fSubstates.size(); ss++){
					if(fSubstates.at(ss).GetM() < 0){
						int stateindex = fSubstates.at(ss).GetStateIndex();
						double par = pow(-1,fNucleus.GetLevelP()[0]-fNucleus.GetLevelP()[stateindex]-fNucleus.GetLevelJ()[stateindex]); 
						for(int i = 0; i < LMax; i++){
							Amplitude.at(0)[ss][i] = Amplitude.at(0)[fSubstates.at(ss).GetMirrorIndex()][i] * par; 
							Amplitude.at(1)[ss][i] = Amplitude.at(1)[fSubstates.at(ss).GetMirrorIndex()][i] * par; 
						}
					}
				}
			}		

			// Increment Omega
			Omega = Omega + dOmega; 
			
			if(fTrack){
				// Check excitation probabilities
				SubStateProbabilities.Clear();
				SubStateProbabilities.ResizeTo(fSubstates.size());
				for(unsigned int i=0;i<fSubstates.size();i++)
					for(int j=0;j<LMax;j++)
						SubStateProbabilities[i] = SubStateProbabilities[i] + (pow(Amplitude.at(0)[i][j],2) + pow(Amplitude.at(1)[i][j],2));

				StateProbabilities.Clear();
				StateProbabilities.ResizeTo(fNucleus.GetNstates());
				for(size_t i=0;i<fSubstates.size();i++)
					StateProbabilities[fSubstates.at(i).GetStateIndex()] = StateProbabilities[fSubstates.at(i).GetStateIndex()] + SubStateProbabilities[i];

				fOmegaTracking.push_back(Omega);
				std::vector<double>	tmpProb;
				for(int i=0;i<StateProbabilities.GetNrows();i++){
					tmpProb.push_back(StateProbabilities[i]);
				}
				fStateProbTracking.push_back(tmpProb);	
			}
	
			AmplitudeDerivative = ComputeAmpDerivativeMatrices(Amplitude, Epsilon, Omega); 

			Amplitude.at(0) = Amplitude.at(0) + 0.5857864 * (dOmega*AmplitudeDerivative.at(0) - Q1_Matrix.at(0)); 
			Amplitude.at(1) = Amplitude.at(1) + 0.5857864 * (dOmega*AmplitudeDerivative.at(1) - Q1_Matrix.at(1)); 
			Q1_Matrix.at(0) = 0.5857864 * dOmega * AmplitudeDerivative.at(0) + 0.1213204 * Q1_Matrix.at(0); 
			Q1_Matrix.at(1) = 0.5857864 * dOmega * AmplitudeDerivative.at(1) + 0.1213204 * Q1_Matrix.at(1);
			if(fNucleus.GetLevelJ()[0] == 0 && UseSymmetry()){
				for(size_t ss = 0; ss < fSubstates.size(); ss++){
					if(fSubstates.at(ss).GetM() < 0){
						int stateindex = fSubstates.at(ss).GetStateIndex();
						double par = pow(-1,fNucleus.GetLevelP()[0]-fNucleus.GetLevelP()[stateindex]-fNucleus.GetLevelJ()[stateindex]); 
						for(int i = 0; i < LMax; i++){
							Amplitude.at(0)[ss][i] = Amplitude.at(0)[fSubstates.at(ss).GetMirrorIndex()][i] * par; 
							Amplitude.at(1)[ss][i] = Amplitude.at(1)[fSubstates.at(ss).GetMirrorIndex()][i] * par; 
						}
					}
				}
			}				

			// Determine amplitude derivative at omega based on updated initial amplitude
			AmplitudeDerivative = ComputeAmpDerivativeMatrices(Amplitude, Epsilon, Omega); 
			
			Amplitude.at(0) = Amplitude.at(0) + 3.414214 * (dOmega * AmplitudeDerivative.at(0) - Q1_Matrix.at(0)); 
			Amplitude.at(1) = Amplitude.at(1) + 3.414214 * (dOmega * AmplitudeDerivative.at(1) - Q1_Matrix.at(1)); 
			Q1_Matrix.at(0) = 3.414214 * dOmega * AmplitudeDerivative.at(0) - 4.1213204 * Q1_Matrix.at(0); 
			Q1_Matrix.at(1) = 3.414214 * dOmega * AmplitudeDerivative.at(1) - 4.1213204 * Q1_Matrix.at(1); 
			if(fNucleus.GetLevelJ()[0] == 0 && UseSymmetry()){
				for(size_t ss = 0; ss < fSubstates.size(); ss++){
					if(fSubstates.at(ss).GetM() < 0){
						int stateindex = fSubstates.at(ss).GetStateIndex();
						double par = pow(-1,fNucleus.GetLevelP()[0]-fNucleus.GetLevelP()[stateindex]-fNucleus.GetLevelJ()[stateindex]); 
						for(int i = 0; i < LMax; i++){
							Amplitude.at(0)[ss][i] = Amplitude.at(0)[fSubstates.at(ss).GetMirrorIndex()][i] * par; 
							Amplitude.at(1)[ss][i] = Amplitude.at(1)[fSubstates.at(ss).GetMirrorIndex()][i] * par; 
						}
					}
				}
			}				

			Omega = Omega + dOmega;

			// Determine amplitude derivative at new omega
			AmplitudeDerivative = ComputeAmpDerivativeMatrices(Amplitude, Epsilon, Omega);

			Amplitude.at(0) = Amplitude.at(0) + (1./3.) * dOmega * AmplitudeDerivative.at(0) - (2./3.) * Q1_Matrix.at(0); 
			Amplitude.at(1) = Amplitude.at(1) + (1./3.) * dOmega * AmplitudeDerivative.at(1) - (2./3.) * Q1_Matrix.at(1);
			if(fNucleus.GetLevelJ()[0] == 0 && UseSymmetry()){
				for(size_t ss = 0; ss < fSubstates.size(); ss++){
					if(fSubstates.at(ss).GetM() < 0){
						int stateindex = fSubstates.at(ss).GetStateIndex();
						double par = pow(-1,fNucleus.GetLevelP()[0]-fNucleus.GetLevelP()[stateindex]-fNucleus.GetLevelJ()[stateindex]); 
						for(int i = 0; i < LMax; i++){
							Amplitude.at(0)[ss][i] = Amplitude.at(0)[fSubstates.at(ss).GetMirrorIndex()][i] * par; 
							Amplitude.at(1)[ss][i] = Amplitude.at(1)[fSubstates.at(ss).GetMirrorIndex()][i] * par; 
						}
					}
				}
			}				
			
			if(fTrack){
				// Check excitation probabilities
				SubStateProbabilities.Clear();
				SubStateProbabilities.ResizeTo(fSubstates.size());
				for(unsigned int i=0;i<fSubstates.size();i++)
					for(int j=0;j<LMax;j++)
						SubStateProbabilities[i] = SubStateProbabilities[i] + (pow(Amplitude.at(0)[i][j],2) + pow(Amplitude.at(1)[i][j],2));

				StateProbabilities.Clear();
				StateProbabilities.ResizeTo(fNucleus.GetNstates());
				for(size_t i=0;i<fSubstates.size();i++)
					StateProbabilities[fSubstates.at(i).GetStateIndex()] = StateProbabilities[fSubstates.at(i).GetStateIndex()] + SubStateProbabilities[i];

				fOmegaTracking.push_back(Omega);
				std::vector<double>	tmpProb;
				for(int i=0;i<StateProbabilities.GetNrows();i++){
					tmpProb.push_back(StateProbabilities[i]);
				}
				fStateProbTracking.push_back(tmpProb);	
			}
 
			// Determine amplitude derivative at omega based on updated initial amplitude
			AmplitudeDerivative = ComputeAmpDerivativeMatrices(Amplitude, Epsilon, Omega);

			RealAmpF.at(n) = AmplitudeDerivative.at(0);
			ImagAmpF.at(n) = AmplitudeDerivative.at(1);

		} 

		// The Runge-Kutta has given us our initial points, we can use the fast 
		// Adams-Moulton method to get the rest done: 
		// https://en.wikipedia.org/wiki/Linear_multistep_method#Adams-Moulton_methods 

		// If our step size has to change then we'll need to go back to the 
		// Runge-Kutta to start again 
		while( Omega < Up && dOmegaFlag) 
		{

			AmplitudeP.at(0) = Amplitude.at(0) + (dOmega/12.0) * (55.0 * RealAmpF.at(3) - 59.0 * RealAmpF.at(2) + 37.0 * RealAmpF.at(1) - 9.0 * RealAmpF.at(0));
			AmplitudeP.at(1) = Amplitude.at(1) + (dOmega/12.0) * (55.0 * ImagAmpF.at(3) - 59.0 * ImagAmpF.at(2) + 37.0 * ImagAmpF.at(1) - 9.0 * ImagAmpF.at(0)); 
			if(fNucleus.GetLevelJ()[0] == 0 && UseSymmetry()){
				for(size_t ss = 0; ss < fSubstates.size(); ss++){
					if(fSubstates.at(ss).GetM() < 0){
						int stateindex = fSubstates.at(ss).GetStateIndex();
						double par = pow(-1,fNucleus.GetLevelP()[0]-fNucleus.GetLevelP()[stateindex]-fNucleus.GetLevelJ()[stateindex]); 
						for(int i = 0; i < LMax; i++){
							AmplitudeP.at(0)[ss][i] = AmplitudeP.at(0)[fSubstates.at(ss).GetMirrorIndex()][i] * par; 
							AmplitudeP.at(1)[ss][i] = AmplitudeP.at(1)[fSubstates.at(ss).GetMirrorIndex()][i] * par; 
						}
					}
				}
			}
	
			Omega = Omega + dOmega + dOmega;

			// Determine amplitude derivative at new omega
			AmplitudeDerivative = ComputeAmpDerivativeMatrices(AmplitudeP, Epsilon, Omega);

			Amplitude.at(0) = Amplitude.at(0) + (dOmega/12.0) * (9.0 * AmplitudeDerivative.at(0) + 19.0 * RealAmpF.at(3) - 5.0 * RealAmpF.at(2) + RealAmpF.at(1));
			Amplitude.at(1) = Amplitude.at(1) + (dOmega/12.0) * (9.0 * AmplitudeDerivative.at(1) + 19.0 * ImagAmpF.at(3) - 5.0 * ImagAmpF.at(2) + ImagAmpF.at(1));
			if(fNucleus.GetLevelJ()[0] == 0 && UseSymmetry()){
				for(size_t ss = 0; ss < fSubstates.size(); ss++){
					if(fSubstates.at(ss).GetM() < 0){
						int stateindex = fSubstates.at(ss).GetStateIndex();
						double par = pow(-1,fNucleus.GetLevelP()[0]-fNucleus.GetLevelP()[stateindex]-fNucleus.GetLevelJ()[stateindex]); 
						for(int i = 0; i < LMax; i++){
							Amplitude.at(0)[ss][i] = Amplitude.at(0)[fSubstates.at(ss).GetMirrorIndex()][i] * par; 
							Amplitude.at(1)[ss][i] = Amplitude.at(1)[fSubstates.at(ss).GetMirrorIndex()][i] * par; 
						}
					}
				}
			}
	
			// Determine amplitude derivative at omega based on updated initial amplitude
			AmplitudeDerivative = ComputeAmpDerivativeMatrices(Amplitude, Epsilon, Omega);

			// Shift all of the differential points along by one - Adams-Moulton only 
			// requires three points:
			RealAmpF.at(0) = RealAmpF.at(1);
			ImagAmpF.at(0) = ImagAmpF.at(1);
			RealAmpF.at(1) = RealAmpF.at(2);
			ImagAmpF.at(1) = ImagAmpF.at(2);	
			RealAmpF.at(2) = RealAmpF.at(3);
			ImagAmpF.at(2) = ImagAmpF.at(3);	
			RealAmpF.at(3) = AmplitudeDerivative.at(0);
			ImagAmpF.at(3) = AmplitudeDerivative.at(1);	

			// Check amplitudes to see if the stepsize needs changing
			if(Omega + dOmega <= Up && !UseFixedStep())
			{
				double FF=0;
				for(int i=0;i<LMax;i++)
				{
					for(unsigned int j=0;j<fSubstates.size();j++)
					{
						double FZR = AmplitudeP.at(0)[j][i] - Amplitude.at(0)[j][i];
						double FZI = AmplitudeP.at(1)[j][i] - Amplitude.at(1)[j][i];
						double FZ = TMath::Sqrt(pow(FZR,2) + pow(FZI,2)) * 19. / 270.;
						if(FZ > FF) FF = FZ;

					}
				}

				if(FF <= Accuracy_50)
				{
					dOmegaFlag = false;
					dOmega = 1.5*dOmega;
					if(verbose){
						std::cout	<< "At Omega: " << std::setw(12) << std::left << Omega
								<< "Stepwidth doubled to: " << std::setw(12) << std::left << 2*dOmega
								<< std::endl;
					}
				}
				if(FF > Accuracy)
				{
					dOmegaFlag = false;
					dOmega = (1./1.5)*dOmega;
					if(verbose){
						std::cout	<< "At Omega: " << std::setw(12) << std::left << Omega
								<< "Stepwidth halved to:  " << std::setw(12) << std::left << 2*dOmega
								<< std::endl;
					}
				}
								
			}

			// Check excitation probabilities
			SubStateProbabilities.Clear();
			SubStateProbabilities.ResizeTo(fSubstates.size());
			for(unsigned int i=0;i<fSubstates.size();i++)
				for(int j=0;j<LMax;j++)
					SubStateProbabilities[i] = SubStateProbabilities[i] + (pow(Amplitude.at(0)[i][j],2) + pow(Amplitude.at(1)[i][j],2));

			StateProbabilities.Clear();
			StateProbabilities.ResizeTo(fNucleus.GetNstates());
			for(size_t i=0;i<fSubstates.size();i++)
				StateProbabilities[fSubstates.at(i).GetStateIndex()] = StateProbabilities[fSubstates.at(i).GetStateIndex()] + SubStateProbabilities[i];

			if(fTrack){
				fOmegaTracking.push_back(Omega);
				std::vector<double>	tmpProb;
				for(int i=0;i<StateProbabilities.GetNrows();i++){
					tmpProb.push_back(StateProbabilities[i]);
				}
				fStateProbTracking.push_back(tmpProb);	
			}

			TotalProbability = 0;
			for(int i=0;i<fNucleus.GetNstates();i++)
				TotalProbability = TotalProbability + StateProbabilities[i];

			if(fabs(TotalProbability - 1) > fabs(ABW))
				ABW = TotalProbability - 1.0;

		}

	}

	// Now the integration over all omega is done we can put the results in our 
	// storage objects

	fSubStateProbabilities.ResizeTo(fSubstates.size());
	FinalRealAmplitude.ResizeTo(fSubstates.size(),LMax);
	FinalImagAmplitude.ResizeTo(fSubstates.size(),LMax);

	fSubStateProbabilities = SubStateProbabilities;
	FinalRealAmplitude = Amplitude.at(0);
	FinalImagAmplitude = Amplitude.at(1);

	if(verbose){
		StateProbabilities.Print();
		fSubStateProbabilities.Print();
		FinalRealAmplitude.Print();
		FinalImagAmplitude.Print();
	}

	return StateProbabilities;


}

//****************************************************************************************************//
//	This function evaluates the Amplitude Derivatives at a given point in the integration 
//	process. The amplitude derivative function is defined in the GOSIA manual, Eq. 3.20
//****************************************************************************************************//

std::vector<TMatrixD> PointCoulEx::ComputeAmpDerivativeMatrices(std::vector<TMatrixD> AmplitudeIn, double Epsilon, double Omega)
{

	double	omega_max[7] = {13.12,7.11,5.14,4.17,3.59,3.26,7.11};

	// Calculate the exponent from Eq. 3.20 in the GOSIA manual
	double RAlfa = Epsilon * TMath::SinH(Omega) + Omega; 

	// Clear all vectors and resize them appropriately:
	std::vector<TMatrixD> AmpDot; 
	AmpDot.resize(2);
	AmpDot.at(0).ResizeTo(fSubstates.size(),LMax);
	AmpDot.at(1).ResizeTo(fSubstates.size(),LMax);
	for(size_t f = 0; f < fSubstates.size(); f++){
		for(int l = 0; l < LMax; l++){
			AmpDot.at(0)[f][l] = 0;
			AmpDot.at(1)[f][l] = 0;
		}
	}


	for(int l=0; l < fNucleus.GetMaxLambda(); l++){

		// No need to continue with this lambda if there aren't any matrix elements
		if(fabs(MiscFunctions::GetMaxMatrix(fNucleus.GetMatrixElements().at(l))) < 1e-6)
			continue;

		if(UseFixedStep() && fabs(Omega) > omega_max[l])
			continue;
			
		// Grab the collision functions for this lambda, omega, epsilon combination
		std::vector< std::complex<double> > Q_Matrix = CollisionFunction((l+1),Epsilon,Omega);

		for(size_t fs = 0; fs < fSubstates.size(); fs++){ 	// Loop over all substates - this is the "final state" of this step
			if(UseSymmetry() && fNucleus.GetLevelJ()[0] == 0 && fSubstates.at(fs).GetM() < 0)
				continue;
			for(size_t c = 0; c < fSubstates.at(fs).GetNconnections(); c++){ 	// Loop over all connections to this final state
				Connection tmpConnection = fSubstates.at(fs).GetConnection(c);	// Grab the connection	
				int is = tmpConnection.GetConnectedState();			

				// If this connection for this lambda was not set earlier (i.e. didn't
				// satisfy angular moment conservation or similar) then we don't need
				// to go any further. This is important where a state is connected by 
				// multiple connections of different multipolarities (e.g. E2/M1)
				if(!tmpConnection.IsSet(l))
					continue;

				int finalState 		= fSubstates.at(fs).GetStateIndex();
				int initialState	= fSubstates.at(is).GetStateIndex();

				if(fNucleus.GetMatrixElements().at(l)[finalState][initialState] != 0){	// If there is a non-zero matrix element
			
					// Determine the real and imaginary parts of the exponent:
					double	RealEx	= TMath::Cos(tmpConnection.GetXi(l) * RAlfa);
					double	ImagEx	= TMath::Sin(tmpConnection.GetXi(l) * RAlfa);
				
					// Difference in magnetic substates for selection of the 
					// correct collision function:
					double	MuAbs	= (int)(fabs(fSubstates.at(is).GetM() - fSubstates.at(fs).GetM()));

					// Combine all of the variables and the above information
					// to give real and imaginary amplitudes for this substate
					// from this connection:
					double	RealRC	= (std::real(Q_Matrix.at(MuAbs))*RealEx	- std::imag(Q_Matrix.at(MuAbs))*ImagEx) * tmpConnection.GetZeta(l) * fNucleus.GetMatrixElements().at(l)[finalState][initialState];
					double	ImagRC	= (std::real(Q_Matrix.at(MuAbs))*ImagEx	+ std::imag(Q_Matrix.at(MuAbs))*RealEx) * tmpConnection.GetZeta(l) * fNucleus.GetMatrixElements().at(l)[finalState][initialState];
			
					// If we're dealing with magnetic transitions we also have
					// to include an additional element due to the asymmetry in
					// M-state selection:
					if(l > 5){
							if((fSubstates.at(fs).GetM() - fSubstates.at(is).GetM()) < 0){
								RealRC *= -1;
								ImagRC *= -1;
							}
					}

					// Add the newly determined amplitude for the final state, multiplied
					// by the amplitude in the initial state:
					for(int n=0;n<LMax;n++){
						double 	RRC = RealRC * AmplitudeIn.at(0)[is][n] - ImagRC * AmplitudeIn.at(1)[is][n];
						double	IRC = RealRC * AmplitudeIn.at(1)[is][n] + ImagRC * AmplitudeIn.at(0)[is][n];
						AmpDot.at(0)[fs][n]	+= IRC;
						AmpDot.at(1)[fs][n]	-= RRC;
					}
				}
			}
		} 
	}

	return AmpDot;

}
//****************************************************************************************************//
// 	Calculates collision functions, dependent upon lambda, epsilon and omega for the possible 
// 	values of mu. This is returned as a vector of complex numbers. Even mu values correspond 
// 	to real collision functions, odd mu values correspond to imaginary collision functions.
//****************************************************************************************************//

std::vector< std::complex<double> > PointCoulEx::CollisionFunction(int Lambda, double Epsilon, double Omega)
{

	// The vector of complex numbers which we will be filling
	// and returning
	std::vector< std::complex<double> > output; 

	// Define the Epsilon and Omega dependent variables, as described 
	// in Table 2.1a of the GOSIA manual
	double a = TMath::CosH(Omega) + Epsilon; 
	double b = Epsilon * TMath::CosH(Omega) + 1.;
	double c = TMath::Sqrt(pow(Epsilon,2) - 1.) * TMath::SinH(Omega);

	// Temporary real and imaginary values
	double Real; 
	double Imag;

	// Polarisation factor to account for dipole polarisation contributions
	double PolFactor = (1. - ElectricDipoleNormalization()/b);

	// E1 Collision Functions:
	if(Lambda == 1)
	{
		
		output.resize(2);
		Real = 0.5 * (a / pow(b,2) );
		Imag = 0.0;		
		output.at(0) = std::complex<double>(Real,Imag);

		Real = 0.0;
		Imag = -(1.0 / (2.0 * TMath::Sqrt(2) ) ) * c / pow(b,2);
		output.at(1.0) = std::complex<double>(Real,Imag);
		
		return output;

	}

	// E2 Collision Functions:
	if(Lambda == 2)
	{
		
		output.resize(3);
		Real = (3./4.) * (2.*pow(a,2) - pow(c,2))/pow(b,4) * PolFactor;
		Imag = 0.;
		output.at(0) = std::complex<double>(Real,Imag);

		Real = 0.;
		Imag = -((3. * TMath::Sqrt(3))/(2. * TMath::Sqrt(2) )) * (a * c)/pow(b,4) * PolFactor;
		output.at(1) = std::complex<double>(Real,Imag);

		Real = -((3. * TMath::Sqrt(3)) / (4. * TMath::Sqrt(2) )) * pow(c,2)/pow(b,4) * PolFactor;
		Imag = 0.;		
		output.at(2) = std::complex<double>(Real,Imag);
		
		return output;

	}

	// E3 Collision Functions:
	if(Lambda == 3)
	{
		
		output.resize(4);
		Real = (15./8.) * a*(2*pow(a,2) - 3*pow(c,2))/pow(b,6);
		Imag = 0.;		
		output.at(0) = std::complex<double>(Real,Imag);

		Real = 0.;
		Imag = -((15. * TMath::Sqrt(3))/16.) * c * (4 * pow(a,2) - pow(c,2))/pow(b,6);
		output.at(1) = std::complex<double>(Real,Imag);

		Real = -((15. * TMath::Sqrt(15)) / (8 * TMath::Sqrt(2) )) * a*pow(c,2)/pow(b,6);
		Imag = 0.;		
		output.at(2) = std::complex<double>(Real,Imag);

		Real = 0.;
		Imag = ((15. * TMath::Sqrt(5))/16.) * pow(c,3) / pow(b,6);
		output.at(3) = std::complex<double>(Real,Imag);
	
		return output;

	}	

	// E4 Collision Functions:
	if(Lambda == 4)
	{
		
		output.resize(5);
		Real = (35/32) * (8*pow(a,4) - 24*pow(a,2)*pow(c,2) + 3*pow(c,4)) / pow(b,8);
		Imag = 0.;		
		output.at(0) = std::complex<double>(Real,Imag);

		Real = 0.;
		Imag = -((35 * TMath::Sqrt(5))/16) * a*c * (4 * pow(a,2) - 3*pow(c,2))/pow(b,8);
		output.at(1) = std::complex<double>(Real,Imag);

		Real = -((35 * TMath::Sqrt(5)) / (16 * TMath::Sqrt(2) )) * pow(c,2) * (6*pow(a,2) - pow(c,2)) / pow(b,8);
		Imag = 0.;		
		output.at(2) = std::complex<double>(Real,Imag);

		Real = 0.;
		Imag = ((35 * TMath::Sqrt(35))/16) * a * pow(c,3) / pow(b,8);
		output.at(3) = std::complex<double>(Real,Imag);

		Real = (35 * TMath::Sqrt(35))/(32 * TMath::Sqrt(2)) * pow(c,4)/pow(b,8);
		Imag = 0.;		
		output.at(4) = std::complex<double>(Real,Imag);
		
		return output;

	}	

	if(Lambda == 7)
	{

		output.resize(2);
		Real = 0.;
		Imag = 0.;
		output.at(0) = std::complex<double>(Real,Imag);

		Real = -0.3535533905 * TMath::Sqrt(TMath::Power(Epsilon,2) - 1) / TMath::Power(b,2);
		Imag = 0.;
		output.at(1) = std::complex<double>(Real,Imag);		

		return output;
	}

	// This will be extended. In addition, magnetic transitions will eventually be included.
	std::cout << "Error: Cannot process electric multipolarities greater than E4 or magnetic multipolarities other than M1" << std::endl;

	return output;

}

//****************************************************************************************************//
//	Calculate Electric Dipole Normalization value this accounts for E1 strength. See GOSIA 
//	manual Eqs. 3.27 and 3.28 for description.
//****************************************************************************************************//

double PointCoulEx::ElectricDipoleNormalization()
{

	double factor;
	if(projectileExcitation){
		factor = Reaction::dipole * (fReaction.GetLabEnergy() * fReaction.GetBeamA()) / (pow(fReaction.GetBeamZ(),2) * (1. + (double)fReaction.GetBeamA()/(double)fReaction.GetTargetA()));
	}
	else{
		factor = Reaction::dipole * (fReaction.GetLabEnergy() * fReaction.GetTargetA()) / (pow(fReaction.GetTargetZ(),2) * (1 + (double)fReaction.GetBeamA()/(double)fReaction.GetTargetA()));
	}
	return factor;

}

//************************************ INTEGRATION SECTION ENDS **************************************//

//*********************************** TENSOR CALCULATION BEGINS **************************************//
//****************************************************************************************************//
//	Here we will calculate the statistical tensors which define the gamma-ray angular
//	distribution based on the M-state amplitudes
//	
//	The logic behind these calculations is found in the GOSIA manual
//
//	First stage, for each state:
//	rho_kX(I) = sqrt(2*I + 1)/(2*I_0 + 1) * SUM(M_0,M,M'){
//		(-1)^(I - M') * ThreeJ(I,k,I,-M',X,M) *
//			amplitude*(IM')[M0] * amplitude(IM)[M0]	
//	}
//	k 	= 0, 2, 4 		(as with orders of Legendre polynomials)
//	X 	= 0,1,....,k-1,k 	(kappa in the below code)
//	I 	= state spin
//	I_0	= ground-state spin
//	M_0	= ground-state m (= 0 for J=0 g.s.)
//	M	= substate 1
//	M'	= substate 2
//
//****************************************************************************************************//

void PointCoulEx::CalculateTensors(){
	
	CalculateTensorsExcitationFrame();
	CalculateTensorsLabFrame();

}

void PointCoulEx::CalculateTensorsLabFrame(){

	StatisticalTensor statTensors;

	double beta	= (TMath::Pi() + TMath::DegToRad()*fTheta)/2.;

	for(size_t s1  = 1; s1  < fStates.size(); s1++){
	

		int kmax = 6;
		if(2*fStates.at(s1).GetJ() < 6)
			kmax = fStates.at(s1).GetJ()*2;

		StateTensor tmpTensor;
		tmpTensor.SetState(s1);


		//	Loop over k = 0,2,4,6 (Legendre polynomial coefficients, basically)
		for(int k=0; k<=kmax; k+=2){ 

			if(k==0){
				int index = fTensorsB.GetStateTensor(s1-1).IndexFromKkappa(0,0);
				tmpTensor.AddElement(fTensorsB.GetStateTensor(s1-1).GetTensor(index),0,0);	
			}
			else{ 
				//	Loop over kappa
				for(int kappa = 0; kappa<= k; ++kappa){		
					double tmpProb = 0;
					for(int k_kappa = 0; k_kappa <= k; ++k_kappa){
						int phase 	= TMath::Power(-1,(kappa + (int)(k_kappa/2.)));	
						int index	= fTensorsB.GetStateTensor(s1-1).IndexFromKkappa(k,k_kappa);
						double tensorB  = fTensorsB.GetStateTensor(s1-1).GetTensor(index);
						double DJMM	= MiscFunctions::RotationFunction(beta,k,k_kappa,kappa);
						tmpProb += tensorB * phase * DJMM;
						if(k_kappa != 0){
							phase 		= TMath::Power(-1,(kappa - (int)((k_kappa+1)/2.)));	
							DJMM		= MiscFunctions::RotationFunction(beta,k,-k_kappa,kappa);
							tmpProb 	+= tensorB * phase * DJMM;
						}
					}
					tmpTensor.AddElement(tmpProb,k,kappa);
				}
			}
		}
		statTensors.AddStateTensor(tmpTensor);
	}

	fTensors	= statTensors;

	if(verbose){
		std::cout	<< std::setw(16) << std::left << "Lab Frame:"
				<< std::endl;
		std::cout	<< std::setw(10) << std::left << "INDEX:"
				<< std::setw(10) << std::left << "KA:"
				<< std::setw(10) << std::left << "KAPPA:"
				<< std::setw(14) << std::left << "RHOC:"
				<< std::endl;
		for(int i=0;i<fTensors.GetNstates();i++){	
			StateTensor tmpTensor = fTensors.GetStateTensor(i);
			for(int j=0;j<tmpTensor.GetNelements();j++){
				std::cout 	<< std::setw(10) << std::left << tmpTensor.GetState()
						<< std::setw(10) << std::left << tmpTensor.GetK(j)
						<< std::setw(10) << std::left << tmpTensor.GetKappa(j)
						<< std::setw(14) << std::left << tmpTensor.GetTensor(j)
						<< std::endl;
			}
		}
	}

}

void PointCoulEx::CalculateTensorsExcitationFrame(){

	StatisticalTensor statTensors;

	for(size_t s1  = 1; s1  < fStates.size(); s1++){
		

		bool flag 		= true;
		int init_substate 	= 0;
		int final_substate 	= 0;
		for(size_t ss1 = 0; ss1 < fSubstates.size(); ss1++){
			if(fSubstates.at(ss1).GetStateIndex()==(int)s1){
				if(flag){
					init_substate = ss1;
					flag = false;
				}
				final_substate = ss1;
			}
		}  

		double leadingFactor = TMath::Sqrt(2*fStates.at(s1).GetJ() + 1) / (2 * fStates.at(0).GetJ()+1);

		StateTensor tmpTensor;
		tmpTensor.SetState(s1);

		int kmax = 6;
		if(2*fStates.at(s1).GetJ() < 6)
			kmax = fStates.at(s1).GetJ()*2;

		//	Loop over k = 0,2,4,6 (Legendre polynomial coefficients, basically)
		for(int k=0; k<=kmax; k+=2){ 
			//	Loop over kappa
			for(int kappa = 0; kappa<= k; ++kappa){
				double tmpRho  = 0;
				double tmpProb = 0;

				for(size_t ss1 = init_substate; ss1 <= (size_t)final_substate; ss1++){
					if(((int)ss1 - kappa) < init_substate)
						continue;
					size_t ss2 	= ss1 - kappa;

					int phase 	= TMath::Power(-1,(fStates.at(s1).GetJ() - fSubstates.at(ss2).GetM()));

					double J2 	= 2 * fStates.at(s1).GetJ();
					double M1 	= 2 * fSubstates.at(ss1).GetM();
					double M2 	= 2 * fSubstates.at(ss2).GetM();

					double threeJ 	= fReaction.ThreeJ(J2,2*k,J2,-M1,2*kappa,M2);
				
					if( (kappa % 2) == 0)
						tmpProb = (FinalRealAmplitude[ss1][0] * FinalRealAmplitude[ss2][0] + FinalImagAmplitude[ss1][0] * FinalImagAmplitude[ss2][0]);
					else
						tmpProb = (FinalRealAmplitude[ss2][0] * FinalImagAmplitude[ss1][0] - FinalImagAmplitude[ss2][0] * FinalRealAmplitude[ss1][0]);

					tmpRho += tmpProb * phase * threeJ;

				}
				tmpRho *= leadingFactor;
				tmpTensor.AddElement(tmpRho,k,kappa);
			}
		}
		statTensors.AddStateTensor(tmpTensor);
	}

	fTensorsB	= statTensors;

	if(verbose){
		std::cout	<< std::setw(16) << std::left << "Excitation Frame:"
				<< std::endl;
		std::cout	<< std::setw(10) << std::left << "INDEX:"
				<< std::setw(10) << std::left << "KA:"
				<< std::setw(10) << std::left << "KAPPA:"
				<< std::setw(14) << std::left << "RHOB:"
				<< std::endl;
		for(int i=0;i<fTensorsB.GetNstates();i++){	
			StateTensor tmpTensor = fTensorsB.GetStateTensor(i);
			for(int j=0;j<tmpTensor.GetNelements();j++){
				std::cout 	<< std::setw(10) << std::left << tmpTensor.GetState()
						<< std::setw(10) << std::left << tmpTensor.GetK(j)
						<< std::setw(10) << std::left << tmpTensor.GetKappa(j)
						<< std::setw(14) << std::left << tmpTensor.GetTensor(j)
						<< std::endl;
			}
		}
	}
}

double	PointCoulEx::GetATS(int Ne){

	if ( Ne <= 0 || Ne > 96 ) {
    		return 0.0;
	}

	int m = Ne/2 + 1;
	if (Ne%2) {

		if ( 	m==1 || m==2 || m==3 || m==6 || 
 			m==7 || m==10 || m==15 || m==16 || 
 			m==19 || m==24 || m==25 || m==28 || 
 			m==31 || m==35 || m==37 || m==40 || 
 			m==41 || m==44 ) {
			return 0.5;
		}
		else if ( m==4 || m==5 || m==8 || m==9 || 
      			m==11 || m==17 || m==18 || m==20 ||
      			m==26 || m==27 || m==36 || m==42 ||
      			m==43 || m==45 ) {
			return 1.5;
		}
		else if ( m==12 || m==14 || m==21 || m==23 ||
      			m==32 || m==39 ) {
			return 2.5;
		}
		else if ( m==13 || m==22 || m==38 ) {
			return 4.5;
		}
		else if ( m==29 || m==30 || m==48 ) {
			return 3.5;
		}
		else if ( m==33 ) {
			return 7.5;
		}
		else if ( m==34 ) {
			return 6.5;
		}
		else if ( m==46 || m==47 ) {
			return 5.5;
		}
	}
   
	m -= 1;
	if ( 	m==4 || m==8 || m==17 || m==26 || 
		m==28 || m==30 || m==32 || m==42 || 
		m==45 || m==48 ) {
		return 2.0;
	}
	else if ( m==10 || m==36 ) {
		return 1.0;
	}
	else if ( m==12 || m==21 || m==37 ) {
		return 3.0;
	}
	else if ( m==13 || m==22 || m==29 || m==31 || 
	    	m==34 || m==38 || m==47 ) {
		return 4.0;
	}
	else if ( m==33 ) {
		return 8.0;
	}
	else if ( m==46 ) {
		return 6.0;
	}

	return 0.0;

}

void PointCoulEx::XSTATIC(int iz, double beta, int& id, int& iu, double& qcen, double& dq, double& xnor) {

	double h = 1.0/(1.0 + TMath::Power(TMath::Power(iz,0.45)*0.012008/beta,5.0/3.0));
	qcen = iz*TMath::Power(h,0.6);
	dq = TMath::Sqrt(qcen*(1.0-h))/2.0;

	iu = int(qcen + 3.0*dq + 0.5);
	id = int(qcen - 3.0*dq - 0.5);

	if(iu > iz) {
		iu = iz;
	}
	if(id < 1) {
		id = 1;
	}

	xnor = 0.0;
	for(int i=id;i<iu+1;i++) {
		xnor += TMath::Exp(-TMath::Power((qcen-i)/dq,2)/2.0);
	}

	return;
  
}

std::array<double,7> PointCoulEx::GKK(const int iz, const int ia, const double beta, const double spin, const double time) {
  
	//model parameters
	const double Avji = 3.0; 				//Average atomic spin
	const double Gam = 0.02; 				//FWHM of frequency distribution (ps^-1)
	const double Xlamb = 0.0345; 				//Fluctuating state to static state transition rate (ps^-1)
	const double TimeC = 3.5; 				//Mean time between random reorientations of fluctuating state  (ps)
	const double Gfac = iz/double(ia); 			//Nuclear gyromagnetic factor
	const double Field = 6.0*TMath::Power(10.0,-6.0); 	//Hyperfine field coefficient (600 T)
	const double Power = 0.6; 				//Hyperfine field exponent

	int  inq, ifq;
	double qcen, dq, xnor;
	XSTATIC(iz,beta,inq,ifq,qcen,dq,xnor);

	double aks[6] = {0.0,0.0,0.0,0.0,0.0,0.0}; 		//alpha_k
	double sum[3] = {0.0,0.0,0.0}; 				//stores sum over 6j symbols

	for(int j=inq;j<ifq+1;j++) {
		int nz = iz - j;
		double xji = GetATS(nz);
		double sm = spin;
		if(spin > xji) {
			sm = xji;
		}

		int ncoup = int(2.0*sm + 0.5) + 1;
		sum[0] = 0.0;
		sum[1] = 0.0;
		sum[2] = 0.0;
	
		double valmi = spin - xji;
		if(valmi < 0.0) {
			valmi *= -1;
		}
	
		for(int m=0;m<ncoup;m++) {
				double f = valmi + m;
			for(int k=0;k<3;k++) {
				double rk = 2.0*k + 2.0;
				int if2 = int(2.0*f + 0.0001);
				int irk2 = int(2.0*rk + 0.0001);
				int ispin2 = int(2.0*spin + 0.0001);
				int ixji2 = int(2.0*xji + 0.0001);
				sum[k] += TMath::Power((2.0*f + 1.0)*ROOT::Math::wigner_6j(if2,if2,irk2,ispin2,ispin2,ixji2),2.0)/(2.0*xji + 1.0);
			}
		}
	
		for(int k=0;k<3;k++) {
			int k1 = 2*k;
			aks[k1] += sum[k]*TMath::Exp(-TMath::Power((qcen-j)/dq,2)/2.0)/xnor;
		}

	} //end loop over j (charge states)

	double xji = Avji;
	double sm = spin;
	if(spin > xji) {
		sm = xji;
	}

	int ncoup = int(2.0*sm + 0.5) + 1;
	sum[0] = 0.0;
	sum[1] = 0.0;
	sum[2] = 0.0;

	double valmi = spin - xji;
	if(valmi < 0.0) {
		valmi *= -1;
	}

	for(int m=0;m<ncoup;m++) {
		double f = valmi + m;
		for(int k=0;k<3;k++) {
			double rk = 2.0*k + 2.0;
			int if2 = int(2.0*f + 0.0001);
			int irk2 = int(2.0*rk + 0.0001);
			int ispin2 = int(2.0*spin + 0.0001);
			int ixji2 = int(2.0*xji + 0.0001);
			sum[k] += TMath::Power((2.0*f+1.0)*ROOT::Math::wigner_6j(if2,if2,irk2,ispin2,ispin2,ixji2),2.0)/(2.0*xji+1.0);
		}
	}

	for(int k=0;k<3;k++) {   
		int k1 = 2*k + 1;
		aks[k1] += sum[k];
	}

	double hmean = Field*iz*TMath::Power(beta,Power);
	double wsp = 4789.0*Gfac*hmean/Avji; // 4789 is mu_N
	wsp *= TimeC;
	wsp *= (wsp*Avji*(Avji+1.0)/3.0);

	std::array<double,7> Gk = {1.0,1.0,1.0,1.0,1.0,1.0,1.0};
	for(int k=0;k<3;k++) {

		int k2 = 2*k + 2;
		int k1 = 2*k + 1;

		double wrt = wsp*k2*(k2+1);
		double w2 = wrt;

		wrt *= (-1.0/(1.0 - aks[k2-1]));

		double xlam = (1.0 - aks[k2-1])*(1.0 - TMath::Exp(wrt))/TimeC;
		double up = (Gam*time*aks[k1-1] + 1.0)/(time*Gam + 1.0);
		up *= (Xlamb*time);
		up += 1.0;
	
		double down = time*(xlam+Xlamb) + 1.0;
		Gk.at(k2) = up/down;

		double alp = TMath::Sqrt(9.0*xlam*xlam + 8.0*xlam*TimeC*(w2 - xlam*xlam)) - 3.0*xlam;
		alp /= (4.0*xlam*TimeC);
	
		double upc = xlam*time*(down - 2.0*alp*alp*time*TimeC);
		double dwc = (down + alp*time)*(down + 2.0*alp*time);
		double ccf = 1.0 + upc/dwc;
		Gk.at(k2) *= ccf;

	}

	return Gk;

}

//************************************ TENSOR CALCULATIONS END ***************************************//

//************************************* PRINTING SECTION BEGINS **************************************//

//****************************************************************************************************//
//	Really simple printing function to output to the terminal some of the more pertient
//	calculation details.
//****************************************************************************************************//
void PointCoulEx::Print(){

	std::cout << "\nTheta [CM]: " << fTheta << " degrees" << std::endl;
	std::cout << "Theta [Lab]: " << fReaction.ConvertThetaCmToLab(fTheta * TMath::DegToRad(),2) * TMath::RadToDeg() << " degrees" << std::endl;
	std::cout << std::setw(20) << std::left << "Substate: " <<  std::setw(20) << std::left << "Substate prob." << std::left << std::setw(20) << "Real Amp." << std::left << std::setw(20) << "Imag. Amp." << std::endl;
	for(int i=0;i<fSubStateProbabilities.GetNrows();i++){
		std::cout 	<< std::setw(20) << std::left << fSubstates.at(i).GetM() 	
				<< std::setw(20) << std::left << fSubStateProbabilities[i] 
				<< std::setw(20) << std::left << FinalRealAmplitude[i][0] 
				<< std::setw(20) << std::left << FinalImagAmplitude[i][0] << std::endl;
	}
	std::cout << "\n State Probabilities:" << std::endl;
	for(int i=0;i<Probabilities.GetNrows();i++)
		std::cout << Probabilities[i] << std::endl;
	
}

//****************************************************************************************************//
//	Print detailed summary of the calculation to a file, contains information on coupling
//	constants, probabilities, amplitudes (real and imaginary) and more.
//****************************************************************************************************//
void PointCoulEx::WriteDetailsToFile(const char* outfilename){

        std::ofstream outfile;
        outfile.open(outfilename);	

	int nPart = 2;
	if(bTargetDetection)
		nPart = 3;

	outfile << std::setw(10) << std::left << "Proj. Z:" 
		<< std::setw(10) << std::left << "Proj. A:" 
		<< std::setw(10) << std::left << "Tar. Z:" 
		<< std::setw(10) << std::left << "Tar. A:"
		<< std::setw(14) << std::left << "E Lab [MeV]:"
		<< std::setw(14) << std::left << "E CM [MeV]:"
		<< std::setw(14) << std::left << "Rutherford:"
		<< std::endl;

	outfile << std::setw(10) << std::left << fReaction.GetBeamZ() 
		<< std::setw(10) << std::left << fReaction.GetBeamA() 
		<< std::setw(10) << std::left << fReaction.GetTargetZ()  
		<< std::setw(10) << std::left << fReaction.GetTargetA()  
		<< std::setw(14) << std::left << fReaction.GetLabEnergy() 
		<< std::setw(14) << std::left << fReaction.GetCMEnergy()
		<< std::setw(14) << std::left << fReaction.RutherfordCM(fTheta,nPart)
		<< std::endl;

	outfile << std::endl
		<< "Accuracy: " << fAccuracy
		<< std::endl;

	outfile << "\nTheta [lab]: " << fReaction.ConvertThetaCmToLab(fTheta * TMath::DegToRad(),nPart) * TMath::RadToDeg() << " (degrees), Theta [CM]: " << fTheta  << " (degrees)" << std::endl;

	outfile << "\nClosest separation distance (fm): " << fReaction.ClosestApproach() << std::endl;

	outfile << "\nSommerfield Parameter (g.s.): " << fReaction.EtaCalc(0) << std::endl;

	outfile << "\nElectric dipole polarization factor: " << ElectricDipoleNormalization() << std::endl;

	outfile << std::endl;

	std::string mult[7] = {"E1","E2","E3","E4","E5","E6","M1"};

	for(int i=0;i<fNucleus.GetMaxLambda();i++){
		if(MiscFunctions::GetMaxAbsMatrix(fNucleus.GetMatrixElements().at(i)) == 0)
			continue;
		outfile 	<< mult[i]
				<< std::endl;
		MiscFunctions::WriteMatrixNucleus(outfile,fNucleus.GetMatrixElements().at(i),fNucleus);
	}

	outfile << std::endl;

	WriteConnections(outfile);

	outfile << "Final amplitudes:"
		<< std::endl;
	outfile << std::setw(10) << std::left << "SPIN:" 
		<< std::setw(10) << std::left << "M:" 
		<< std::setw(15) << std::left << "REAL AMP.:" 
		<< std::setw(15) << std::left << "IMAG AMP.:"
		<< std::endl;
	for(size_t ss = 0; ss < fSubstates.size(); ss++){
		outfile << std::setw(10) << std::left << fNucleus.GetLevelJ().at(fSubstates.at(ss).GetStateIndex())  
			<< std::setw(10) << std::left << fSubstates.at(ss).GetM()
			<< std::setw(15) << std::left << FinalRealAmplitude[ss][0] 
			<< std::setw(15) << std::left << FinalImagAmplitude[ss][0] << "\n";
	}

	double Amplitudes[fSubstates.size()];
	for(size_t j=0;j<fSubstates.size();j++)
		Amplitudes[j] = fSubStateProbabilities[j];
	

	outfile << "\nSubstate probabilities:\n\n";
	for(size_t ss = 0; ss < fSubstates.size(); ss++){	
		outfile 	<< std::setw(6) << fNucleus.GetLevelJ().at(fSubstates.at(ss).GetStateIndex()) 
				<< std::setw(6) << fSubstates.at(ss).GetM() 
				<< std::setw(15) << Amplitudes[ss] 
				<< std::endl;
	}

	outfile << "\n\n" 
		<< std::setw(8) << std::left << "State:" 
		<< std::setw(16) << std::left << "Prob.:"
		<< std::setw(16) << std::left << "CS:"
		<< std::endl;
	for(int i=0;i<fNucleus.GetNstates();i++)
		outfile << std::setw(8) << std::left << fNucleus.GetLevelJ().at(i) 
			<< std::setw(16) << std::left << Probabilities[i] 
			<< std::setw(16) << std::left << Probabilities[i] * fReaction.RutherfordCM(fTheta,nPart) * TMath::Sin(fReaction.ConvertThetaCmToLab(fTheta * TMath::DegToRad(),2)) << "\n";

	outfile		<< std::endl
			<< "STATISTICAL TENSORS:" 
			<< std::endl;

	outfile		<< std::endl
			<< "EXCITATION FRAME:" 
			<< std::endl;

	outfile		<< std::setw(10) << std::left << "INDEX:"
			<< std::setw(10) << std::left << "KA:"
			<< std::setw(10) << std::left << "KAPPA:"
			<< std::setw(14) << std::left << "RHOB:"
			<< std::endl;
	for(int i=0;i<fTensorsB.GetNstates();i++){	
		StateTensor tmpTensor = fTensorsB.GetStateTensor(i);
		for(int j=0;j<tmpTensor.GetNelements();j++){
			outfile 	<< std::setw(10) << std::left << tmpTensor.GetState()
					<< std::setw(10) << std::left << tmpTensor.GetK(j)
					<< std::setw(10) << std::left << tmpTensor.GetKappa(j)
					<< std::setw(14) << std::left << tmpTensor.GetTensor(j)
					<< std::endl;
		}
	}

	outfile		<< std::endl
			<< "LAB FRAME:" 
			<< std::endl;

	outfile		<< std::setw(10) << std::left << "INDEX:"
			<< std::setw(10) << std::left << "KA:"
			<< std::setw(10) << std::left << "KAPPA:"
			<< std::setw(14) << std::left << "RHOB:"
			<< std::endl;
	for(int i=0;i<fTensors.GetNstates();i++){	
		StateTensor tmpTensor = fTensors.GetStateTensor(i);
		for(int j=0;j<tmpTensor.GetNelements();j++){
			outfile 	<< std::setw(10) << std::left << tmpTensor.GetState()
					<< std::setw(10) << std::left << tmpTensor.GetK(j)
					<< std::setw(10) << std::left << tmpTensor.GetKappa(j)
					<< std::setw(14) << std::left << tmpTensor.GetTensor(j)
					<< std::endl;
		}
	}
	outfile.close();

}

//****************************************************************************************************//
//	Sub-function specifically to write connection information to file in a vaguely
//	sensible format.	
//****************************************************************************************************//
void PointCoulEx::WriteConnections(std::ofstream& outfile){

	if(!outfile.is_open())
		return;

	std::string mult[7] = {"E1","E2","E3","E4","E5","E6","M1"};
	for(int i=0;i<fNucleus.GetMaxLambda();i++){

		if(MiscFunctions::GetMaxAbsMatrix(fNucleus.GetMatrixElements().at(i)) == 0)
			continue;

		outfile		<< "Multipolarity: " << mult[i] << std::endl;

		outfile 	<< "Substates: " << fSubstates.size() << std::endl;

		outfile 	<< std::setw(18) << std::right << "State Index F:"
				<< std::setw(15) << std::right << "State J F:"
				<< std::setw(15) << std::right << "Substate F:"	
				<< std::setw(18) << std::right << "State Index I:"
				<< std::setw(15) << std::right << "State I:"
				<< std::setw(15) << std::right << "Substate I:"
				<< std::setw(12) << std::right << "Xi:"
				<< std::setw(12) << std::right << "Psi:"
				<< std::setw(12) << std::right << "Zeta:"
				<< std::endl;

		int counter = 0;

		for(size_t ss = 0; ss < fSubstates.size(); ss++){
			bool first = true;
			for(size_t ss2 = 0; ss2 < fSubstates.at(ss).GetNconnections(); ss2++){
				counter++;
				int substate 	= fSubstates.at(ss).GetConnection(ss2).GetConnectedState();
				if(first){
					outfile		<< std::setw(18) << std::right << fSubstates.at(ss).GetStateIndex()
							<< std::setw(15) << std::right << fStates.at(fSubstates.at(ss).GetStateIndex()).GetJ()
							<< std::setw(15) << std::right << fSubstates.at(ss).GetM()
							<< std::setw(18) << std::right << fSubstates.at(substate).GetStateIndex()
							<< std::setw(15) << std::right << fStates.at(fSubstates.at(substate).GetStateIndex()).GetJ()
							<< std::setw(15) << std::right << fSubstates.at(substate).GetM()
							<< std::setw(12) << std::right << fSubstates.at(ss).GetConnection(ss2).GetXi(i)
							<< std::setw(12) << std::right << fSubstates.at(ss).GetConnection(ss2).GetPsi(i)
							<< std::setw(12) << std::right << fSubstates.at(ss).GetConnection(ss2).GetZeta(i)	
							<< std::endl;
					first 		= false;
				}
				else{
					outfile		<< std::setw(18) << std::right << " "
							<< std::setw(15) << std::right << " "
							<< std::setw(15) << std::right << " "
							<< std::setw(18) << std::right << fSubstates.at(substate).GetStateIndex()
							<< std::setw(15) << std::right << fStates.at(fSubstates.at(substate).GetStateIndex()).GetJ()
							<< std::setw(15) << std::right << fSubstates.at(substate).GetM()
							<< std::setw(12) << std::right << fSubstates.at(ss).GetConnection(ss2).GetXi(i)
							<< std::setw(12) << std::right << fSubstates.at(ss).GetConnection(ss2).GetPsi(i)
							<< std::setw(12) << std::right << fSubstates.at(ss).GetConnection(ss2).GetZeta(i)
							<< std::endl;
				}
			}
		}

	}

}

