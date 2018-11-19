CYGNUS
=============================
A CODE FOR COULOMB-EXCITATION ANALYSIS IN C++
==================================================================

Note that CYGNUS is still in a testing phase and should not be used for analysis in isolation, but should be used alongside a traditional code such as GOSIA. 

Note that kinematics calculations differ from those in GOSIA. This difference is typically less than 0.1 degrees and results in a difference in excitation probability of less than 0.1%. Typically this will lie beyond experimental sensitivities, however it should be considered when comparing to GOSIA calculations.

Requirements
-----------------------

- Cygnus requires C++11 and ROOT v6 with the MathMore libraries
- In order to run compiled scripts, the Cygnus/bin directory must be added to LD_LIBRARY_PATH

Running Cygnus
-----------------------

In a ROOT terminal:
- Run rootstart.C on root startup (loads libraries)

A compiled code:
- An example (Fitting.cxx) and makefile are included in the scripts directory

Data format
------------------------

Cygnus ignores any line containing an exclamation mark (!)

Nucleus file (an example is included):

1st line: A, Z, NStates, MaxLambda
- A and Z are mass number and proton number of the nucleus
- NStates is the number of states in the nucleus
- MaxLambda is the maximum lambda (1-6 = E1-E6, 7 = M1) recommended to leave MaxLambda = 7

Lines 2 -> 2 + NStates: Index, E, J, P
- Index is the state index, ordered by energy (0 = ground state)
- E is the energy, in units of MeV
- J and P are the spin and parity of the state

Lines 2 + NStates -> End: Index_I, Index_F, ME, Lambda
- Index_I is the initial state (by convention, lower energy)
- Index_F is the final state (higher energy) 
- ME is the matrix element in eb units (note: Not efm)
- Lambda is the multipolarity of the transition (same convention as previously stated) 

Data file (an example is included):

For every experiment:

Line 1: EXPT, Index, Energy, Theta_min, Theta_max
- EXPT (text: "EXPT")
- Index is the experiment number
- Energy is the beam energy - not presently used
- Theta_min is the minimum theta - not presently used
- Theta_max is the maximum theta - not presently used

Line 2 -> 2+nData: Index_I, Index_F, Yield, Uncertainty 
- Index_I is initial state. Note that now, the initial state is the higher energy state.
- Index_F is the final state. Note that now, the final state is the lower energy state.
- Yield is the experimental counts
- Uncertainty is the uncertainty on the experimental data

Basic Usage
------------------------

Class: NucleusReader

Loads data from a nucleus file (formatted as above) and puts it into a Nucleus format.

Class: Nucleus

Class holding the nucleus information, as loaded from the nucleus data file. Can be grabbed from NucleusReader:
> Nucleus *nucl = NucleusReader->GetNucleus();

Class: Reaction

Deals with the reaction kinematics of the experiment:
> Reaction *reac = new Reaction(aB,zB,aT,zT,eBeam)
- aB and zB are mass and proton number of the beam
- aT and zT are mass and proton number of the target
- eBeam is the beam energy in the lab (MeV)

Can set actual nuclear masses using:
> reac->SetMasses(mB,mT)
- mB and mT are beam and target masses (amu)

Class: PointCoulEx

Performs a Coulomb excitation calculation for a single energy at a single theta value. Requires a Nucleus and Reaction:
> PointCoulEx *poin = new PointCoulEx(nucl,reac)

To run the calculation:
> poin->CalculatePointProbabilities(theta)
- theta in the center of mass frame, in units of degrees

To get the excitation probabilities:
> TVectorD vec = point->GetProbabilitiesVector()

Class: StoppingPower

Holder class for stopping power information. Performs a cubic fit to stopping power data. 
> StoppingPower dEdX; 
> dEdX.AddStoppingPower(E,SP)
- E and SP are the energy (MeV) and stopping power (MeV/mg/cm2) meshpoint
> dEdX.FitStoppingPowers()
- Fits stopping power data with cubic polynomial

Class: Experiments

Performs integrations of Coulomb excitation values and determines the correction for a point calculation to a full integral over theta and energy. Requires Nucleus and Reaction:
> Experiments *expt = new Experiments(nucl,reac);

To define a new experimental range:
> expts->NewExperimentRange(tmin,tmax,nt,emin,emax,ne,tarDet)
- tmin and tmax are the min and max theta in degrees, in the lab frame
- nt is the number of theta meshpoints to use
- emin and emax are th min and max energies, in MeV
- ne is the number of energy meshpoints to use
- tarDet is a bool which is TRUE for target detection and FALSE for projectile detection
> expts->SetStopping(dEdX)
- THIS IS A NECESSARY STEP FOR THE INTEGRATION PROCESS TO WORK - without the stopping powers, the integration is invalid
> expts->PointCorrections
- Determines the point correction factors

Fitting Data
------------------------

Cygnus uses the ROOT::Minimizer package. This comes with a number of minimization options. The user can define a minimizer choice and select an algorithm for that minimizer. The options are:

| MINIMIZER: |	ALGORITHM:      |
|------------|------------------|
|Minuit2     |	Migrad          |
|	     |	Simplex         |
|	     |	Combined        |
|	     |	Scan            |
|	     |	Fumili2         |
|Genetic     |                  -	
|GSLMultiMin |	ConjugateFR     |
|	     |	ConjugatePR     |
|            |	BFGS            |
|            |	BFGS2           |
|            |	SteepestDescent |
|GSLMultiFit |                  |
|GSLSimAn    |	                |

License
-----------------------

CYGNUS is distributed under the GPL v3 license.

``LLNL-CODE-761182``
