void MacroExample(){

	NucleusReader *reader_mg = new NucleusReader("scripts/NucleusFile_Mg22.txt");
	Nucleus *nucl_mg = reader_mg->GetNucleus();
	Reaction *reac_mg = new Reaction(22,12,110,46,(double)78.);
	reac_mg->SetMass(21.999574,109.90515);
	reac_mg->SetGOSIAKinematics();
	reac_mg->SetExcitationEnergy(1.274);

	PointCoulEx *point_mg = new PointCoulEx(nucl_mg,reac_mg);
	point_mg->SetProjectileExcitation(true);
	point_mg->SetAccuracy(1e-8);
	point_mg->CalculatePointProbabilities(reac_mg->ConvertThetaLabToCm(26.2270 * TMath::DegToRad(),2) * TMath::RadToDeg());

	point_mg->Print();
	
}
