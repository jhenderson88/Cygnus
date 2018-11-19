void MacroExample(){

	NucleusReader *reader_mg = new NucleusReader("NucleusFile_Mg22.txt");
	Nucleus *nucl_mg = reader_mg->GetNucleus();
	NucleusReader *reader_pd = new NucleusReader("NucleusFile_Pd110.txt");
	Nucleus *nucl_pd = reader_pd->GetNucleus();
	Reaction *reac_mg = new Reaction(22,12,110,46,(double)78.);
	reac_mg->SetMass(21.999574,109.90515);
	Reaction *reac_pd = new Reaction(22,12,110,46,(double)78.);
	reac_pd->SetMass(21.999574,109.90515);
	reac_mg->SetExcitationEnergy(1.274);
	reac_pd->SetExcitationEnergy(0.3738);

	PointCoulEx *point_mg = new PointCoulEx(nucl_mg,reac_mg);
	point_mg->SetProjectileExcitation(true);
	point_mg->CalculatePointProbabilities(reac_mg->ConvertThetaLabToCm(130 * TMath::DegToRad(),2) * TMath::RadToDeg());
	point_mg->SetAccuracy(1e-7);
	PointCoulEx *point_pd = new PointCoulEx(nucl_pd,reac_pd);
	point_pd->SetProjectileExcitation(false);
	point_pd->CalculatePointProbabilities(reac_mg->ConvertThetaLabToCm(130 * TMath::DegToRad(),2) * TMath::RadToDeg());
	point_pd->SetAccuracy(1e-7);

	point_mg->WriteDetailsToFile("CoulExCalc_mg.txt");
	point_pd->WriteDetailsToFile("CoulExCalc_pd.txt");
}
