#include "ParticleDetectorS3.h"

ParticleDetectorS3::ParticleDetectorS3(int id, double z, double x_off, double y_off, int ring_in, int ring_out, double t_min, double t_max){

	char hname[128];
	sprintf(hname,"Expt_%i_ThetaPhi_S3_Ring_%i_to_Ring_%i",id,ring_in,ring_out);
	hThetaPhi = new TH2F(hname,hname,720,-180,180,3600,t_min,t_max);	

	double	radius1	= 10. + ring_in;
	double 	radius2	= 11. + ring_out;

	TGraph *g1 = new TGraph();
        TGraph *g2 = new TGraph();
        for(int i=0;i<360;i++){

                double  phi     = (double)i-180.;
                double  phi_r   = phi * TMath::DegToRad();
                double  theta   = TMath::ATan2(radius1,z) * TMath::RadToDeg();
                double  theta_r = TMath::ATan2(radius1,z);

                double  r       = z / TMath::Cos(theta_r);
                double  x       = r * TMath::Sin(theta_r) * TMath::Cos(phi_r) + x_off;
                double  y       = r * TMath::Sin(theta_r) * TMath::Sin(phi_r) + y_off;
                r               = TMath::Sqrt(TMath::Power(z,2)+TMath::Power(x,2)+TMath::Power(y,2));


                theta_r         = TMath::ACos(z/r);
                phi_r           = TMath::ATan2(y,x);
                theta           = theta_r * TMath::RadToDeg();
                phi             = phi_r * TMath::RadToDeg();

                g1->SetPoint(i,phi,theta);

        }
        for(int i=0;i<360;i++){

                double  phi     = (double)i-180.;
                double  phi_r   = phi * TMath::DegToRad();
                double  theta   = TMath::ATan2(radius2,z) * TMath::RadToDeg();
                double  theta_r = TMath::ATan2(radius2,z);

                double  r       = z / TMath::Cos(theta_r);
                double  x       = r * TMath::Sin(theta_r) * TMath::Cos(phi_r) + x_off;
                double  y       = r * TMath::Sin(theta_r) * TMath::Sin(phi_r) + y_off;
                r               = TMath::Sqrt(TMath::Power(z,2)+TMath::Power(x,2)+TMath::Power(y,2));


                theta_r         = TMath::ACos(z/r);
                phi_r           = TMath::ATan2(y,x);
                theta           = theta_r * TMath::RadToDeg();
                phi             = phi_r * TMath::RadToDeg();

                g2->SetPoint(i,phi,theta);

        }

	for(int i=1;i<=hThetaPhi->GetNbinsX();i++){
		double	phi = hThetaPhi->GetXaxis()->GetBinCenter(i);
		for(int j=1;j<=hThetaPhi->GetNbinsY();j++){
			double	theta = hThetaPhi->GetYaxis()->GetBinCenter(j);
			if(theta > g1->Eval(phi) && theta < g2->Eval(phi))
				hThetaPhi->SetBinContent(i,j,1/(double)hThetaPhi->GetNbinsX());
		}
	}

	sprintf(hname,"Expt_%i_ThetaEff_S3_Ring_%i_to_Ring_%i",id,ring_in,ring_out);
	hThetaEff = (TH1F*)hThetaPhi->ProjectionY(hname);

	setFlag	= true;

}
