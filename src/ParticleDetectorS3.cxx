#include "ParticleDetectorS3.h"

ParticleDetectorS3::ParticleDetectorS3(TGraph2D* g2, TGraph* g){

	SetThetaPhiMap(g2);
	SetThetaEfficiencyGraph(g);
	setFlag = true;	

}
 
ParticleDetectorS3::ParticleDetectorS3(int id, double z, double x_off, double y_off, int ring_in, int ring_out, double t_min, double t_max){

	char gname[128];
	sprintf(gname,"Expt_%i_ThetaPhi_S3_Ring_%i_to_Ring_%i",id,ring_in,ring_out);
	gThetaPhi = new TGraph2D();
	gThetaPhi->SetName(gname);	

	double	radius1	= 10. + ring_in;
	double 	radius2	= 11. + ring_out;

	// Define two TGraphs which define the limits of coverage in theta-phi space
	// to be used to construct the full theta-phi
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

	gThetaEff = new TGraph();
	int c = 0;
	double t_step = (t_max - t_min)/2000.;
	for(int i=0;i<2001;i++){
		double theta = t_min + i*t_step;
		double integral = 0;
		for(int j=0;j<360;j++){
			double phi = (double)j/10.-180.;
			if(theta > g1->Eval(phi) && theta < g2->Eval(phi)){
				gThetaPhi->SetPoint(c,theta,phi,1);
				integral += 1;
			}
			else
				gThetaPhi->SetPoint(c,theta,phi,0);
			c++;
		}
		gThetaEff->SetPoint(i,theta,integral/360.);	
	}
	sprintf(gname,"Expt_%i_ThetaEff_S3_Ring_%i_to_Ring_%i",id,ring_in,ring_out);
	gThetaEff->SetName(gname);

	setFlag	= true;

}
