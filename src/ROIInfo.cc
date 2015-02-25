#include "UserCode/HGCanalysis/interface/ROIInfo.h"
#include "UserCode/HGCanalysis/interface/HGCAnalysisTools.h"

using namespace std;

ClassImp(ROIInfo)

//
ROIInfo::ROIInfo() { }

//
ROIInfo::ROIInfo(float eta, float phi) 
{ 
  reset(); 
  center_eta=eta; center_phi=phi; 
  dr[0]=0.05; dr[1]=0.1; dr[2]=0.15; dr[3]=0.2; dr[4]=0.25; dr[5]=0.3; dr[6]=0.4; dr[7]=0.5; 
}

//
ROIInfo::ROIInfo(const ROIInfo &other)
{
  center_eta=other.center_eta;
  center_phi=other.center_phi;
  for(size_t i=0; i<MAXDRINROI; i++)
    {
      dr[i]=other.dr[i];
      
      //sub-detector quantities
      for(size_t j=0; j<3; j++)
	{
	  for(size_t k=0; k<3; k++)	
	    {
	      en[i][j][k]=other.en[i][j][k]; 
	      nhits[i][j][k]=other.nhits[i][j][k]; 
	      nhits5mip[i][j][k]=other.nhits5mip[i][j][k]; 
	      et2[i][j][k]=other.et2[i][j][k];
	      eta[i][j][k]=other.eta[i][j][k];
	      phi[i][j][k]=other.phi[i][j][k];
	      shh[i][j][k]=other.shh[i][j][k]; 
	      shp[i][j][k]=other.shp[i][j][k]; 
	      spp[i][j][k]=other.spp[i][j][k];
	      width[i][j][k]=other.width[i][j][k];
	      totalVolume[i][j][k]=other.totalVolume[i][j][k];
	    }
	}
	  
      //layer quantities
      for(size_t j=0; j<54; j++)
	{
	  for(size_t k=0; k<3; k++)
	    {
	      en_lay[i][j][k] = other.en_lay[i][j][k];
	      x[i][k][j]      = other.x[i][j][k]; 
	      y[i][k][j]      = other.y[i][j][k]; 
	      if(i==0 && j==0) z[k] = other.z[k];
	      rho[i][k][j]  = other.rho[i][j][k]; 
	      rho2[i][k][j] = other.rho2[i][j][k]; 
	      area[i][k][j] = other.area[i][j][k];
	    }
	}
    }
}

//
void ROIInfo::add(float hit_en, float hit_x, float hit_y, float hit_z, float hit_eta, float hit_phi, float hit_weight, float hit_weight2, int hit_layer)
{
  {

    //distance to center eta-phi
    float hit_dR=deltaR(hit_eta,hit_phi,center_eta,center_phi);
    float hit_deta=hit_eta-center_eta;
    float hit_dphi=deltaPhi(hit_phi,center_phi);
 
    //distance to center cartesian
    float refRho=TMath::Abs(hit_z/TMath::SinH(center_eta));
    float hit_dx=hit_x-refRho*TMath::Cos(center_phi);
    float hit_dy=hit_y-refRho*TMath::Sin(center_phi);
    float hit_rho=sqrt(hit_dx*hit_dx+hit_dy*hit_dy);
    
    //hardcoded scale factors
    size_t subDet(0);
    float m(0),k(0),mu(0); 
    if (hit_layer==1)                       { subDet=0; m=0.2339; k=0.1778; mu=1.0;            }
    else if(hit_layer>=2 && hit_layer<=11)  { subDet=0; m=0.2339; k=0.1778; mu=1.0;            }
    else if(hit_layer>=12 && hit_layer<=21) { subDet=0; m=0.2339; k=0.1778; mu=1.0;            }
    else if(hit_layer>=22 && hit_layer<=30) { subDet=0; m=0.2339; k=0.1778; mu=1.0;            }
    else if(hit_layer==31)                  { subDet=1; m=0.1828; k=0.9601; mu=1.0535;         }
    else if(hit_layer>=32 && hit_layer<=42) { subDet=1; m=0.1828; k=0.9601; mu=1.0535;         }
    else if(hit_layer>42)                   { subDet=2; m=0.2471; k=1.4563; mu=1.29208*1.0535; }
    
    //convert from MIP to e.m. scale
    float em_en = (hit_en*m+k)*mu*getLambdaForHGCLayer(hit_layer);
    
    for(size_t i=0; i<MAXDRINROI; i++)
      {
	//neglect if out of cone
	if(hit_dR>dr[i]) continue;

	//construct sums depending on energy estimator
	for(size_t k=0; k<3; k++)
	  {
	    float theWeight(1.0);
	    if(k==1) theWeight=hit_weight;
	    if(k==2) theWeight=hit_weight2;
	    float en_k=em_en*theWeight;
	    float et_k(en_k/TMath::CosH(hit_eta));
	    en[i][subDet][k]          += en_k;
	    nhits[i][subDet][k]++;
	    en_lay[i][hit_layer-1][k] += en_k;
	    eta[i][subDet][k]         += en_k*hit_eta;
	    phi[i][subDet][k]         += en_k*TVector2::Phi_mpi_pi(hit_phi);
	    et2[i][subDet][k]         += pow(et_k,2);
	    shh[i][subDet][k]         += pow(et_k*hit_deta,2);
	    shp[i][subDet][k]         += -pow(et_k,2)*hit_deta*hit_dphi;
	    spp[i][subDet][k]         += pow(et_k*hit_dphi,2);
	    x[i][hit_layer-1][k]      += en_k*hit_x;
	    y[i][hit_layer-1][k]      += en_k*hit_y;
	    z[hit_layer-1]             = hit_z;
	    if(en_k>5)
	      {
		nhits5mip[i][subDet][k]++;
		rho[i][hit_layer-1][k]    += hit_rho;
		rho2[i][hit_layer-1][k]   += pow(hit_rho,2);
	      }
	  }
      }
  }
}

//
void ROIInfo::reset()
{
  for(size_t i=0; i<MAXDRINROI; i++)
    {
      //sub-detector quantities
      for(size_t j=0; j<3; j++)
	{
	  for(size_t k=0; k<3; k++)	
	    {
	      en[i][j][k]=0; 
	      nhits[i][j][k]=0; 
	      nhits5mip[i][j][k]=0; 
	      et2[i][j][k]=0;
	      eta[i][j][k]=0;
	      phi[i][j][k]=0;
	      shh[i][j][k]=0; 
	      shp[i][j][k]=0; 
	      spp[i][j][k]=0;
	      width[i][j][k]=0;
	      totalVolume[i][j][k]=0;
	    }
	}

      //layer quantities
      for(size_t j=0; j<54; j++)
	{
	  if(i==0) z[j]=0;
	    
	  for(size_t k=0; k<3; k++)
	    {
	      en_lay[i][j][k]=0;
	      x[i][j][k]=0; 
	      y[i][j][k]=0; 
	      rho[i][j][k]=0; 
	      rho2[i][j][k]=0; 
	      area[i][j][k]=0;
	    }
	}
    }
}

//
void ROIInfo::finalize()
{
  for(size_t i=0; i<MAXDRINROI; i++)
    {
      //sub-detector level
      for(size_t j=0; j<3; j++)
	{
	  for(size_t k=0; k<3; k++)
	    {
	      if(en[i][j][k]<=0) continue;
	      eta[i][j][k]/=en[i][j][k];
	      phi[i][j][k]/=en[i][j][k];

	      if(et2[i][j][k]>0)
		{
		  shh[i][j][k]/=et2[i][j][k];
		  shp[i][j][k]/=et2[i][j][k];
		  spp[i][j][k]/=et2[i][j][k];

		  double mvals[MAXDRINROI] = {
		    shh[i][j][k], shp[i][j][k],
		    shp[i][j][k], spp[i][j][k]
		  };
		    
		  TMatrixDSym m(2,mvals);
		  TMatrixDSymEigen me(m);
		  TVectorD eigenval = me.GetEigenValues();
		  width[i][j][k] = sqrt(pow(eigenval[0],2)+pow(eigenval[1],2));
		  if(width[i][j][k]>dr[i])
		    {
		      m.Print();
		      std::cout << et2[i][j][k] << std::endl;
		      std::cout << eigenval[0] << " " << eigenval[1] << " " << width[i][j][k] << std::endl;
		    }
		}
	    }
	}
	
      //layer level
      for(size_t j=0; j<54; j++)
	for(size_t k=0; k<3; k++)
	  {
	    if(en_lay[i][j][k]<=0) continue;
	    x[i][j][k]    /= en_lay[i][j][k];
	    y[i][j][k]    /= en_lay[i][j][k];
	    if(nhits5mip[i][j][k]>0){
	      rho[i][j][k]  /= nhits5mip[i][j][k];
	      rho2[i][j][k] /= nhits5mip[i][j][k];
	    }
	    float sigma=rho2[i][j][k]-pow(rho[i][j][k],2);
	    //if only one hit this may happen
	    if(sigma<0) {
	      sigma=0.0;
	      //std::cout << "Negative sigma found, set to 0" << std::endl;
	    }
	    area[i][j][k] = TMath::Pi()*sigma;

	    int subDet(0);
	    if(j>30) subDet=1;
	    if(j>42) subDet=2;
	    totalVolume[i][subDet][k] += area[i][j][k]*getLambdaForHGCLayer(j+1)*fabs(TMath::TanH(center_eta));
	  }
    }
}
