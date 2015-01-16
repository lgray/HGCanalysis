#ifndef _roiinfo_h_
#define _roiinfo_h_

#include "DataFormats/Math/interface/deltaR.h"

#include "TMath.h"
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TVectorD.h>

/**
   @short indices are cone size, sub-detector or layer, energy estimator
 */
class ROIInfo
{
 public:

  float center_eta, center_phi;

  float dr[4];
  float en[4][3][2],     et2[4][3][2], eta[4][3][2], phi[4][3][2],  shh[4][3][2],   shp[4][3][2], spp[4][3][2], width[4][3][2], totalVolume[4][3][2];
  float en_lay[4][54][2],x[4][54][2],  y[4][54][2],  rho[4][54][2], rho2[4][54][2], area[4][54][2];
  float z[54];


  /**
     @short CTOR
   */
  ROIInfo(float eta, float phi) { reset(); center_eta=eta; center_phi=phi; dr[0]=0.1; dr[1]=0.2; dr[2]=0.3; dr[3]=0.4; }
  ROIInfo(const ROIInfo &other)
    {
      center_eta=other.center_eta;
      center_phi=other.center_phi;
      dr[0]=other.dr[0]; dr[1]=other.dr[1]; dr[2]=other.dr[2]; dr[3]=other.dr[3];
      for(size_t i=0; i<4; i++)
	{
	  //sub-detector quantities
	  for(size_t j=0; j<3; j++)
	    {
	      for(size_t k=0; k<2; k++)	
		{
		  en[i][j][k] =other.en[i][j][k]; 
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
	      for(size_t k=0; k<2; k++)
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

  float getLambda(int hit_layer)
  {
    if (hit_layer==1)                       return 0.01; 
    else if(hit_layer>=2 && hit_layer<=11)  return 0.036;
    else if(hit_layer>=12 && hit_layer<=21) return 0.043;
    else if(hit_layer>=22 && hit_layer<=30) return 0.056;
    else if(hit_layer==31)                  return 0.338;
    else if(hit_layer>=32 && hit_layer<=42) return 0.273;
    else if(hit_layer>42)                   return 0.475;
    return 0;
  }

  /**
     @short adds information
   */
  void add(float hit_en, float hit_x, float hit_y, float hit_z, float hit_eta, float hit_phi, float hit_weight, int hit_layer)
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
    else if(hit_layer==31)                  { subDet=1; m=0.1828; k=0.9601; mu=1.29208;        }
    else if(hit_layer>=32 && hit_layer<=42) { subDet=1; m=0.1828; k=0.9601; mu=1.29208;        }
    else if(hit_layer>42)                   { subDet=2; m=0.2471; k=1.4563; mu=1.29208*1.0535; }
    
    //convert from MIP to e.m. scale
    float em_en = (hit_en*m+k)*mu*getLambda(hit_layer);
    
    for(size_t i=0; i<4; i++)
      {
	//neglect if out of cone
	if(hit_dR>dr[i]) continue;

	//construct sums depending on energy estimator
	for(size_t k=0; k<2; k++)
	  {
	    float en_k=em_en*(k==0 ? 1.0 : hit_weight);
	    float et_k(en_k/TMath::CosH(hit_eta));
	    en[i][subDet][k]      += en_k;
	    en_lay[i][hit_layer-1][k] += en_k;
	    eta[i][subDet][k]     += en_k*hit_eta;
	    phi[i][subDet][k]     += en_k*hit_phi;
	    et2[i][subDet][k]     += pow(et_k,2);
	    shh[i][subDet][k]     += pow(en_k*hit_deta,2);
	    shp[i][subDet][k]     += -pow(en_k,2)*hit_deta*hit_dphi;
	    spp[i][subDet][k]     += pow(en_k*hit_dphi,2);
	    x[i][hit_layer-1][k]      += en_k*hit_x;
	    y[i][hit_layer-1][k]      += en_k*hit_y;
	    z[hit_layer-1]             = hit_z;
	    rho[i][hit_layer-1][k]    += en_k*hit_rho;
	    rho2[i][hit_layer-1][k]   += en_k*pow(hit_rho,2);
	  }
      }
  }

  /**
     @short sets information to 0
   */
  void reset()
  {
    for(size_t i=0; i<4; i++)
      {
	//sub-detector quantities
	for(size_t j=0; j<3; j++)
	  {
	    for(size_t k=0; k<2; k++)	
	      {
		en[i][j][k]=0; 
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
	    
	    for(size_t k=0; k<2; k++)
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

  /**
     @short finalizes averages
   */
  void finalize()
  {
    for(size_t i=0; i<4; i++)
      {
	//sub-detector level
	for(size_t j=0; j<3; j++)
	  {
	    for(size_t k=0; k<2; k++)
	      {
		if(en[i][j][k]<=0) continue;
		eta[i][j][k]/=en[i][j][k];
		phi[i][j][k]/=en[i][j][k];
	   
		if(et2[i][j][k]>0)
		  {
		    double mvals[4] = {
		      shh[i][j][k], shp[i][j][k],
		      shp[i][j][k], spp[i][j][k]
		    };
		    
		    TMatrixDSym m(2,mvals);
		    TMatrixDSymEigen me(m);
		    TVectorD eigenval = me.GetEigenValues();
		    width[i][j][k] = sqrt((pow(eigenval[0],2)+pow(eigenval[1],2))/et2[i][j][k]);
		  }
	      }
	  }
	
	//layer level
	for(size_t j=0; j<54; j++)
	  for(size_t k=0; k<2; k++)
	    {
	      if(en_lay[i][j][k]<=0) continue;
	      x[i][j][k]    /= en_lay[i][j][k];
	      y[i][j][k]    /= en_lay[i][j][k];
	      rho[i][j][k]  /= en_lay[i][j][k];
	      rho2[i][j][k] /= en_lay[i][j][k];
	      area[i][j][k] = TMath::Pi()*(rho2[i][j][k]-pow(rho[i][j][k],2));

	      int subDet(0);
	      if(j>30) subDet=1;
	      if(j>42) subDet=2;
	      totalVolume[i][subDet][k] += area[i][j][k]*getLambda(j+1)*fabs(TMath::TanH(center_eta));
	    }
      }
  }

  ~ROIInfo() {}


};


#endif
