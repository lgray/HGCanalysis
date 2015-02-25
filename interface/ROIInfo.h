#ifndef _roiinfo_h_
#define _roiinfo_h_

#include <iostream>

#include "DataFormats/Math/interface/deltaR.h"

#include "TMath.h"
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TVectorD.h>
#include "TVector2.h"

#define MAXDRINROI 8

/**
   @short indices are cone size, sub-detector or layer, energy estimator
 */
class ROIInfo : public TObject {
  
 public:
  Float_t center_eta, center_phi;
  Float_t dr[MAXDRINROI];
  Float_t en[MAXDRINROI][3][3],     et2[MAXDRINROI][3][3], eta[MAXDRINROI][3][3], phi[MAXDRINROI][3][3],  shh[MAXDRINROI][3][3],   shp[MAXDRINROI][3][3], spp[MAXDRINROI][3][3], width[MAXDRINROI][3][3], totalVolume[MAXDRINROI][3][3];
  Float_t nhits[MAXDRINROI][3][3];
  Float_t en_lay[MAXDRINROI][54][3],x[MAXDRINROI][54][3],  y[MAXDRINROI][54][3],  rho[MAXDRINROI][54][3], rho2[MAXDRINROI][54][3], area[MAXDRINROI][54][3];
  Float_t z[54];

  /**
     @short CTOR
   */
  ROIInfo();
  ROIInfo(float eta, float phi);
  ROIInfo(const ROIInfo &other);


  /**
     @short adds information
   */
  void add(float hit_en, float hit_x, float hit_y, float hit_z, float hit_eta, float hit_phi, float hit_weight, float hit_weight2, int hit_layer);

  /**
     @short sets information to 0
   */
  void reset();

  /**
     @short finalizes averages
   */
  void finalize();

  virtual ~ROIInfo() { }
   
  ClassDef(ROIInfo,1) 
};








#endif
