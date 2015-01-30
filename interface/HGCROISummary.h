#ifndef _hgcroisummary_h_
#define _hgcroisummary_h_

#include "TTree.h"

struct  HGCROISummary_t
{ 
  //event identifier
  Int_t run,lumi,event;

  Int_t ncandsc;
  
  //roi information
  Int_t roi_nsc;
  Float_t roi_eta, roi_phi;

  //generator level ROI
  Int_t gen_id;
  Float_t gen_pt, gen_eta, gen_phi, gen_en, gen_emfrac, gen_hadfrac, gen_invfrac;
  
  float en[6][3], eta[6][3], phi[6][3], width[6][3], totalVolume[6][3];
  int nhits[6][3];
  float x[6][54], y[6][54];
  float z[54];

  float wgt_en[6][3], wgt_eta[6][3], wgt_phi[6][3], wgt_width[6][3], wgt_totalVolume[6][3];
  int wgt_nhits[6][3];
  float wgt_x[6][54], wgt_y[6][54];
};


void initHGCROITree(TTree *t,HGCROISummary_t &simEvt);
void attachHGCROITree(TTree *t,HGCROISummary_t &simEvt);

#endif
