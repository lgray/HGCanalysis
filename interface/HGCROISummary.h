#ifndef _hgcroisummary_h_
#define _hgcroisummary_h_

#include "TTree.h"

struct  HGCROISummary_t
{ 
  //event identifier
  Int_t run,lumi,event;

  //generator level ROI
  Int_t gen_id;
  Float_t gen_pt, gen_eta, gen_phi, gen_en, gen_emfrac, gen_hadfrac, gen_invfrac;
  
  float en[4][3], eta[4][3], phi[4][3], width[4][3], totalVolume[4][3];
  float x[4][54], y[4][54];
  float z[54];

  float wgt_en[4][3], wgt_eta[4][3], wgt_phi[4][3], wgt_width[4][3], wgt_totalVolume[4][3];
  float wgt_x[4][54], wgt_y[4][54];
};


void initHGCROITree(TTree *t,HGCROISummary_t &simEvt);
void attachHGCROITree(TTree *t,HGCROISummary_t &simEvt);

#endif
