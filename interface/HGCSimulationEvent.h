#ifndef _hgcsimulationevent_h_
#define _hgcsimulationevent_h_

#include "TTree.h"

#define MAXGENPEREVENT 100
#define MAXHGCHITSPEREVENT 1000000

struct  HGCSimEvent_t
{ 
  //event identifier
  Int_t run,lumi,event;

  //generator level particles
  Short_t ngen;
  Int_t gen_id[MAXGENPEREVENT];
  Float_t gen_pt[MAXGENPEREVENT], gen_eta[MAXGENPEREVENT], gen_phi[MAXGENPEREVENT], gen_en[MAXGENPEREVENT];

  //generator level jets
  Short_t njgen;
  Float_t genj_pt[MAXGENPEREVENT], genj_eta[MAXGENPEREVENT], genj_phi[MAXGENPEREVENT], genj_en[MAXGENPEREVENT];
  Float_t genj_emfrac[MAXGENPEREVENT], genj_hadfrac[MAXGENPEREVENT], genj_invfrac[MAXGENPEREVENT];
  
  //sim hits and ADC counts
  Int_t nhits;
  Short_t hit_layer[MAXHGCHITSPEREVENT];
  Float_t hit_x[MAXHGCHITSPEREVENT],hit_y[MAXHGCHITSPEREVENT],hit_z[MAXHGCHITSPEREVENT];
  Float_t hit_eta[MAXHGCHITSPEREVENT],hit_phi[MAXHGCHITSPEREVENT];
  Float_t hit_edep[MAXHGCHITSPEREVENT];
  
  HGCSimEvent_t()
  {
    ngen=0;
    njgen=0;
    nhits=0;
  }
};


void initHGCSimulationEventTree(TTree *t,HGCSimEvent_t &simEvt);


#endif
