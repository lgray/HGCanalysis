#ifndef _hgcsimulationevent_h_
#define _hgcsimulationevent_h_

#include "TTree.h"

#define MAXGENPEREVENT 10
#define MAXG4PEREVENT 1500
#define MAXTKSPEREVENT 1000
#define MAXHGCHITSPEREVENT 1000000
#define MAXLAYERSINGEO 100
#define MAXSDINGEO 5

struct  HGCSimEvent_t
{ 
  //event identifier
  Int_t run, lumi,event;

  //geometry information
  Short_t nsd,nlay[MAXSDINGEO];

  //generator level
  Short_t ngen;
  Int_t gen_id[MAXGENPEREVENT];
  Float_t gen_pt[MAXGENPEREVENT], gen_eta[MAXGENPEREVENT], gen_phi[MAXGENPEREVENT], gen_en[MAXGENPEREVENT];
  
  //geant4 information
  Short_t   ng4;
  Int_t g4_id[MAXG4PEREVENT];
  Float_t g4_vtx[MAXG4PEREVENT],g4_vty[MAXG4PEREVENT],g4_vtz[MAXG4PEREVENT];
  Float_t g4_en[MAXG4PEREVENT], g4_eta[MAXG4PEREVENT],g4_phi[MAXG4PEREVENT];
  Float_t g4_dEnInTracker,g4_dEnIonInTracker;
  
  //tracks and their extrapolation to HGCal
  Int_t ntk;
  Short_t tk_nhits[MAXTKSPEREVENT];
  Float_t tk_pt[MAXTKSPEREVENT],tk_eta[MAXTKSPEREVENT],tk_phi[MAXTKSPEREVENT],tk_chi2[MAXTKSPEREVENT];
  Float_t tk_extrapol_x[MAXTKSPEREVENT][66],tk_extrapol_y[MAXTKSPEREVENT][66];

  //sim hits and ADC counts
  Int_t nhits;
  Short_t hit_type[MAXHGCHITSPEREVENT], hit_layer[MAXHGCHITSPEREVENT];
  Float_t hit_x[MAXHGCHITSPEREVENT],hit_y[MAXHGCHITSPEREVENT],hit_z[MAXHGCHITSPEREVENT];
  Float_t hit_eta[MAXHGCHITSPEREVENT],hit_phi[MAXHGCHITSPEREVENT];
  Float_t hit_edep[MAXHGCHITSPEREVENT],hit_edep_sample[MAXHGCHITSPEREVENT][5];
  Short_t digi_adc[MAXHGCHITSPEREVENT];
  
  HGCSimEvent_t()
  {
    ng4=0;
    ngen=0;
    ntk=0;
    nsd=0;
    nhits=0;
  }
};


void initHGCSimulationEventTree(TTree *t,HGCSimEvent_t &simEvt);


#endif
