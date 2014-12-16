#include "UserCode/HGCanalysis/interface/HGCSimulationEvent.h"

void initHGCSimulationEventTree(TTree *t,HGCSimEvent_t &simEvt)
{
  t->Branch("run",       &simEvt.run,        "run/I");
  t->Branch("lumi",      &simEvt.lumi,       "lumi/I");
  t->Branch("event",     &simEvt.event,      "event/I");

  t->Branch("ngen",      &simEvt.ngen,       "ngen/S");
  t->Branch("gen_id",     simEvt.gen_id,     "gen_id[ngen]/I");
  t->Branch("gen_pt",     simEvt.gen_pt,     "gen_pt[ngen]/F");
  t->Branch("gen_eta",    simEvt.gen_eta,    "gen_eta[ngen]/F");
  t->Branch("gen_phi",    simEvt.gen_phi,    "gen_phi[ngen]/F");
  t->Branch("gen_en",     simEvt.gen_en,     "gen_en[ngen]/F"); 

  t->Branch("njgen",       &simEvt.njgen,        "njgen/S");
  t->Branch("genj_pt",      simEvt.genj_pt,      "genj_pt[njgen]/F");
  t->Branch("genj_eta",     simEvt.genj_eta,     "genj_eta[njgen]/F");
  t->Branch("genj_phi",     simEvt.genj_phi,     "genj_phi[njgen]/F");
  t->Branch("genj_en",      simEvt.genj_en,      "genj_en[njgen]/F"); 
  t->Branch("genj_emfrac",  simEvt.genj_emfrac,  "genj_emfrac[njgen]/F"); 
  t->Branch("genj_hadfrac", simEvt.genj_hadfrac, "genj_hadfrac[njgen]/F"); 
  t->Branch("genj_invfrac", simEvt.genj_invfrac, "genj_invfrac[njgen]/F"); 

  t->Branch("nhits",     &simEvt.nhits,      "nhits/I");
  t->Branch("hit_layer",  simEvt.hit_layer,  "hit_layer[nhits]/S");
  t->Branch("hit_edep",   simEvt.hit_edep,   "hit_edep[nhits]/F");
  t->Branch("hit_x",      simEvt.hit_x,      "hit_x[nhits]/F");
  t->Branch("hit_y",      simEvt.hit_y,      "hit_y[nhits]/F");
  t->Branch("hit_z",      simEvt.hit_z,      "hit_z[nhits]/F");
  t->Branch("hit_eta",    simEvt.hit_eta,    "hit_eta[nhits]/F");
  t->Branch("hit_phi",    simEvt.hit_phi,    "hit_phi[nhits]/F");
}
