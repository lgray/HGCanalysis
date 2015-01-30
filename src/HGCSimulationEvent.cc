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
  t->Branch("hit_wgt",    simEvt.hit_wgt,    "hit_wgt[nhits]/F");
  t->Branch("hit_wgt_t",  simEvt.hit_wgt_t,  "hit_wgt_t[nhits]/F");
  t->Branch("hit_x",      simEvt.hit_x,      "hit_x[nhits]/F");
  t->Branch("hit_y",      simEvt.hit_y,      "hit_y[nhits]/F");
  t->Branch("hit_z",      simEvt.hit_z,      "hit_z[nhits]/F");
  t->Branch("hit_eta",    simEvt.hit_eta,    "hit_eta[nhits]/F");
  t->Branch("hit_phi",    simEvt.hit_phi,    "hit_phi[nhits]/F");
}

//
void attachHGCSimulationEventTree(TTree *t,HGCSimEvent_t &simEvt)
{
  t->SetBranchAddress("run",       &simEvt.run);
  t->SetBranchAddress("lumi",      &simEvt.lumi);
  t->SetBranchAddress("event",     &simEvt.event);

  t->SetBranchAddress("ngen",      &simEvt.ngen);
  t->SetBranchAddress("gen_id",     simEvt.gen_id);
  t->SetBranchAddress("gen_pt",     simEvt.gen_pt);
  t->SetBranchAddress("gen_eta",    simEvt.gen_eta);
  t->SetBranchAddress("gen_phi",    simEvt.gen_phi);
  t->SetBranchAddress("gen_en",     simEvt.gen_en);

  t->SetBranchAddress("njgen",       &simEvt.njgen);
  t->SetBranchAddress("genj_pt",      simEvt.genj_pt);
  t->SetBranchAddress("genj_eta",     simEvt.genj_eta);
  t->SetBranchAddress("genj_phi",     simEvt.genj_phi);
  t->SetBranchAddress("genj_en",      simEvt.genj_en);
  t->SetBranchAddress("genj_emfrac",  simEvt.genj_emfrac);
  t->SetBranchAddress("genj_hadfrac", simEvt.genj_hadfrac);
  t->SetBranchAddress("genj_invfrac", simEvt.genj_invfrac);

  t->SetBranchAddress("nhits",     &simEvt.nhits);
  t->SetBranchAddress("hit_layer",  simEvt.hit_layer);
  t->SetBranchAddress("hit_edep",   simEvt.hit_edep);
  t->SetBranchAddress("hit_wgt",    simEvt.hit_wgt);
  t->SetBranchAddress("hit_wgt_t",  simEvt.hit_wgt_t);
  t->SetBranchAddress("hit_x",      simEvt.hit_x);
  t->SetBranchAddress("hit_y",      simEvt.hit_y);
  t->SetBranchAddress("hit_z",      simEvt.hit_z);
  t->SetBranchAddress("hit_eta",    simEvt.hit_eta);
  t->SetBranchAddress("hit_phi",    simEvt.hit_phi);
}
