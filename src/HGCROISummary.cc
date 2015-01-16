#include "UserCode/HGCanalysis/interface/HGCROISummary.h"

void initHGCROITree(TTree *t,HGCROISummary_t &simEvt)
{
  t->Branch("run",       &simEvt.run,        "run/I");
  t->Branch("lumi",      &simEvt.lumi,       "lumi/I");
  t->Branch("event",     &simEvt.event,      "event/I");

  t->Branch("gen_id",    &simEvt.gen_id,     "gen_id/I");
  t->Branch("gen_pt",    &simEvt.gen_pt,     "gen_pt/F");
  t->Branch("gen_eta",   &simEvt.gen_eta,    "gen_eta/F");
  t->Branch("gen_phi",   &simEvt.gen_phi,    "gen_phi/F");
  t->Branch("gen_en",    &simEvt.gen_en,     "gen_en/F"); 
  t->Branch("gen_emfrac",    &simEvt.gen_emfrac,     "gen_emfrac/F"); 
  t->Branch("gen_hadfrac",    &simEvt.gen_hadfrac,     "gen_hadfrac/F"); 
  t->Branch("gen_invfrac",    &simEvt.gen_invfrac,     "gen_invfrac/F"); 

  t->Branch("en", simEvt.en, "en[4][3]/F");
  t->Branch("nhits", simEvt.nhits, "nhits[4][3]/I");
  t->Branch("eta", simEvt.eta, "eta[4][3]/F");
  t->Branch("phi", simEvt.phi, "phi[4][3]/F");
  t->Branch("width", simEvt.width, "width[4][3]/F");
  t->Branch("totalVolume", simEvt.totalVolume, "totalVolume[4][3]/F");
  t->Branch("x", simEvt.x, "x[4][54]/F");
  t->Branch("y", simEvt.y, "y[4][54]/F");

  t->Branch("wgt_en", simEvt.wgt_en, "wgt_en[4][3]/F");
  t->Branch("wgt_nhits", simEvt.wgt_nhits, "wgt_nhits[4][3]/I");
  t->Branch("wgt_eta", simEvt.wgt_eta, "wgt_eta[4][3]/F");
  t->Branch("wgt_phi", simEvt.wgt_phi, "wgt_phi[4][3]/F");
  t->Branch("wgt_width", simEvt.wgt_width, "wgt_width[4][3]/F");
  t->Branch("wgt_totalVolume", simEvt.wgt_totalVolume, "wgt_totalVolume[4][3]/F");
  t->Branch("wgt_x", simEvt.wgt_x, "wgt_x[4][54]/F");
  t->Branch("wgt_y", simEvt.wgt_y, "wgt_y[4][54]/F");

  t->Branch("z", simEvt.z, "z[54]/F");

}

//
void attachHGCSimulationEventTree(TTree *t,HGCROISummary_t &simEvt)
{
  t->SetBranchAddress("run",       &simEvt.run);
  t->SetBranchAddress("lumi",      &simEvt.lumi);
  t->SetBranchAddress("event",     &simEvt.event);

  t->SetBranchAddress("gen_id",    &simEvt.gen_id);
  t->SetBranchAddress("gen_pt",    &simEvt.gen_pt);
  t->SetBranchAddress("gen_eta",   &simEvt.gen_eta);
  t->SetBranchAddress("gen_phi",   &simEvt.gen_phi);
  t->SetBranchAddress("gen_en",    &simEvt.gen_en);
  t->SetBranchAddress("gen_emfrac",  &simEvt.gen_emfrac);
  t->SetBranchAddress("gen_hadfrac", &simEvt.gen_hadfrac);
  t->SetBranchAddress("gen_invfrac", &simEvt.gen_invfrac);

  t->SetBranchAddress("en", simEvt.en);
  t->SetBranchAddress("nhits", simEvt.nhits);
  t->SetBranchAddress("eta", simEvt.eta);
  t->SetBranchAddress("phi", simEvt.phi);
  t->SetBranchAddress("width", simEvt.width);
  t->SetBranchAddress("totalVolume", simEvt.totalVolume);
  t->SetBranchAddress("x", simEvt.x);
  t->SetBranchAddress("y", simEvt.y);

  t->SetBranchAddress("wgt_en", simEvt.wgt_en);
  t->SetBranchAddress("wgt_nhits", simEvt.wgt_nhits);
  t->SetBranchAddress("wgt_eta", simEvt.wgt_eta);
  t->SetBranchAddress("wgt_phi", simEvt.wgt_phi);
  t->SetBranchAddress("wgt_width", simEvt.wgt_width);
  t->SetBranchAddress("wgt_totalVolume", simEvt.wgt_totalVolume);
  t->SetBranchAddress("wgt_x", simEvt.wgt_x);
  t->SetBranchAddress("wgt_y", simEvt.wgt_y);

  t->SetBranchAddress("z", simEvt.z);

}
