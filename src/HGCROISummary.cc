#include "UserCode/HGCanalysis/interface/HGCROISummary.h"

void initHGCROITree(TTree *t,HGCROISummary &simEvt)
{
  t->Branch("run",       &simEvt.run,        "run/I");
  t->Branch("lumi",      &simEvt.lumi,       "lumi/I");
  t->Branch("event",     &simEvt.event,      "event/I");

  t->Branch("gen_x",         &simEvt.gen_x,          "gen_x/F");
  t->Branch("gen_y",         &simEvt.gen_y,          "gen_y/F");
  t->Branch("gen_z",         &simEvt.gen_z,          "gen_z/F");

  t->Branch("nvtx",          &simEvt.nvtx,           "nvtx/I");
  t->Branch("vtx_ntk",        simEvt.vtx_ntk,        "vtx_ntk[nvtx]/I");
  t->Branch("vtx_rank",       simEvt.vtx_rank,       "vtx_rank[nvtx]/I");
  t->Branch("vtx_x",          simEvt.vtx_x,          "vtx_x[nvtx]/F");
  t->Branch("vtx_y",          simEvt.vtx_y,          "vtx_y[nvtx]/F");
  t->Branch("vtx_z",          simEvt.vtx_z,          "vtx_z[nvtx]/F");
  t->Branch("vtx_pt",         simEvt.vtx_pt,         "vtx_pt[nvtx]/F");
  t->Branch("vtx_normchi2",   simEvt.vtx_normchi2,   "vtx_normchi2[nvtx]/F");

  t->Branch("ntkj",          &simEvt.ntkj,           "ntkj/I");
  t->Branch("tkj_ntk",        simEvt.tkj_ntk,        "tkj_ntk[ntkj]/I");
  t->Branch("tkj_vtxIdx",     simEvt.tkj_vtxIdx,     "tkj_vtxIdx[ntkj]/I");
  t->Branch("tkj_genIdx",     simEvt.tkj_genIdx,     "tkj_genIdx[ntkj]/I");
  t->Branch("tkj_pt",         simEvt.tkj_pt,         "tkj_pt[ntkj]/F");
  t->Branch("tkj_eta",        simEvt.tkj_eta,        "tkj_eta[ntkj]/F");
  t->Branch("tkj_phi",        simEvt.tkj_phi,        "tkj_phi[ntkj]/F");

  t->Branch("ngen",           &simEvt.ngen,         "ngen/I");
  t->Branch("gen_id",          simEvt.gen_id,       "gen_id[ngen]/I");
  t->Branch("gen_status",      simEvt.gen_status,   "gen_status[ngen]/I");
  t->Branch("gen_pt",          simEvt.gen_pt,       "gen_pt[ngen]/F");
  t->Branch("gen_eta",         simEvt.gen_eta,      "gen_eta[ngen]/F");
  t->Branch("gen_phi",         simEvt.gen_phi,      "gen_phi[ngen]/F");
  t->Branch("gen_emfrac",      simEvt.gen_emfrac,   "gen_emfrac[ngen]/F"); 
  t->Branch("gen_hadfrac",     simEvt.gen_hadfrac,  "gen_hadfrac[ngen]/F"); 
  t->Branch("gen_invfrac",     simEvt.gen_invfrac,  "gen_invfrac[ngen]/F"); 
  t->Branch("gen_mjj",        &simEvt.gen_mjj,      "gen_mjj/F");
  t->Branch("gen_detajj",     &simEvt.gen_detajj,   "gen_detajj/F");

  t->Branch("roi_info",          "std::vector<ROIInfo>",      &simEvt.roi_info);
}

//
void attachHGCROITree(TTree *t,HGCROISummary &simEvt)
{
  t->SetBranchAddress("run",       &simEvt.run);
  t->SetBranchAddress("lumi",      &simEvt.lumi);
  t->SetBranchAddress("event",     &simEvt.event);

  t->SetBranchAddress("gen_x",         &simEvt.gen_x);
  t->SetBranchAddress("gen_y",         &simEvt.gen_y);
  t->SetBranchAddress("gen_z",         &simEvt.gen_z);

  t->SetBranchAddress("nvtx",          &simEvt.nvtx);
  t->SetBranchAddress("vtx_ntk",        simEvt.vtx_ntk);
  t->SetBranchAddress("vtx_rank",       simEvt.vtx_rank);
  t->SetBranchAddress("vtx_x",          simEvt.vtx_x);
  t->SetBranchAddress("vtx_y",          simEvt.vtx_y);
  t->SetBranchAddress("vtx_z",          simEvt.vtx_z);
  t->SetBranchAddress("vtx_pt",         simEvt.vtx_pt);
  t->SetBranchAddress("vtx_normchi2",   simEvt.vtx_normchi2);

  t->SetBranchAddress("ntkj",          &simEvt.ntkj);
  t->SetBranchAddress("tkj_ntk",        simEvt.tkj_ntk);
  t->SetBranchAddress("tkj_vtxIdx",     simEvt.tkj_vtxIdx);
  t->SetBranchAddress("tkj_genIdx",     simEvt.tkj_genIdx);
  t->SetBranchAddress("tkj_pt",         simEvt.tkj_pt);
  t->SetBranchAddress("tkj_eta",        simEvt.tkj_eta);
  t->SetBranchAddress("tkj_phi",        simEvt.tkj_phi);

  t->SetBranchAddress("ngen",           &simEvt.ngen);
  t->SetBranchAddress("gen_id",          simEvt.gen_id);
  t->SetBranchAddress("gen_status",      simEvt.gen_status);
  t->SetBranchAddress("gen_pt",          simEvt.gen_pt);
  t->SetBranchAddress("gen_eta",         simEvt.gen_eta);
  t->SetBranchAddress("gen_phi",         simEvt.gen_phi);
  t->SetBranchAddress("gen_emfrac",      simEvt.gen_emfrac);
  t->SetBranchAddress("gen_hadfrac",     simEvt.gen_hadfrac);
  t->SetBranchAddress("gen_invfrac",     simEvt.gen_invfrac);
  t->SetBranchAddress("gen_mjj",        &simEvt.gen_mjj);
  t->SetBranchAddress("gen_detajj",     &simEvt.gen_detajj);

  t->SetBranchAddress("roi_info",          &simEvt.roi_info);

}
