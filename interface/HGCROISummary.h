#ifndef _hgcroisummary_h_
#define _hgcroisummary_h_

#include "UserCode/HGCanalysis/interface/ROIInfo.h"

#include "TTree.h"

#define MAXVTX  500
#define MAXNTKJ 500
#define MAXGENJ 500
#define MAXROI  500

class  HGCROISummary
{ 
 public:

  //event identifier
  Int_t run,lumi,event;

  //gen jet information
  Int_t ngen, gen_id[MAXGENJ],gen_status[MAXGENJ];
  Float_t gen_pt[MAXGENJ], gen_eta[MAXGENJ], gen_phi[MAXGENJ], gen_emfrac[MAXGENJ], gen_hadfrac[MAXGENJ], gen_invfrac[MAXGENJ];
  Float_t gen_mjj, gen_detajj;

  //vertex information
  Float_t gen_x,gen_y,gen_z;
  Int_t nvtx, vtx_rank[MAXVTX],vtx_ntk[MAXVTX];
  Float_t vtx_x[MAXVTX], vtx_y[MAXVTX], vtx_z[MAXVTX], vtx_pt[MAXVTX],vtx_normchi2[MAXVTX]; 

  //track jet information
  Int_t ntkj, tkj_ntk[MAXNTKJ], tkj_vtxIdx[MAXNTKJ], tkj_genIdx[MAXNTKJ];
  Float_t tkj_pt[MAXNTKJ], tkj_eta[MAXNTKJ], tkj_phi[MAXNTKJ];

  //regions of interest
  std::vector<ROIInfo> *roi_info;

  HGCROISummary() :  roi_info(new std::vector<ROIInfo>) { }
  inline void newEvent(int newRun, int newLumi,int newEvent) 
  { 
    run=newRun; lumi=newLumi; event=newEvent;
    nvtx=0; ntkj=0; ngen=0; 
    gen_x=0; gen_y=0; gen_z=0;
    gen_mjj=0; gen_detajj=0;
    roi_info->clear();
  }
  inline void addGenJet(int id,float pt, float eta, float phi, float emfrac, float hadfrac, float invfrac)
  {
    if(ngen>=MAXGENJ) return;
    gen_id[ngen]=id;
    gen_pt[ngen]=pt;         gen_eta[ngen]=eta;         gen_phi[ngen]=phi;
    gen_emfrac[ngen]=emfrac; gen_hadfrac[ngen]=hadfrac; gen_invfrac[ngen]=invfrac;
    ngen++;
  }
  inline void addGenVertex(float x, float y, float z)
  {
    gen_x=x; gen_y=y; gen_z=z;
  }
  inline void addGenJet(int id,int status,float pt, float eta, float phi, float emfrac, float hadfrac, float invfrac)
  {
    if(ngen>=MAXGENJ) return;
    gen_id[ngen]=id;
    gen_status[ngen]=status;
    gen_pt[ngen]=pt;
    gen_eta[ngen]=eta;
    gen_phi[ngen]=phi;
    gen_emfrac[ngen]=emfrac;
    gen_hadfrac[ngen]=hadfrac;
    gen_invfrac[ngen]=invfrac;
    ngen++;
  }  
  inline void addRecoVertex(int rank,int ntk, float x, float y, float z, float pt, float normchi2)
  {
    if(nvtx>MAXVTX) return;
    vtx_ntk[nvtx]=ntk;
    vtx_rank[nvtx]=rank;
    vtx_x[nvtx]=x;
    vtx_y[nvtx]=y;
    vtx_z[nvtx]=z;
    vtx_pt[nvtx]=pt;
    vtx_normchi2[nvtx]=normchi2;
    nvtx++;
  }
  inline void addTrackJet(int pvIdx,int genIdx, float pt, float eta, float phi, int ntk)
  {
    if(ntkj>=MAXNTKJ) return;
    tkj_ntk[ntkj]=ntk;
    tkj_vtxIdx[ntkj]=pvIdx;
    tkj_genIdx[ntkj]=genIdx;
    tkj_pt[ntkj]=pt;
    tkj_eta[ntkj]=eta;
    tkj_phi[ntkj]=phi;
    ntkj++;
  }

  ~HGCROISummary() { }
};


void initHGCROITree(TTree *t,HGCROISummary &simEvt);
void attachHGCROITree(TTree *t,HGCROISummary &simEvt);

#endif
