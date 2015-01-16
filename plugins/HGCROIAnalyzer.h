#ifndef _HGCROIAnalyzer_h_
#define _HGCROIAnalyzer_h_

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "SimG4CMS/Calo/interface/HGCNumberingScheme.h"
#include "Geometry/FCalGeometry/interface/HGCalGeometry.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "UserCode/HGCanalysis/interface/ROIManager.h"
#include "UserCode/HGCanalysis/interface/HGCSimulationEvent.h"
#include "UserCode/HGCanalysis/interface/HGCROISummary.h"

#include "TH2F.h"
#include "TH1F.h"
#include "TTree.h"

#include <string>

/**
   @class HGCROIAnalyzer
   @author P. Silva (CERN)
*/

class HGCROIAnalyzer : public edm::EDAnalyzer 
{
  
 public:
  
  explicit HGCROIAnalyzer( const edm::ParameterSet& );
  ~HGCROIAnalyzer();
  virtual void analyze( const edm::Event&, const edm::EventSetup& );

 private:

  virtual void endJob() ;

  int evtCtr_;

  //tree for hit summary
  TTree *hitT_;
  HGCSimEvent_t simEvt_;
  bool saveHitTree_;

  //ROI stuff
  ROIManager roi_;
  HGCROISummary_t roiEvt_;

  TTree *roiT_;
  TH2F *regsH_;
  Int_t nLayerBins_, nEtaBins_;

  TH2F *csidrH_,*csitdrH_;
  Int_t   ndRbins_,       nCsiBins_;
  Float_t drMin_, drMax_, csiMin_,csiMax_;

  TH2F *medianPU_csiH_,  *widthPU_csiH_;
  TH2F *medianPU_csitH_, *widthPU_csitH_;

  //
  bool taggingMode_;
  
  //
  edm::FileInPath roipuParamFile_;

  //gen level
  std::string genSource_, genJetsSource_;
  
  //hgcal
  std::vector<std::string> geometrySource_;
  std::vector<std::string> hitCollections_;
  std::vector<double> mipEn_;
};
 

#endif
