#ifndef HGCSimpleHitAnalyzer_h
#define HGCSimpleHitAnalyzer_h

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h" 
#include "DataFormats/Candidate/interface/Candidate.h"

#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"

#include "DetectorDescription/OfflineDBLoader/interface/GeometryInfoDump.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/FCalGeometry/interface/HGCalGeometry.h" 
#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"
#include "Geometry/HcalCommonData/interface/HcalDDDRecConstants.h"
#include "Geometry/Records/interface/HcalRecNumberingRecord.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1F.h"
#include "TH2F.h"

#include <iostream>
#include <vector>
#include <map>
#include <string>

struct GangedHitInfo_t
{
  int layer,subdet;
  float energy,time;
  float eta,etagen,phi,phigen;
  float q;
};

class HGCSimpleHitAnalyzer : public edm::EDAnalyzer 
{
 public:
  
  explicit HGCSimpleHitAnalyzer(const edm::ParameterSet&);
  ~HGCSimpleHitAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  void analyzeHits (std::vector<PCaloHit>& hits);
  
private:
  void countHit(int layer,float en, float time, float eta,float geneta, float phi,float genphi,float q);
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  bool defineGeometry(edm::ESTransientHandle<DDCompactView> &ddViewH);
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  
  // ----------member data ---------------------------
  std::vector<std::string> hitCollections_, geometrySource_;       
  std::map<TString, std::map<int,TH1F *> > histos_;
  std::map<TString, std::map<int,TH2F *> > histos2D_;

  TH2F *sdH_;
  const HcalDDDRecConstants *hcalDDD_;
};
#endif
