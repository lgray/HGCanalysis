#include "UserCode/HGCanalysis/plugins/HGCSimHitsAnalyzer.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DetectorDescription/OfflineDBLoader/interface/GeometryInfoDump.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/FCalGeometry/interface/HGCalGeometry.h"

#include "SimG4CMS/Calo/interface/CaloHitID.h"

#include "DetectorDescription/Core/interface/DDFilter.h"
#include "DetectorDescription/Core/interface/DDFilteredView.h"
#include "DetectorDescription/Core/interface/DDSolid.h"

#include "CLHEP/Units/GlobalSystemOfUnits.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <iostream>

using namespace std;

//
HGCSimHitsAnalyzer::HGCSimHitsAnalyzer( const edm::ParameterSet &iConfig )
{
  //configure analyzer
  hitCollections_   = iConfig.getUntrackedParameter< std::vector<std::string> >("hitCollections");
  geometrySource_   = iConfig.getUntrackedParameter< std::vector<std::string> >("geometrySource");
  mipEn_            = iConfig.getUntrackedParameter< std::vector<double> >("mipEn");
  thrList_          = iConfig.getUntrackedParameter< std::vector<double> >("thrList");

  //init tree
  edm::Service<TFileService> fs;
  t_=fs->make<TTree>("HGC","Event Summary");
  t_->Branch("genId",  &genId_,  "genId/I");
  t_->Branch("genEn",  &genEn_,  "genEn/F");  
  t_->Branch("genEta", &genEta_, "genEta/F");
  t_->Branch("genPhi", &genPhi_, "genPhi/F");
  t_->Branch("nlay",   &nlay_,   "nlay/I");
  for(size_t i=0; i<thrList_.size(); i++)
    {
      TString pf("_thr"); pf+=i;
      edeps_.push_back( new Float_t[100] ); 
      nhits_.push_back( new Int_t[100] );
      t_->Branch("edep"+pf,  edeps_[i], "edeps"+pf+"[nlay]/F");
      t_->Branch("nhits"+pf, nhits_[i], "nhits"+pf+"[nlay]/I");
    }
}

//
HGCSimHitsAnalyzer::~HGCSimHitsAnalyzer()
{
}

//
void HGCSimHitsAnalyzer::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  //generator level particles
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel("genParticles", genParticles);
  if(genParticles->size()>1) std::cout << "[Warning] found more than 1 gen particle, will save only first one" << std::endl;
  const reco::GenParticle & p = (*genParticles)[0];
  genId_  = p.pdgId();
  genEn_  = p.energy();
  genEta_ = p.eta();
  genPhi_ = p.phi();
    
  //hits + geometry
  std::map<int, std::vector<float> > eventMaxEdep;
  std::map<int, std::vector<std::vector<float> > > eventEdeps;
  std::map<int, std::vector<std::vector<int> > > eventNhits;
  for(size_t i=0; i<hitCollections_.size(); i++)
    {
      uint32_t mySubDet(ForwardSubdetector::HGCEE);
      if(i==1) mySubDet=ForwardSubdetector::HGCHEF;
      else if(i==2) mySubDet=ForwardSubdetector::HGCHEB;

      //get hits+geometry
      edm::Handle<edm::PCaloHitContainer> caloHits;
      iEvent.getByLabel(edm::InputTag("g4SimHits",hitCollections_[i]),caloHits); 
      edm::ESHandle<HGCalGeometry> geom;
      iSetup.get<IdealGeometryRecord>().get(geometrySource_[i],geom);
      const HGCalTopology &topo=geom->topology();
      const HGCalDDDConstants &dddConst=topo.dddConstants();

      //prepare counters
      std::vector< std::vector<float> > layerEdeps( thrList_.size() );
      std::vector< std::vector<int> > layerNhits( thrList_.size() );
      for(size_t ithr=0; ithr<thrList_.size(); ithr++)
	{
	  layerEdeps[ithr]=std::vector<float>(dddConst.layers(true),0.);
	  layerNhits[ithr]=std::vector<int>(dddConst.layers(true),0.);
	}
      
      for(edm::PCaloHitContainer::const_iterator hit_it = caloHits->begin(); hit_it != caloHits->end(); ++hit_it) 
	{
	  //gang SIM->RECO cells to get final layer assignment
	  HGCalDetId simId(hit_it->id());
	  int layer(simId.layer()),cell(simId.cell());
	  std::pair<int,int> recoLayerCell=dddConst.simToReco(cell,layer,topo.detectorType());
	  cell  = recoLayerCell.first;
	  layer = recoLayerCell.second;
	  if(layer<0) continue;

	  //get global position
	  uint32_t recoDetId( (i==0) ?
			      (uint32_t)HGCEEDetId(ForwardSubdetector(mySubDet),simId.zside(),layer,simId.sector(),simId.subsector(),cell) :
			      (uint32_t)HGCHEDetId(ForwardSubdetector(mySubDet),simId.zside(),layer,simId.sector(),simId.subsector(),cell)
			      );
	  const GlobalPoint pos( std::move( geom->getPosition(recoDetId) ) );	  
	  float hitEta(pos.eta());
	  if(hitEta*genEta_<0) continue;    
	  
	  //save it if interesting
	  float hitEnInMIPs(hit_it->energy()*TMath::TanH(hitEta)/mipEn_[i]);
	  if(hitEnInMIPs<0.5) continue;
	  for(size_t ithr=0; ithr<thrList_.size(); ithr++)
	    {
	      if(hitEnInMIPs<thrList_[ithr]) continue;
	      layerEdeps[ithr][layer-1]+=hitEnInMIPs;
	      layerNhits[ithr][layer-1]++;
	    }
	}
      eventEdeps[i]=layerEdeps;
      eventNhits[i]=layerNhits;
    }
  
  //fill the tree (reset first for all layers)
  for(size_t ithr=0; ithr<thrList_.size(); ithr++)
    for(size_t ilay=0; ilay<100; ilay++)
      {
	edeps_[ithr][ilay]=0;
	nhits_[ithr][ilay]=0;
      }
  nlay_=0;
  for(std::map<int, std::vector<std::vector<float> > >::iterator it = eventEdeps.begin(); it!=eventEdeps.end(); it++)
    {
      for(size_t ilay=0; ilay<(it->second)[0].size(); ++ilay)
	{
	  for(size_t ithr=0; ithr<thrList_.size(); ithr++)
	    {
	      edeps_[ithr][nlay_]          = (it->second)[ithr][ilay];
	      nhits_[ithr][nlay_]          = eventNhits[it->first][ithr][ilay];
	    }
	  nlay_++;
	}
    }
  
  //all done, save information
  t_->Fill();
}


//define this as a plug-in
DEFINE_FWK_MODULE(HGCSimHitsAnalyzer);
