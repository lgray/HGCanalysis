#include "UserCode/HGCanalysis/plugins/HGCHitsAnalyzer.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "SimG4CMS/Calo/interface/CaloHitID.h"

#include "DetectorDescription/OfflineDBLoader/interface/GeometryInfoDump.h"
#include "DetectorDescription/Core/interface/DDFilter.h"
#include "DetectorDescription/Core/interface/DDFilteredView.h"
#include "DetectorDescription/Core/interface/DDSolid.h"

#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

#include "CLHEP/Units/GlobalSystemOfUnits.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TVector2.h"

#include <iostream>

using namespace std;

//
HGCHitsAnalyzer::HGCHitsAnalyzer( const edm::ParameterSet &iConfig )
{
  //configure analyzer
  genSource_        = iConfig.getUntrackedParameter<std::string>("genSource");
  genJetsSource_    = iConfig.getUntrackedParameter<std::string>("genJetsSource");
  geometrySource_   = iConfig.getUntrackedParameter< std::vector<std::string> >("geometrySource");
  hitCollections_   = iConfig.getUntrackedParameter< std::vector<std::string> >("hitCollections");
  mipEn_            = iConfig.getUntrackedParameter< std::vector<double> >("mipEn");

  //init tree
  edm::Service<TFileService> fs;
  t_=fs->make<TTree>("HGC","Event Summary");
  initHGCSimulationEventTree(t_,simEvt_);
}

//
HGCHitsAnalyzer::~HGCHitsAnalyzer()
{
}

//
void HGCHitsAnalyzer::analyzeGenParticles(edm::Handle<edm::View<reco::Candidate> > &genParticles, edm::Handle<reco::GenJetCollection> &genJets)
{
  //store a summary of stable gen particles
  simEvt_.ngen=0;
  if(genParticles.isValid())
    {
      for(size_t i = 0; i < genParticles->size(); ++ i)
	{
	  const reco::GenParticle & p = dynamic_cast<const reco::GenParticle &>( (*genParticles)[i] );
	  
	  //hard process
	  bool isOutgoingNeutrino( p.status()==1 && (abs(p.pdgId())==12||abs(p.pdgId())==14||abs(p.pdgId())==16) );
	  if(!isOutgoingNeutrino && p.status()!=3) continue;
	  if(abs(p.pdgId())==2212) continue;
	  
	  simEvt_.gen_id[simEvt_.ngen]=p.pdgId();
	  simEvt_.gen_pt[simEvt_.ngen]=p.pt();
	  simEvt_.gen_eta[simEvt_.ngen]=p.eta();
	  simEvt_.gen_phi[simEvt_.ngen]=p.phi();
	  simEvt_.gen_en[simEvt_.ngen]=p.energy();
	  simEvt_.ngen++;

	  if(simEvt_.ngen>=MAXGENPEREVENT) break;
	}
    }
     

  //store a summary of gen jets
  simEvt_.njgen=0;
  if(genJets.isValid())
    {
      for(size_t k=0; k<genJets->size(); ++k)
	{
	  const reco::GenJet & j =(*genJets)[k];
	  if(j.pt()<20 || fabs(j.eta())>4.7) continue;
	  simEvt_.genj_en[simEvt_.njgen]=j.energy();
	  simEvt_.genj_eta[simEvt_.njgen]=j.eta();
	  simEvt_.genj_phi[simEvt_.njgen]=j.phi();
	  simEvt_.genj_pt[simEvt_.njgen]=j.pt();
	  simEvt_.genj_emfrac[simEvt_.njgen]=j.emEnergy()/j.energy();
	  simEvt_.genj_hadfrac[simEvt_.njgen]=j.hadEnergy()/j.energy();
	  simEvt_.genj_invfrac[simEvt_.njgen]=j.invisibleEnergy()/j.energy();
	  simEvt_.njgen++;
	  if(simEvt_.njgen>=MAXGENPEREVENT) break;
	}
    }
  cout << simEvt_.ngen << " " << simEvt_.njgen << endl;
}
  
//
void HGCHitsAnalyzer::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  //event header
  simEvt_.run    = iEvent.id().run();
  simEvt_.lumi   = iEvent.luminosityBlock();
  simEvt_.event  = t_->GetEntriesFast()+1;
  
  //generator level particles
  edm::Handle<edm::View<reco::Candidate> > genParticles;
  iEvent.getByLabel(edm::InputTag(genSource_), genParticles);
  edm::Handle<reco::GenJetCollection> genJets;
  iEvent.getByLabel(edm::InputTag(genJetsSource_), genJets);
  analyzeGenParticles(genParticles,genJets);
    
  //read rec hits and geometry
  int layerCtrOffset(1);
  simEvt_.nhits=0;
  for(size_t i=0; i<geometrySource_.size(); i++)
    {
      edm::Handle<HGCRecHitCollection> recHits;
      iEvent.getByLabel(edm::InputTag("HGCalRecHit",hitCollections_[i]),recHits);
      
      edm::ESHandle<HGCalGeometry> geomH;
      iSetup.get<IdealGeometryRecord>().get(geometrySource_[i],geomH);
      const HGCalGeometry *geom=geomH.product();
      
      for(HGCRecHitCollection::const_iterator hit_it=recHits->begin(); hit_it!=recHits->end(); hit_it++)
	{
	  uint32_t recoDetId(hit_it->id());

	  //convert energy to keV
	  float hitEn(hit_it->energy()*1e6/mipEn_[i]);
	  if(hitEn<0.5) continue;
	  
	  //decode position
	  const GlobalPoint refPos( std::move( geom->getPosition(recoDetId) ) );
	  int layer( ((recoDetId >> 19) & 0x1f) + layerCtrOffset-1 );

	  simEvt_.hit_layer[simEvt_.nhits]=layer;
	  simEvt_.hit_x[simEvt_.nhits] = refPos.x();
	  simEvt_.hit_y[simEvt_.nhits] = refPos.y();
	  simEvt_.hit_z[simEvt_.nhits] = refPos.z();
	  simEvt_.hit_eta[simEvt_.nhits] = refPos.eta();
	  simEvt_.hit_phi[simEvt_.nhits] = refPos.phi();
	  simEvt_.hit_edep[simEvt_.nhits] = hitEn;
	  simEvt_.nhits++;
	}
      
      const HGCalTopology &topo=geom->topology();
      const HGCalDDDConstants &dddConst=topo.dddConstants();
      layerCtrOffset +=dddConst.layers(true);
    }

  
  
  //fill tree
  t_->Fill();
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCHitsAnalyzer);
