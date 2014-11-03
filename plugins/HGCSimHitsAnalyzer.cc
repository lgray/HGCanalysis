#include "UserCode/HGCanalysis/plugins/HGCSimHitsAnalyzer.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFraction.h"

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
  hitCollections_        = iConfig.getUntrackedParameter< std::vector<std::string> >("hitCollections");
  recHitCollections_     = iConfig.getUntrackedParameter< std::vector<std::string> >("recHitCollections");
  pfClustersCollections_ = iConfig.getUntrackedParameter< std::vector<std::string> >("pfClustersCollections");
  geometrySource_        = iConfig.getUntrackedParameter< std::vector<std::string> >("geometrySource");
  mipEn_                 = iConfig.getUntrackedParameter< std::vector<double> >("mipEn");
  thrList_               = iConfig.getUntrackedParameter< std::vector<double> >("thrList");

  //init tree
  edm::Service<TFileService> fs;
  t_=fs->make<TTree>("HGC","Event Summary");
  t_->Branch("genId",  &genId_,  "genId/I");
  t_->Branch("genEn",  &genEn_,  "genEn/F");  
  t_->Branch("genEta", &genEta_, "genEta/F");
  t_->Branch("genPhi", &genPhi_, "genPhi/F");
  t_->Branch("nlay",   &nlay_,   "nlay/I");
  for(size_t istep=0; istep<3; istep++)
    {
      TString key("sim");
      if(istep==1) key="rec";
      if(istep==2) key="clus";
      std::vector<Float_t *> templVF;
      edeps_[key]    = templVF;
      emeanPhi_[key] = templVF;
      emeanEta_[key] = templVF;
      sihih_[key]    = templVF;
      sipip_[key]    = templVF;
      sipih_[key]    = templVF;
      std::vector<Int_t *> templVI;
      nhits_[key]    = templVI;
      for(size_t i=0; i<thrList_.size(); i++)
	{
	  TString pf("_thr"); pf+=i;
	  edeps_[key].push_back( new Float_t[100] ); 
	  emeanPhi_[key].push_back( new Float_t[100] ); 
	  emeanEta_[key].push_back( new Float_t[100] ); 
	  sihih_[key].push_back( new Float_t[100] ); 
	  sipip_[key].push_back( new Float_t[100] ); 
	  sipih_[key].push_back( new Float_t[100] ); 
	  nhits_[key].push_back( new Int_t[100] );
	  t_->Branch("edep_"+key+pf,      edeps_[key][i],    "edeps_"+key+pf+"[nlay]/F");
	  t_->Branch("emeanPhi_"+key+pf,  emeanPhi_[key][i], "emeanPhi_"+key+pf+"[nlay]/F");
	  t_->Branch("emeanEta_"+key+pf,  emeanEta_[key][i], "emeanEta_"+key+pf+"[nlay]/F");
	  t_->Branch("sihih_"+key+pf,     sihih_[key][i],    "sihih_"+key+pf+"[nlay]/F");
	  t_->Branch("sipip_"+key+pf,     sipip_[key][i],    "sipip_"+key+pf+"[nlay]/F");
	  t_->Branch("sipih_"+key+pf,     sipih_[key][i],    "sipih_"+key+pf+"[nlay]/F");
	  t_->Branch("nhits_"+key+pf,     nhits_[key][i],    "nhits_"+key+pf+"[nlay]/I");
	}
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
  size_t maxGenParts(genParticles->size());
  if(maxGenParts>2) maxGenParts=2;
  if(genParticles->size()>maxGenParts) std::cout << "[Warning] found more than " << maxGenParts << " gen particles, will save only first " << maxGenParts << std::endl;
  
  //loop over generator level particles
  for(size_t igen=0; igen<maxGenParts; igen++)
    {
      //start new particle
      resetCounters();

      const reco::GenParticle & p = (*genParticles)[igen];
      genId_  = p.pdgId();
      genEn_  = p.energy();
      genEta_ = p.eta();
      genPhi_ = p.phi();

      std::cout << genId_ << " " << genEn_ << " " << genEta_ << " " << genPhi_ << std::endl;

      //hits and clusters
      std::map<uint32_t,float> templEdep;
      std::map<TString, std::map<uint32_t,float> > allEdeps;
      std::map<TString, std::pair<uint32_t,float> > maxEdep;
      for(size_t i=0; i<hitCollections_.size(); i++)
	{
	  uint32_t mySubDet(ForwardSubdetector::HGCEE);
	  if(i==1) mySubDet=ForwardSubdetector::HGCHEF;
	  else if(i==2) mySubDet=ForwardSubdetector::HGCHEB;
	  
	  //get geometry
	  edm::ESHandle<HGCalGeometry> geom;
	  iSetup.get<IdealGeometryRecord>().get(geometrySource_[i],geom);
	  const HGCalTopology &topo=geom->topology();
	  const HGCalDDDConstants &dddConst=topo.dddConstants();

	  //SIM HITS
	  TString key("sim");
	  allEdeps[key] = templEdep;
	  maxEdep[key]=std::pair<uint32_t,float>(0,-1.0);
	  edm::Handle<edm::PCaloHitContainer> caloHits;
	  iEvent.getByLabel(edm::InputTag("g4SimHits",hitCollections_[i]),caloHits); 
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
	      uint32_t recoDetId = ( (i==0) ?
				     (uint32_t)HGCEEDetId(ForwardSubdetector(mySubDet),simId.zside(),layer,simId.sector(),simId.subsector(),cell) :
				     (uint32_t)HGCHEDetId(ForwardSubdetector(mySubDet),simId.zside(),layer,simId.sector(),simId.subsector(),cell)
				     );
	      
	      const GlobalPoint pos( std::move( geom->getPosition(recoDetId) ) );	  
	      float hitEta(pos.eta());
	      if(hitEta*genEta_<0) continue;    
	      
	      //save it if interesting (convert energy to keV)
	      float hitEnInMIPs(hit_it->energy()*1e6/mipEn_[i]);
	      if(hitEnInMIPs<0.5) continue;
	      if(allEdeps[key].find(recoDetId)==allEdeps[key].end()) allEdeps[key][recoDetId] = 0;
	      allEdeps[key][recoDetId] += hitEnInMIPs;
	      
	      //check if maximum found
	      if(maxEdep[key].second>allEdeps[key][recoDetId]) continue;
	      maxEdep[key].first=recoDetId;
	      maxEdep[key].second=allEdeps[key][recoDetId];
	    }
	  cout << "\t" << key << " " << maxEdep[key].first << " " << maxEdep[key].second << endl;

	  //RECO: save rec hits where sim hits exist
	  key="rec";
	  allEdeps[key] = templEdep;
	  maxEdep[key]  = std::pair<uint32_t,float>(0,-1.0);
	  edm::Handle<HGCRecHitCollection> recHits;
	  iEvent.getByLabel(edm::InputTag("HGCalRecHit",recHitCollections_[i]),recHits);
	  for(HGCRecHitCollection::const_iterator hit_it=recHits->begin(); hit_it!=recHits->end(); hit_it++)
	    {
	      uint32_t recoDetId(hit_it->id());
	      if(allEdeps["sim"].find(recoDetId)==allEdeps["sim"].end()) continue;

	      //convert energy to keV
	      float hitEn(hit_it->energy()*1e6/mipEn_[i]);
	      if(allEdeps[key].find(recoDetId)==allEdeps[key].end()) allEdeps[key][recoDetId] = 0;
	      allEdeps[key][recoDetId] += hitEn;
	      
	      //check if maximum found
	      if(maxEdep[key].second>allEdeps[key][recoDetId]) continue;
	      maxEdep[key].first=recoDetId;
	      maxEdep[key].second=allEdeps[key][recoDetId];
	    } 
	  cout << "\t" << key << " " << maxEdep[key].first << " " << maxEdep[key].second << endl;
	

	  //CLUSTERS
	  key="clus";
	  allEdeps[key] = templEdep;
	  maxEdep[key]=std::pair<uint32_t,float>(0,-1.0);
	  edm::Handle<reco::PFClusterCollection> pfClusters;
	  iEvent.getByLabel(edm::InputTag(pfClustersCollections_[i],""),pfClusters);
	  
	  //get all within DR=0.4
	  int nClusters(0);
	  for(reco::PFClusterCollection::const_iterator c_it=pfClusters->begin(); 
	      c_it!=pfClusters->end(); 
	      c_it++)
	    {
	      float dR=deltaR(c_it->position(),p);
	      if(dR>0.4) continue;
	      nClusters++;
	      for( const auto& rhf : c_it->recHitFractions() ) {
		const reco::PFRecHit& hit = *(rhf.recHitRef());
		uint32_t recoDetId(hit.detId());
		
		//rec hits : convert energy to keV
		float eclus(hit.energy()*1e6/mipEn_[i]);
		eclus *= rhf.fraction();
		if(allEdeps[key].find(recoDetId)==allEdeps[key].end())  allEdeps[key][recoDetId] = 0;
		allEdeps[key][recoDetId] += eclus;
		
		//check if maximum found
		if(maxEdep[key].second>allEdeps[key][recoDetId]) continue;
		maxEdep[key].first=recoDetId;
		maxEdep[key].second=allEdeps[key][recoDetId];
	      }
	    }
	  cout << "\t" << key << " " << maxEdep[key].first << " " << maxEdep[key].second << endl;
	  
	  
	  //now compute the relevant variables for the regression
	  /*
	  for(std::map<TString, std::map<uint32_t,float> >::iterator stepIt=allEdeps.begin();
	      stepIt!=allEdeps.end();
	      stepIt++)
	    {
	      TString key(stepIt->first);

	      const GlobalPoint refPos( std::move( geom->getPosition(maxEdep[key].first) ) );
	      float refEta=refPos.eta();
	      float refPhi=refPos.phi();

	      for(std::map<uint32_t,float>::iterator detIt=stepIt->second.begin();
		  detIt!=stepIt->second.end();
		  detIt++)
		{
		  const GlobalPoint pos( std::move( geom->getPosition(detIt->first) ) );
		  float hitPhi=pos.phi();
		  float hitEta=pos.eta();
		  int layer((detIt->first >> 19) & 0x1f);
		  float en(detIt->second);
		  for(size_t ithr=0; ithr<thrList_.size(); ithr++)
		    {
		      if(en<thrList_[ithr]) continue;
		      nhits_[key][ithr][layer]    ++;
		      edeps_[key][ithr][layer]    += en;
		      emeanPhi_[key][ithr][layer] += en*hitPhi;
		      emeanEta_[key][ithr][layer] += en*hitEta;
		      sihih_[key][ithr][layer]    += en*pow(hitEta-refEta,2);
		      sipip_[key][ithr][layer]    += en*pow(hitPhi-refPhi,2);
		      sipih_[key][ithr][layer]    += en*(hitEta-refEta)*(hitPhi-refPhi);
		    }
		}
	    }
	  */
	}
    }

  /*
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
  */
}


//define this as a plug-in
DEFINE_FWK_MODULE(HGCSimHitsAnalyzer);
