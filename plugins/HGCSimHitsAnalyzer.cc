#include "UserCode/HGCanalysis/plugins/HGCSimHitsAnalyzer.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFraction.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"

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
#include <unordered_map>

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

  //init tree
  edm::Service<TFileService> fs;
  t_=fs->make<TTree>("HGC","Event Summary");
  t_->Branch("genId",  &genId_,  "genId/I");
  t_->Branch("genEn",  &genEn_,  "genEn/F");  
  t_->Branch("genEta", &genEta_, "genEta/F");
  t_->Branch("genPhi", &genPhi_, "genPhi/F");
  t_->Branch("nlay",   &nlay_,   "nlay/I");
  for(size_t istep=0; istep<4; istep++)
    {
      TString key("sim");
      if(istep==1) key="rec";
      if(istep==2) key="clus";
      if(istep==3) key="pf";

      nhits_[key]    = new Int_t[100];
      t_->Branch("nhits_"+key,     nhits_[key],    "nhits_"+key+"[nlay]/I");

      edeps_[key]    = new Float_t[100];
      t_->Branch("edep_"+key,      edeps_[key],    "edeps_"+key+"[nlay]/F");
      
      edeps3x3_[key]    = new Float_t[100];
      t_->Branch("edep3x3_"+key,      edeps3x3_[key],    "edeps3x3_"+key+"[nlay]/F");

      edeps5x5_[key]    = new Float_t[100];
      t_->Branch("edep5x5_"+key,      edeps5x5_[key],    "edeps5x5_"+key+"[nlay]/F");

      emeanPhi_[key] = new Float_t[100];
      t_->Branch("emeanPhi_"+key,  emeanPhi_[key], "emeanPhi_"+key+"[nlay]/F");
      
      emeanEta_[key] = new Float_t[100];
      t_->Branch("emeanEta_"+key,  emeanEta_[key], "emeanEta_"+key+"[nlay]/F");
      
      sihih_[key]    = new Float_t[100];
      t_->Branch("sihih_"+key,     sihih_[key],    "sihih_"+key+"[nlay]/F");
      
      sipip_[key]    = new Float_t[100];
      t_->Branch("sipip_"+key,     sipip_[key],    "sipip_"+key+"[nlay]/F");

      sipih_[key]    = new Float_t[100];
      t_->Branch("sipih_"+key,     sipih_[key],    "sipih_"+key+"[nlay]/F");

      if(istep>=2)
	{	
	  nClusters_[key] = new Int_t[3];
	  t_->Branch("nClusters_"+key,          nClusters_[key],        "nClusters_"+key+"[3]/I");
	  
	  nHitsInClusters_[key] = new Int_t[3];
	  t_->Branch("nHitsInClusters_"+key,    nHitsInClusters_[key],  "nHitsInClusters_"+key+"[3]/I");
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

      //particle flow candidates
      edm::Handle<reco::PFCandidateCollection> pflow;
      iEvent.getByLabel("particleFlow",pflow);
      std::vector<int> selPFs;
      for(size_t ipf=0; ipf<pflow->size(); ipf++)
	{
	  const reco::PFCandidate &cand=(*pflow)[ipf];
	  float dr=deltaR(p,cand);
	  if(dr>0.1) continue;
	  selPFs.push_back(ipf);
	}

      //hits and clusters
      std::map<uint32_t,float> templEdep;
      std::map<TString, std::map<uint32_t,float> > allEdeps;
      std::map<TString, std::pair<uint32_t,float> > maxEdep;
      uint32_t baseLayerIdx(0);
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

	  //RECO: save rec hits where sim hits exist
	  key="rec";
	  allEdeps[key] = templEdep;
	  maxEdep[key]  = std::pair<uint32_t,float>(0,-1.0);
	  edm::Handle<HGCRecHitCollection> recHits;
	  iEvent.getByLabel(edm::InputTag("HGCalRecHit",recHitCollections_[i]),recHits);
	  std::unordered_map<uint32_t,uint32_t> recoDetIdMap;
	  uint32_t recHitCtr(0);
	  for(HGCRecHitCollection::const_iterator hit_it=recHits->begin(); hit_it!=recHits->end(); hit_it++,recHitCtr++)
	    {
	      uint32_t recoDetId(hit_it->id());
	      recoDetIdMap[recoDetId]=recHitCtr;
	      if(allEdeps["sim"].find(recoDetId)==allEdeps["sim"].end()) continue;

	      //convert energy to keV
	      float hitEn(hit_it->energy()*1e6/mipEn_[i]);
	      if(hitEn<0.5) continue;
	      if(allEdeps[key].find(recoDetId)==allEdeps[key].end()) allEdeps[key][recoDetId] = 0;
	      allEdeps[key][recoDetId] += hitEn;
	      
	      //check if maximum found
	      if(maxEdep[key].second>allEdeps[key][recoDetId]) continue;
	      maxEdep[key].first=recoDetId;
	      maxEdep[key].second=allEdeps[key][recoDetId];
	    } 

	  //CLUSTERS
	  key="clus";
	  allEdeps[key] = templEdep;
	  maxEdep[key]=std::pair<uint32_t,float>(0,-1.0);
	  edm::Handle<reco::PFClusterCollection> pfClusters;
	  iEvent.getByLabel(edm::InputTag(pfClustersCollections_[i],""),pfClusters);
	  
	  //get all within DR=0.4
	  for(reco::PFClusterCollection::const_iterator c_it=pfClusters->begin();   
	      c_it!=pfClusters->end(); 
	      c_it++)
	    {
	      float dR=deltaR(c_it->position(),p);
	      if(dR>0.4) continue;
	      
	      nClusters_[key][i]++;
	      for( const auto& rhf : c_it->hitsAndFractions() ) {
		uint32_t recoDetId( rhf.first.rawId() );
		float recEnFracClustered=rhf.second;
		if(recoDetIdMap.find(recoDetId)==recoDetIdMap.end()) continue;
		nHitsInClusters_[key][i]++;
		float recHitEn=(*recHits)[ recoDetIdMap[recoDetId] ].energy();
		
		//rec hits : convert energy to keV
		float eclus(recHitEn*recEnFracClustered*1e6/mipEn_[i]);
		if(allEdeps[key].find(recoDetId)==allEdeps[key].end())  allEdeps[key][recoDetId] = 0;
		allEdeps[key][recoDetId] += eclus;
		
		//check if maximum found
		if(maxEdep[key].second>allEdeps[key][recoDetId]) continue;
		maxEdep[key].first=recoDetId;
		maxEdep[key].second=allEdeps[key][recoDetId];
	      }
	    }
	  
	  //IN PF CANDIDATES
	  key="pf";
	  for(size_t ipf=0; ipf<selPFs.size(); ipf++)
	    {
	      const reco::PFCandidate &cand=(*pflow)[ selPFs[ipf] ];
	      const reco::PFCandidate::ElementsInBlocks&einb=cand.elementsInBlocks();
	      for(size_t ieleinb=0; ieleinb<einb.size(); ieleinb++)
		{
		  const reco::PFBlockRef blockRef = einb[ieleinb].first;
		  const edm::OwnVector< reco::PFBlockElement > &eleList=blockRef->elements();
		  for(unsigned int iEle=0; iEle<eleList.size(); iEle++)
		    {
		      //11-14 (EE/HEF/HEB)
		      reco::PFBlockElement::Type eletype = eleList[iEle].type();
		      if(i==0 && eletype!=reco::PFBlockElement::HGC_ECAL)  continue;
		      if(i==1 && eletype!=reco::PFBlockElement::HGC_HCALF) continue;
		      if(i==2 && eletype!=reco::PFBlockElement::HGC_HCALB) continue;
		      
		      //get the cluster
		      const reco::PFBlockElementCluster *sc = dynamic_cast<const reco::PFBlockElementCluster*>(&(eleList[iEle]));
		      nClusters_[key][i]++;
		      for( const auto& rhf : sc->clusterRef()->hitsAndFractions() )
			{
			  uint32_t recoDetId( rhf.first.rawId() );
			  float recEnFracClustered=rhf.second;
			  if(recoDetIdMap.find(recoDetId)==recoDetIdMap.end() ) continue;
			  nHitsInClusters_[key][i]++;
			  float recHitEn=(*recHits)[ recoDetIdMap[recoDetId] ].energy();
			  
			  //rec hits : convert energy to keV
			  float eclus(recHitEn*recEnFracClustered*1e6/mipEn_[i]);
			  if(allEdeps[key].find(recoDetId)==allEdeps[key].end())  allEdeps[key][recoDetId] = 0;
			  allEdeps[key][recoDetId] += eclus;
			  
			  //check if maximum found
			  if(maxEdep[key].second>allEdeps[key][recoDetId]) continue;
			  maxEdep[key].first=recoDetId;
			  maxEdep[key].second=allEdeps[key][recoDetId];
			}
		    }
		}
	    }
	  	    

	  //now compute the relevant variables for the regression
	  for(std::map<TString, std::map<uint32_t,float> >::iterator stepIt=allEdeps.begin();
	      stepIt!=allEdeps.end();
	      stepIt++)
	    {
	      TString key(stepIt->first);

	      //check if there is any reasonable energy here
	      if(maxEdep[key].second<0.5) continue;

	      //this is the reference
	      const GlobalPoint refPos( std::move( geom->getPosition(maxEdep[key].first) ) );
	      float refX=refPos.x();
	      float refY=refPos.y();
	      float refEta=refPos.eta();
	      float refPhi=refPos.phi();

	      //get the cell size for the reference hit
	      int layer((maxEdep[key].first >> 19) & 0x1f); 
	      std::vector<HGCalDDDConstants::hgtrap>::const_iterator recModIt( dddConst.getFirstModule(true) );
	      std::pair<int,int>  simToReco=dddConst.simToReco(1,layer,false);
	      for(int klay=1; klay<simToReco.second; klay++) recModIt++;
	      float cellSize=recModIt->cellSize;

	      //loop over the hits
	      for(std::map<uint32_t,float>::iterator detIt=stepIt->second.begin();
		  detIt!=stepIt->second.end();
		  detIt++)
		{
		  const GlobalPoint pos( std::move( geom->getPosition(detIt->first) ) );
		  float hitEta=pos.eta();
		  float hitPhi=pos.phi();
		  int idx(abs((pos.x()-refX)/cellSize));
		  int idy(abs((pos.y()-refY)/cellSize));
		  int layerIdx((detIt->first >> 19) & 0x1f); 
		  layerIdx+=(baseLayerIdx-1);
		  float en(detIt->second);

		  nhits_[key][layerIdx]    ++;
		  edeps_[key][layerIdx]    += en;
		  if(idx<=1 && idy<=1) edeps3x3_[key][layerIdx] += en;
		  if(idx<=3 && idy<=3) edeps5x5_[key][layerIdx] += en;
		  emeanPhi_[key][layerIdx] += en*hitPhi;
		  emeanEta_[key][layerIdx] += en*hitEta;
		  sihih_[key][layerIdx]    += en*pow(hitEta-refEta,2);
		  sipip_[key][layerIdx]    += en*pow(deltaPhi(hitPhi,refPhi),2);
		  sipih_[key][layerIdx]    += en*(hitEta-refEta)*(deltaPhi(hitPhi,refPhi));
		}
	    }

	  //increment base layer idx
	  baseLayerIdx += dddConst.layers(true); 
	}
      

      //finalize normalizing by total energy in each layer
      nlay_=baseLayerIdx;
      for(std::map<TString, Float_t * >::iterator keyIt=edeps_.begin();
	  keyIt!=edeps_.end();
	  keyIt++)
	{
	  TString key(keyIt->first);
	  for(Int_t ilay=0; ilay<nlay_; ilay++)
	    {
	      float totalEn=edeps_[key][ilay];
	      if(totalEn<=0) continue;
	      emeanPhi_[key][ilay] /= totalEn;	  
	      emeanEta_[key][ilay] /= totalEn;
	      sihih_[key][ilay]    /= totalEn;
	      sipip_[key][ilay]    /= totalEn;
	      sipih_[key][ilay]    /= totalEn;
	    }
	}
      
      t_->Fill();
    }
}


//define this as a plug-in
DEFINE_FWK_MODULE(HGCSimHitsAnalyzer);
