#include "UserCode/HGCanalysis/plugins/HGCSimHitsAnalyzer.h"

#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFraction.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"

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
  hitCollections_           = iConfig.getUntrackedParameter< std::vector<std::string> >("hitCollections");
  recHitCollections_        = iConfig.getUntrackedParameter< std::vector<std::string> >("recHitCollections");
  pfClustersCollection_     = iConfig.getUntrackedParameter< std::string >("pfClustersCollection");
  emPFClustersCollection_   = iConfig.getUntrackedParameter< std::string >("emPFClustersCollection");
  geometrySource_           = iConfig.getUntrackedParameter< std::vector<std::string> >("geometrySource");
  mipEn_                    = iConfig.getUntrackedParameter< std::vector<double> >("mipEn");
  pfCandAssociationCone_    = iConfig.getUntrackedParameter< double >("pfCandAssociationCone");
  pfClusterAssociationCone_ = iConfig.getUntrackedParameter< double >("pfClusterAssociationCone");
  g4TracksSource_           = iConfig.getUntrackedParameter<std::string>("g4TracksSource");
  g4VerticesSource_         = iConfig.getUntrackedParameter<std::string>("g4VerticesSource");

  //init tree
  edm::Service<TFileService> fs;
  t_=fs->make<TTree>("HGC","Event Summary");
  t_->Branch("genId",  &genId_,  "genId/I");
  t_->Branch("genEn",  &genEn_,  "genEn/F");  
  t_->Branch("genEta", &genEta_, "genEta/F");
  t_->Branch("genPhi", &genPhi_, "genPhi/F");
  t_->Branch("genHitX",  &genHitX_,  "genHitX/F");  
  t_->Branch("genHitY",  &genHitY_,  "genHitY/F");  
  t_->Branch("genHitZ",  &genHitZ_,  "genHitZ/F");  
  t_->Branch("hasInteractionBeforeHGC",   &hasInteractionBeforeHGC_, "hasInteractionBeforeHGC/O");
  t_->Branch("hasChargedInteraction",     &hasChargedInteraction_,   "hasChargedInteraction/O");
  t_->Branch("nlay",        &nlay_,       "nlay/I");
  t_->Branch("pfMatchId",   &pfMatchId_,  "pfMatchId/I");

  for(size_t istep=0; istep<4; istep++)
    {
      TString key("sim");
      if(istep==1) key="rec";
      if(istep==2) key="clus";
      if(istep==3) key="pf";

      showerMeanX_[key]=0;
      t_->Branch("showerMeanX_"+key,       &showerMeanX_[key],     "showerMeanX_"+key+"/F");
      showerMeanY_[key]=0;
      t_->Branch("showerMeanY_"+key,       &showerMeanY_[key],     "showerMeanY_"+key+"/F");
      showerMeanZ_[key]=0;
      t_->Branch("showerMeanZ_"+key,       &showerMeanZ_[key],     "showerMeanZ_"+key+"/F");
      showerMeanEta_[key]=0;
      t_->Branch("showerMeanEta_"+key,     &showerMeanEta_[key],   "showerMeanEta_"+key+"/F");
      showerMeanPhi_[key]=0;
      t_->Branch("showerMeanPhi_"+key,     &showerMeanPhi_[key],   "showerMeanPhi_"+key+"/F");

      hitMax_[key]=0;
      t_->Branch("hitMax_"+key,        &hitMax_[key],      "hitMax_"+key+"/F");
      hitMaxLayer_[key]=0;
      t_->Branch("hitMaxLayer_"+key,   &hitMaxLayer_[key], "hitMaxLayer_"+key+"/I");
      hitMaxX_[key]=0;
      t_->Branch("hitMaxX_"+key,       &hitMaxX_[key],     "hitMaxX_"+key+"/F");
      hitMaxY_[key]=0;
      t_->Branch("hitMaxY_"+key,       &hitMaxY_[key],     "hitMaxY_"+key+"/F");
      hitMaxEta_[key]=0;
      t_->Branch("hitMaxEta_"+key,     &hitMaxEta_[key],   "hitMaxEta_"+key+"/F");
      hitMaxPhi_[key]=0;
      t_->Branch("hitMaxPhi_"+key,     &hitMaxPhi_[key],   "hitMaxPhi_"+key+"/F");

      totalE_[key]    = 0;
      t_->Branch("totalE_"+key,      &totalE_[key],       "totalE_"+key+"/F");

      totalX0WgtE_[key]=0;
      t_->Branch("totalX0WgtE_"+key,   &totalX0WgtE_[key],    "totalX0WgtE_"+key+"/F");

      totalLambdaWgtE_[key]=0;
      t_->Branch("totalLambdaWgtE_"+key,   &totalLambdaWgtE_[key],    "totalLambdaWgtE_"+key+"/F");

      totalLength_[key]=0;
      t_->Branch("totalLength_"+key,   &totalLength_[key],    "totalLength_"+key+"/F");

      totalVolume_[key]=0;
      t_->Branch("totalVolume_"+key,   &totalVolume_[key],    "totalVolume_"+key+"/F");

      showerStart_[key]=-1;
      t_->Branch("showerStart_"+key,   &showerStart_[key],    "showerStart_"+key+"/I");

      edeps_[key]    = new Float_t[100];
      t_->Branch("edep_"+key,      edeps_[key],    "edep_"+key+"[nlay]/F");

      if(key=="rec")
	{
	  ctrledeps_[key]    = new Float_t[100];
	  t_->Branch("ctrledep_"+key,      ctrledeps_[key],    "ctrledep_"+key+"[nlay]/F");
	  
	  ctrlnhits_[key]    = new Int_t[100];
	  t_->Branch("ctrlnhits_"+key,     ctrlnhits_[key],    "ctrlnhits_"+key+"[nlay]/I");
	}
      
      edeps3x3_[key]    = new Float_t[100];
      t_->Branch("edep3x3_"+key,      edeps3x3_[key],    "edep3x3_"+key+"[nlay]/F");

      edeps5x5_[key]    = new Float_t[100];
      t_->Branch("edep5x5_"+key,      edeps5x5_[key],    "edep5x5_"+key+"[nlay]/F");

      nhits_[key]    = new Int_t[100];
      t_->Branch("nhits_"+key,     nhits_[key],    "nhits_"+key+"[nlay]/I");

      nhits5mip_[key]    = new Int_t[100];
      t_->Branch("nhits5mip_"+key,     nhits5mip_[key],    "nhits5mip_"+key+"[nlay]/I");

      nhits10mip_[key]    = new Int_t[100];
      t_->Branch("nhits10mip_"+key,     nhits10mip_[key],    "nhits10mip_"+key+"[nlay]/I");

      emeanPhi_[key] = new Float_t[100];
      t_->Branch("emeanPhi_"+key,  emeanPhi_[key], "emeanPhi_"+key+"[nlay]/F");

      emeanEta_[key] = new Float_t[100];
      t_->Branch("emeanEta_"+key,  emeanEta_[key], "emeanEta_"+key+"[nlay]/F");
      
      emeanX_[key] = new Float_t[100];
      t_->Branch("emeanX_"+key,  emeanX_[key], "emeanX_"+key+"[nlay]/F");
      
      emeanY_[key] = new Float_t[100];
      t_->Branch("emeanY_"+key,  emeanY_[key], "emeanY_"+key+"[nlay]/F");

      edepdR_[key]    = new Float_t[100];
      t_->Branch("edepdR_"+key,      edepdR_[key],    "edepdR_"+key+"[nlay]/F");
      
      edepArea_[key]    = new Float_t[100];
      t_->Branch("edepArea_"+key,      edepArea_[key],    "edepArea_"+key+"[nlay]/F");

      sihih_[key]    = new Float_t[100];
      t_->Branch("sihih_"+key,     sihih_[key],    "sihih_"+key+"[nlay]/F");
      
      sipip_[key]    = new Float_t[100];
      t_->Branch("sipip_"+key,     sipip_[key],    "sipip_"+key+"[nlay]/F");
      
      sipih_[key]    = new Float_t[100];
      t_->Branch("sipih_"+key,     sipih_[key],    "sipih_"+key+"[nlay]/F");

      nClusters_[key]=0;
      t_->Branch("nClusters_"+key,      &(nClusters_[key]),        "nClusters_"+key+"/I");
      
      clusterEn_[key] = new Float_t[5];
      t_->Branch("clusterEn_"+key,      clusterEn_[key],        "clusterEn_"+key+"[nClusters_"+key+"]/F");
      
      clusterZ_[key] = new Float_t[5];
      t_->Branch("clusterZ_"+key,      clusterZ_[key],        "clusterZ_"+key+"[nClusters_"+key+"]/F");

      clusterEta_[key] = new Float_t[5];
      t_->Branch("clusterEta_"+key,      clusterEta_[key],        "clusterEta_"+key+"[nClusters_"+key+"]/F");

      clusterPhi_[key] = new Float_t[5];
      t_->Branch("clusterPhi_"+key,      clusterPhi_[key],        "clusterPhi_"+key+"[nClusters_"+key+"]/F");
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

  //Geant4 collections
  edm::Handle<edm::SimTrackContainer> SimTk;
  iEvent.getByLabel(g4TracksSource_,SimTk);
  edm::Handle<edm::SimVertexContainer> SimVtx;
  iEvent.getByLabel(g4VerticesSource_,SimVtx); 
  edm::Handle<std::vector<int> > genBarcodes;
  iEvent.getByLabel("genParticles",genBarcodes);  
  
  //SimHits
  std::vector<edm::Handle<edm::PCaloHitContainer> > caloHits(hitCollections_.size());
  for(size_t i=0; i<hitCollections_.size(); i++) iEvent.getByLabel(edm::InputTag("g4SimHits",hitCollections_[i]),caloHits[i]); 

  //RecHits
  std::vector<edm::Handle<HGCRecHitCollection> > recHits(recHitCollections_.size());
  for(size_t i=0; i<recHitCollections_.size(); i++) iEvent.getByLabel(edm::InputTag("HGCalRecHit",recHitCollections_[i]),recHits[i]);

  //PF clusters and candidates
  edm::Handle<reco::PFClusterCollection> pfClusters;
  iEvent.getByLabel(edm::InputTag(pfClustersCollection_,""),pfClusters);
  edm::Handle<reco::SuperClusterCollection> emPFClusters;
  iEvent.getByLabel(edm::InputTag(emPFClustersCollection_,""),emPFClusters);
  edm::Handle<reco::PFCandidateCollection> pflow;
  iEvent.getByLabel("particleFlow",pflow);

  //Geometry
  std::vector< edm::ESHandle<HGCalGeometry> > geom(geometrySource_.size());
  nlay_=0;
  std::vector<int> layerCtrOffset;
  for(size_t i=0; i<geometrySource_.size(); i++)  {
    iSetup.get<IdealGeometryRecord>().get(geometrySource_[i],geom[i]);
    const HGCalTopology &topo=geom[i]->topology();
    const HGCalDDDConstants &dddConst=topo.dddConstants(); 
    layerCtrOffset.push_back(nlay_);
    nlay_ += dddConst.layers(true); 
  }

  //ready to roll and loop over generator level particles
  for(size_t igen=0; igen<maxGenParts; igen++)
    {
      //start new particle/shower information
      resetCounters();

      //mc truth
      const reco::GenParticle & p = (*genParticles)[igen];
      genId_  = p.pdgId();
      genEn_  = p.energy();
      genEta_ = p.eta();
      genPhi_ = p.phi();

      //match particle flow candidates to gen candidate
      int selPF(-1);
      for(size_t ipf=0; ipf<pflow->size(); ipf++)
	{
	  const reco::PFCandidate &cand=(*pflow)[ipf];
	  float dr=deltaR(p,cand);
	  if(dr>pfCandAssociationCone_) continue;
	  if(selPF<0) selPF=ipf;
	  else if(cand.energy()>(*pflow)[selPF].energy()) selPF=ipf;
	}

      //sim tracks and vertices
      math::XYZVectorD hitPos=getInteractionPosition(p,SimTk,SimVtx,genBarcodes->at(igen),hasChargedInteraction_);
      genHitX_=hitPos.x();
      genHitY_=hitPos.y();
      genHitZ_=hitPos.z();
      hasInteractionBeforeHGC_=(fabs(hitPos.z())<317);
      
      //hits
      std::unordered_map<uint32_t,uint32_t> recHitsIdMap; //needed to decode hits in clusters
      std::map<TString, std::map<uint32_t,float> > allEdeps, allCtrlEdeps;
      allEdeps["sim"]  = std::map<uint32_t,float>();
      allEdeps["rec"]  = std::map<uint32_t,float>();
      allCtrlEdeps["rec"]= std::map<uint32_t,float>();
      allEdeps["clus"] = std::map<uint32_t,float>();
      allEdeps["pf"]   = std::map<uint32_t,float>();
      std::map<TString, std::pair<uint32_t,float> > maxEdep;
      maxEdep["sim"]   = std::pair<uint32_t,float>(0,0);
      maxEdep["rec"]   = std::pair<uint32_t,float>(0,0);
      maxEdep["clus"]  = std::pair<uint32_t,float>(0,0);     
      maxEdep["pf"]    = std::pair<uint32_t,float>(0,0);
      for(size_t i=0; i<hitCollections_.size(); i++)
	{
	  uint32_t mySubDet(ForwardSubdetector::HGCEE);
	  if(i==1) mySubDet=ForwardSubdetector::HGCHEF;
	  else if(i==2) mySubDet=ForwardSubdetector::HGCHEB;

	  const HGCalTopology &topo=geom[i]->topology();
	  const HGCalDDDConstants &dddConst=topo.dddConstants(); 
	 
	  //SIM HITS
	  TString key("sim");
	  for(edm::PCaloHitContainer::const_iterator hit_it = caloHits[i]->begin(); hit_it != caloHits[i]->end(); ++hit_it) 
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
	      
	      //require to be on the same side of the generated particle
	      const GlobalPoint pos( std::move( geom[i]->getPosition(recoDetId) ) );	  
	      float hitEta(pos.eta());
	      if(hitEta*genEta_<0) continue;    
	      
	      //save it if interesting (convert energy to keV)
	      float hitEnInMIPs(hit_it->energy()*1e6/mipEn_[i]);
	      if(hitEnInMIPs<0.5) continue;
	      if(i==2 && hitEnInMIPs<1.0) continue;
	      if(allEdeps[key].find(recoDetId)==allEdeps[key].end()) allEdeps[key][recoDetId] = 0;
	      allEdeps[key][recoDetId] += hitEnInMIPs;
	      
	      //check if maximum found
	      if(maxEdep[key].second>allEdeps[key][recoDetId]) continue;
	      maxEdep[key].first=recoDetId;
	      maxEdep[key].second=allEdeps[key][recoDetId];
	    }

	  //RECO: save rec hits where sim hits exist or after shifting phi
	  key="rec";
	  uint32_t recHitCtr(0);
	  for(HGCRecHitCollection::const_iterator hit_it=recHits[i]->begin(); hit_it!=recHits[i]->end(); hit_it++,recHitCtr++)
	    {
	      //convert energy to keV
	      float hitEn(hit_it->energy()*1e6/mipEn_[i]);
	      if(hitEn<0.5) continue;

	      uint32_t recoDetId(hit_it->id());
	      recHitsIdMap[recoDetId]=recHitCtr;
	      const GlobalPoint pos( std::move( geom[i]->getPosition(recoDetId) ) );	  
	      float dEta(pos.eta()-genEta_);
	      float dPhi(fabs(deltaPhi(pos.phi(),genPhi_)));
	      if(allEdeps["sim"].find(recoDetId)==allEdeps["sim"].end())
		{
		  //control region
		  if(dEta<0.1 && dPhi>TMath::Pi()-0.1 && dPhi<TMath::Pi()+0.1)
		    {
		      if(allCtrlEdeps[key].find(recoDetId)==allCtrlEdeps[key].end()) allCtrlEdeps[key][recoDetId] = 0;
		      allCtrlEdeps[key][recoDetId] += hitEn;
		    }
		  continue;
		}

	      //signal region
	      if(allEdeps[key].find(recoDetId)==allEdeps[key].end()) allEdeps[key][recoDetId] = 0;
	      allEdeps[key][recoDetId] += hitEn;
	      
	      //check if maximum found
	      if(maxEdep[key].second>allEdeps[key][recoDetId]) continue;
	      maxEdep[key].first=recoDetId;
	      maxEdep[key].second=allEdeps[key][recoDetId];
	    }
	}


      //CLUSTERS
      TString key("clus");
      std::vector<const reco::PFCluster *> pToClusters;
      for(reco::PFClusterCollection::const_iterator c_it=pfClusters->begin();   
	  c_it!=pfClusters->end(); 
	  c_it++)
	{
	  float dR=deltaR(c_it->position(),p);
	  if(dR>pfClusterAssociationCone_) continue;
	  pToClusters.push_back(&(*c_it));
	}	 
      
      nClusters_[key]=pToClusters.size();
      sort(pToClusters.begin(),pToClusters.end(),sortClustersByEnergy);
      if(pToClusters.size())
	{
	  for(Int_t iclu=0; iclu<TMath::Min(nClusters_[key],5); iclu++)
	    {
	      clusterEn_[key][iclu]=pToClusters[iclu]->energy();
	      clusterEta_[key][iclu]=pToClusters[iclu]->eta();
	      clusterPhi_[key][iclu]=pToClusters[iclu]->phi();
	      clusterZ_[key][iclu]=pToClusters[iclu]->position().z();
	    }
	  
	  for( const auto& rhf : pToClusters[0]->hitsAndFractions() ) {
	    uint32_t recoDetId( rhf.first.rawId() );
	    float recEnFracClustered=rhf.second;
	    if(recHitsIdMap.find(recoDetId)==recHitsIdMap.end()) continue;
	    
	    int subDetId((recoDetId >>25)&0x7);
	    int subDetCtr(0);
	    if(subDetId==ForwardSubdetector::HGCHEF) subDetCtr=1;
	    if(subDetId==ForwardSubdetector::HGCHEB) subDetCtr=2;
	    
	    //rec hits : convert energy to keV
	    float recHitEn=(*(recHits[subDetCtr]))[ recHitsIdMap[recoDetId] ].energy();
	    float eclus(recHitEn*recEnFracClustered*1e6/mipEn_[subDetCtr]);
	    if(allEdeps[key].find(recoDetId)==allEdeps[key].end())  allEdeps[key][recoDetId] = 0;
	    allEdeps[key][recoDetId] += eclus;
	    
	    //check if maximum found
	    if(maxEdep[key].second>allEdeps[key][recoDetId]) continue;
	    maxEdep[key].first=recoDetId;
	    maxEdep[key].second=allEdeps[key][recoDetId];
	  }
	}
       
      
      //RECHITS FROM SUPERCLUSTERS (e/gamma) or CLUSTERS IN  SELECTED PF CANDIDATE
      key="pf";
      if(abs(genId_)==11 || abs(genId_)==22)
	{
	  //pick the leading supercluster (there should be only one in principle)
	  const reco::SuperCluster * pToCluster=0;
	  for(reco::SuperClusterCollection::const_iterator c_it=emPFClusters->begin();   
	      c_it!=emPFClusters->end(); 
	      c_it++)
	    {
	      float dR=deltaR(c_it->position(),p);
	      if(dR>pfClusterAssociationCone_) continue;
	      if(pToCluster==0) pToCluster=&(*c_it);
	      else if(pToCluster->energy()<c_it->energy()) pToCluster=&(*c_it);
	    }	 

	  if(pToCluster)
	    {
	      nClusters_[key]     = pToCluster->clustersSize();

	      //supercluster
	      clusterEn_[key][0]  = pToCluster->energy();
	      clusterEta_[key][0] = pToCluster->eta();
	      clusterPhi_[key][0] = pToCluster->phi();
	      clusterZ_[key][0]   = pToCluster->position().z();

	      //seed cluster
	      clusterEn_[key][1]  = pToCluster->seed()->energy();
	      clusterEta_[key][1] = pToCluster->seed()->eta();
	      clusterPhi_[key][1] = pToCluster->seed()->phi();
	      clusterZ_[key][1]   = pToCluster->seed()->position().z();
	      
	      for(reco::CaloCluster_iterator cIt = pToCluster->clustersBegin(); cIt!=pToCluster->clustersEnd(); cIt++)
		{
		  for( const auto& rhf : (*cIt)->hitsAndFractions() ) 
		    {
		      uint32_t recoDetId( rhf.first.rawId() );
		      float recEnFracClustered=rhf.second;
		      if(recHitsIdMap.find(recoDetId)==recHitsIdMap.end()) continue;
		      
		      int subDetId((recoDetId >>25)&0x7);
		      int subDetCtr(-1);
		      if(subDetId==ForwardSubdetector::HGCEE) subDetCtr=0;
		      if(subDetId==ForwardSubdetector::HGCHEF) subDetCtr=1;
		      if(subDetId==ForwardSubdetector::HGCHEB) subDetCtr=2;
		      
		      //rec hits : convert energy to keV
		      float recHitEn=(*(recHits[subDetCtr]))[ recHitsIdMap[recoDetId] ].energy();
		      float eclus(recHitEn*recEnFracClustered*1e6/mipEn_[subDetCtr]);
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
      else if(false) //selPF>=0)
	{
	  const reco::PFCandidate &cand=(*pflow)[ selPF ];
	  pfMatchId_=cand.pdgId();
	  
	  const reco::PFCandidate::ElementsInBlocks&einb=cand.elementsInBlocks();
	  for(size_t ieleinb=0; ieleinb<einb.size(); ieleinb++)
	    {
	      const reco::PFBlockRef blockRef = einb[ieleinb].first;
	      const edm::OwnVector< reco::PFBlockElement > &eleList=blockRef->elements();
	      for(unsigned int iEle=0; iEle<eleList.size(); iEle++)
		{
		  //look only at EE/HEF/HEB
		  reco::PFBlockElement::Type eletype = eleList[iEle].type();
		  if(eletype!=reco::PFBlockElement::HGC_ECAL && eletype!=reco::PFBlockElement::HGC_HCALF && eletype!=reco::PFBlockElement::HGC_HCALB) continue;
		  
		  int subDetCtr(0);
		  if(eletype!=reco::PFBlockElement::HGC_HCALF) subDetCtr=1;
		  if(eletype!=reco::PFBlockElement::HGC_HCALB) subDetCtr=2;
		  
		  //get the cluster
		  const reco::PFBlockElementCluster *sc = dynamic_cast<const reco::PFBlockElementCluster*>(&(eleList[iEle]));
   	          nClusters_[key]++;
	          for( const auto& rhf : sc->clusterRef()->hitsAndFractions() )
		    {
		      uint32_t recoDetId( rhf.first.rawId() );
		      float recEnFracClustered=rhf.second;
		      if(recHitsIdMap.find(recoDetId)==recHitsIdMap.end() ) continue;
		      float recHitEn=(*(recHits[subDetCtr]))[ recHitsIdMap[recoDetId] ].energy();
			  
		      //rec hits : convert energy to keV
		      float eclus(recHitEn*recEnFracClustered*1e6/mipEn_[subDetCtr]);
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
      int nFailed(0);
      for(std::map<TString, std::map<uint32_t,float> >::iterator stepIt=allEdeps.begin();
	  stepIt!=allEdeps.end();
	  stepIt++)
	{
	  TString key(stepIt->first);

	  //check if there is any reasonable energy for the max. hit
	  if(maxEdep[key].second>0.5)
	    {
	      int subDetId((maxEdep[key].first >>25)&0x7);
	      int subDetCtr(0);
	      if(subDetId==ForwardSubdetector::HGCHEF) subDetCtr=1;
	      if(subDetId==ForwardSubdetector::HGCHEB) subDetCtr=2;
	      
	      hitMaxLayer_[key] = ((maxEdep[key].first >> 19) & 0x1f) + layerCtrOffset[subDetCtr]-1; 
	      hitMax_[key]      = maxEdep[key].second;
	      try{
		const GlobalPoint refPos( std::move( geom[subDetCtr]->getPosition(maxEdep[key].first) ) );
		hitMaxX_[key]   = refPos.x();
		hitMaxY_[key]   = refPos.y();
		hitMaxEta_[key] = refPos.eta();
		hitMaxPhi_[key] = refPos.phi();
	      }
	      catch(...){
		nFailed++;
	      }
	    }
	  
	  //loop over the hits
	  float totalEn(0),totalX0WgtEn(0), totalLambdaWgtEn(0);
	  for(std::map<uint32_t,float>::iterator detIt=stepIt->second.begin();
	      detIt!=stepIt->second.end();
	      detIt++)
	    {

	      int subDetId((detIt->first>>25)&0x7);
	      int subDetCtr(0);
	      if(subDetId==ForwardSubdetector::HGCHEF) subDetCtr=1;
	      if(subDetId==ForwardSubdetector::HGCHEB) subDetCtr=2;
	      int hitLayer((detIt->first >> 19) & 0x1f);
	      float hitEta(0),hitPhi(0),hitX(0),hitY(0),hitZ(0);
	      try{
		const GlobalPoint pos( std::move( geom[subDetCtr]->getPosition(detIt->first) ) );
		hitX=pos.x();
		hitY=pos.y();
		hitZ=pos.z();
		hitEta=pos.eta();
		hitPhi=pos.phi();
	      }catch(...){
		nFailed++; continue;
	      }
	       
	    
	      int layerIdx(hitLayer+layerCtrOffset[subDetCtr]-1); 
	      float en(detIt->second);
	      totalEn                    += en;
	      totalX0WgtEn               += en*getLayerWeight(layerIdx,true);
	      totalLambdaWgtEn           += en*getLayerWeight(layerIdx,false);
	      edeps_[key][layerIdx]      += en;
	      nhits_[key][layerIdx]      +=1;
	      nhits5mip_[key][layerIdx]  +=1*(en>5);
	      nhits10mip_[key][layerIdx] +=1*(en>10.);
	      emeanPhi_[key][layerIdx] += en*hitPhi;
	      showerMeanPhi_[key]      += en*hitPhi;
	      emeanEta_[key][layerIdx] += en*hitEta;
	      showerMeanEta_[key]      += en*hitEta;
	      emeanX_[key][layerIdx]   += en*hitX;
	      showerMeanX_[key]        += en*hitX;
	      emeanY_[key][layerIdx]   += en*hitY;
	      showerMeanY_[key]        += en*hitY;
	      showerMeanZ_[key]        += en*hitZ;
	    }

	  //compute shower direction (global and local)
	  if(totalEn<=0 || isnan(totalEn) || totalEn>1e10) continue;
	  showerMeanPhi_[key]/=totalEn;
	  showerMeanEta_[key]/=totalEn;
	  showerMeanX_[key]/=totalEn;
	  showerMeanY_[key]/=totalEn;
	  showerMeanZ_[key]/=totalEn;
	  for(Int_t ilay=0; ilay<nlay_; ilay++)
	    {
	      float iTotalEn( edeps_[key][ilay] );
	      if(iTotalEn<=0) continue;
	      emeanPhi_[key][ilay] /= iTotalEn;	  
	      emeanEta_[key][ilay] /= iTotalEn;
	      emeanX_[key][ilay]   /= iTotalEn;	  
	      emeanY_[key][ilay]   /= iTotalEn;
	    }

	  //save energy sums
	  totalE_[key]=totalEn;
	  totalX0WgtE_[key]=totalX0WgtEn;
	  totalLambdaWgtE_[key]=totalLambdaWgtEn;

	  //loop once more over hits to determine distance to shower direction
	  for(std::map<uint32_t,float>::iterator detIt=stepIt->second.begin();
	      detIt!=stepIt->second.end();
	      detIt++)
	    {
	      
	      int subDetId((detIt->first>>25)&0x7);
	      int subDetCtr(0);
	      if(subDetId==ForwardSubdetector::HGCHEF) subDetCtr=1;
	      if(subDetId==ForwardSubdetector::HGCHEB) subDetCtr=2;
	      int hitLayer((detIt->first >> 19) & 0x1f);
	      float hitEta(0),hitPhi(0),hitX(0),hitY(0),hitZ(0);
	      try{
		const GlobalPoint pos( std::move( geom[subDetCtr]->getPosition(detIt->first) ) );
		hitX=pos.x();
		hitY=pos.y();
		hitZ=pos.z();
		hitEta=pos.eta();
		hitPhi=pos.phi();
	      }catch(...){
	      }

	      const HGCalTopology &topo=geom[subDetCtr]->topology();
	      const HGCalDDDConstants &dddConst=topo.dddConstants(); 
	      std::vector<HGCalDDDConstants::hgtrap>::const_iterator recModIt( dddConst.getFirstModule(true) );
	      std::pair<int,int>  simToReco=dddConst.simToReco(1,hitLayer,false);
	      for(int klay=1; klay<simToReco.second; klay++) recModIt++;
	      float cellSize=recModIt->cellSize;
	      int layerIdx(hitLayer+layerCtrOffset[subDetCtr]-1);	      

	      float en(detIt->second);
	      int idx(abs((hitX-emeanX_[key][layerIdx])/cellSize));
	      int idy(abs((hitY-emeanY_[key][layerIdx])/cellSize));
	      if(idx<=1 && idy<=1) edeps3x3_[key][layerIdx] += en;
	      if(idx<=3 && idy<=3) edeps5x5_[key][layerIdx] += en;
	      
	      float refRho(showerMeanEta_[key]!=0 ? fabs(hitZ/TMath::SinH(showerMeanEta_[key])) : 0. );
	      float refX( refRho*TMath::Cos( showerMeanPhi_[key] ) );
	      float refY( refRho*TMath::Sin( showerMeanPhi_[key] ) );
	      
	      float rho=sqrt(pow(hitX-refX,2)+pow(hitY-refY,2));
	      edepdR_[key][layerIdx]   += en*rho;
	      edepArea_[key][layerIdx]  += en*pow(rho,2);
	      sihih_[key][layerIdx]    += en*pow(hitEta-emeanEta_[key][layerIdx],2);
	      sipip_[key][layerIdx]    += en*pow(deltaPhi(hitPhi,emeanPhi_[key][layerIdx]),2);
	      sipih_[key][layerIdx]    += en*(hitEta-emeanEta_[key][layerIdx])*(deltaPhi(hitPhi,emeanPhi_[key][layerIdx]));
	    }
	  
	  //finalize by normalizing
	  int hitLayCtr(0);
	  for(Int_t ilay=0; ilay<nlay_; ilay++)
	    {
	      float iTotalEn( edeps_[key][ilay] );
	      if(iTotalEn<=0) continue;
	      hitLayCtr++;
	      if(showerStart_[key]<0 && hitLayCtr==3) showerStart_[key]=ilay;
	      edepdR_[key][ilay]   /= iTotalEn;
	      edepArea_[key][ilay] /= iTotalEn;
	      edepArea_[key][ilay] = TMath::Pi()*(edepArea_[key][ilay]-pow(edepdR_[key][ilay],2));
	      //if(edepArea_[key][ilay]<0) edepArea_[key][ilay]=0;

	      sihih_[key][ilay]    /= iTotalEn; sihih_[key][ilay]=sqrt(sihih_[key][ilay]);
	      sipip_[key][ilay]    /= iTotalEn; sipip_[key][ilay]=sqrt(sipip_[key][ilay]);
	      sipih_[key][ilay]    /= iTotalEn; sipih_[key][ilay]=sqrt(fabs(sipih_[key][ilay]));
	      
	      float corrOverburden(getLayerWeight(ilay,false)*TMath::TanH(fabs(showerMeanEta_[key])));
	      totalLength_[key] += corrOverburden;

	      //do not use HEB for the volume
	      if(ilay>=41) continue;
	      totalVolume_[key] += (edepArea_[key][ilay]>0 ? edepArea_[key][ilay] : 0 )*corrOverburden;
	    }	 
	}

      //rec hits for control region
      for(std::map<uint32_t,float>::iterator ctrlIt=allCtrlEdeps["rec"].begin();
	  ctrlIt!= allCtrlEdeps["rec"].end();
	  ctrlIt++)
	{
	  int subDetId((ctrlIt->first>>25)&0x7);
	  int subDetCtr(0);
	  if(subDetId==ForwardSubdetector::HGCHEF) subDetCtr=1;
	  if(subDetId==ForwardSubdetector::HGCHEB) subDetCtr=2;
	  int hitLayer((ctrlIt->first >> 19) & 0x1f);
	  int layerIdx(hitLayer+layerCtrOffset[subDetCtr]-1); 
	  float en(ctrlIt->second);
	  ctrlnhits_["rec"][layerIdx]++;
	  ctrledeps_["rec"][layerIdx]+= en;
	}

      //log if failed
      if(nFailed) cout << "Failed to position " << nFailed << " det ids for " << key << " step" << endl;
           
      //that's all folks...
      t_->Fill();
    }
}


//
math::XYZVectorD HGCSimHitsAnalyzer::getInteractionPosition(const reco::GenParticle & genp,
							    edm::Handle<edm::SimTrackContainer> &SimTk,
							    edm::Handle<edm::SimVertexContainer> &SimVtx,
							    int barcode,
							    bool &chargedInteraction)
{
  //loop over vertices
  for (const SimVertex &simVtx : *SimVtx) 
    {
      //require the parent to be the given barcode
      bool noParent( simVtx.noParent() );
      if(noParent) continue;
      int pIdx( simVtx.parentIndex() );
      if( pIdx!=barcode) continue;

      int vtxIdx(simVtx.vertexId());
      int rawTkMult(0),eTkMult(0),gTkMult(0),nTkMult(0),nucleiTkMult(0),pTkMult(0);
      for (const SimTrack &vtxTk : *SimTk)
	{
	  int tkVtxIdx( vtxTk.vertIndex() ); 
	  if(tkVtxIdx!=vtxIdx) continue;
	  
	  int tkType=vtxTk.type();
	  rawTkMult++;
	  eTkMult      += (abs(tkType)==11);
	  gTkMult      += (abs(tkType)==22);
	  nTkMult      += (abs(tkType)==2112 || abs(tkType)==2212);
	  nucleiTkMult += (abs(tkType)>1000000000);
	  pTkMult      += (abs(tkType)==211 || abs(tkType)==111);
	}
      
      if(rawTkMult<2) continue;

      if(rawTkMult==3 && nTkMult==1 && nucleiTkMult==1 && pTkMult==1) chargedInteraction=true;
      return math::XYZVectorD(simVtx.position());

    }

  return math::XYZVectorD(0,0,0);
}


//define this as a plug-in
DEFINE_FWK_MODULE(HGCSimHitsAnalyzer);
