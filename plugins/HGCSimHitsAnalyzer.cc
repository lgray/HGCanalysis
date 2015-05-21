#include "UserCode/HGCanalysis/plugins/HGCSimHitsAnalyzer.h"
#include "UserCode/HGCanalysis/interface/HGCAnalysisTools.h"

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

#include "FWCore/Utilities/interface/RandomNumberGenerator.h"

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

  //init rand
  edm::Service<edm::RandomNumberGenerator> rngs;
  if ( ! rngs.isAvailable() ) {
    throw cms::Exception("Configuration") << "HGCDigitizer requires the RandomNumberGeneratorService - please add this service or remove the modules that require it";
  }

  tdcReso_ = new CLHEP::RandGauss( rngs->getEngine() );

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
  t_->Branch("genVertexX",  &genVertexX_,  "genVertexX/F");  
  t_->Branch("genVertexY",  &genVertexY_,  "genVertexY/F");  
  t_->Branch("genVertexZ",  &genVertexZ_,  "genVertexZ/F");  
  t_->Branch("hasInteractionBeforeHGC",   &hasInteractionBeforeHGC_, "hasInteractionBeforeHGC/O");
  t_->Branch("hasChargedInteraction",     &hasChargedInteraction_,   "hasChargedInteraction/O");
  t_->Branch("layerShowerStart",          &layerShowerStart_,        "layerShowerStart/I");
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

      emeanTime_[key]    = 0;
      t_->Branch("emeanTime_"+key,      &emeanTime_[key],       "emeanTime_"+key+"/F");

      emeanTime20_[key] = 0;
      t_->Branch("emeanTime20_"+key,  &emeanTime20_[key], "emeanTime20_"+key+"/F");

      emeanTime50_[key] = 0;
      t_->Branch("emeanTime50_"+key,  &emeanTime50_[key], "emeanTime50_"+key+"/F");

      emeanTime80_[key] = 0;
      t_->Branch("emeanTime80_"+key,  &emeanTime80_[key], "emeanTime80_"+key+"/F");

      emeanTime100_[key] = 0;
      t_->Branch("emeanTime100_"+key,  &emeanTime100_[key], "emeanTime100_"+key+"/F");

      emeanTime150_[key] = 0;
      t_->Branch("emeanTime150_"+key,  &emeanTime150_[key], "emeanTime150_"+key+"/F");

      emeanTime200_[key] = 0;
      t_->Branch("emeanTime200_"+key,  &emeanTime200_[key], "emeanTime200_"+key+"/F");

      avgEPerHitEE_[key]    = 0;
      t_->Branch("avgEPerHitEE_"+key,      &avgEPerHitEE_[key],       "avgEPerHitEE_"+key+"/F");

      avgEPerHitHEF_[key]    = 0;
      t_->Branch("avgEPerHitHEF_"+key,      &avgEPerHitHEF_[key],       "avgEPerHitHEF_"+key+"/F");

      avgEPerHitHEB_[key]    = 0;
      t_->Branch("avgEPerHitHEB_"+key,      &avgEPerHitHEB_[key],       "avgEPerHitHEB_"+key+"/F");

      totalX0WgtE_[key]=0;
      t_->Branch("totalX0WgtE_"+key,   &totalX0WgtE_[key],    "totalX0WgtE_"+key+"/F");

      totalLambdaWgtE_[key]=0;
      t_->Branch("totalLambdaWgtE_"+key,   &totalLambdaWgtE_[key],    "totalLambdaWgtE_"+key+"/F");

      totalLength_[key]=0;
      t_->Branch("totalLength_"+key,   &totalLength_[key],    "totalLength_"+key+"/F");

      totalVolumeEE_[key]=0;
      t_->Branch("totalVolumeEE_"+key,   &totalVolumeEE_[key],    "totalVolumeEE_"+key+"/F");

      totalVolumeHEF_[key]=0;
      t_->Branch("totalVolumeHEF_"+key,   &totalVolumeHEF_[key],    "totalVolumeHEF_"+key+"/F");

      totalVolumeHEB_[key]=0;
      t_->Branch("totalVolumeHEB_"+key,   &totalVolumeHEB_[key],    "totalVolumeHEB_"+key+"/F");

      showerStart_[key]=-1;
      t_->Branch("showerStart_"+key,   &showerStart_[key],    "showerStart_"+key+"/I");

      edeps_[key]    = new Float_t[100];
      t_->Branch("edep_"+key,      edeps_[key],    "edep_"+key+"[nlay]/F");

      edepstdc_[key]    = new Float_t[100];
      t_->Branch("edepstdc_"+key,      edepstdc_[key],    "edepstdc_"+key+"[nlay]/F");

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

      nhitstdc_[key]    = new Int_t[100];
      t_->Branch("nhitstdc_"+key,     nhitstdc_[key],    "nhitstdc_"+key+"[nlay]/I");

      nhitsavg_[key]    = new Int_t[100];
      t_->Branch("nhitsavg_"+key,     nhitsavg_[key],    "nhitsavg_"+key+"[nlay]/I");

      nhits5mip_[key]    = new Int_t[100];
      t_->Branch("nhits5mip_"+key,     nhits5mip_[key],    "nhits5mip_"+key+"[nlay]/I");

      nhits10mip_[key]    = new Int_t[100];
      t_->Branch("nhits10mip_"+key,     nhits10mip_[key],    "nhits10mip_"+key+"[nlay]/I");

      emeanTimeLayer_[key] = new Float_t[100];
      t_->Branch("emeanTimeLayer_"+key,  emeanTimeLayer_[key], "emeanTimeLayer_"+key+"[nlay]/F");

      emeanTimeLayer20_[key] = new Float_t[100];
      t_->Branch("emeanTimeLayer20_"+key,  emeanTimeLayer20_[key], "emeanTimeLayer20_"+key+"[nlay]/F");

      emeanTimeLayer50_[key] = new Float_t[100];
      t_->Branch("emeanTimeLayer50_"+key,  emeanTimeLayer50_[key], "emeanTimeLayer50_"+key+"[nlay]/F");

      emeanTimeLayer80_[key] = new Float_t[100];
      t_->Branch("emeanTimeLayer80_"+key,  emeanTimeLayer80_[key], "emeanTimeLayer80_"+key+"[nlay]/F");

      emeanTimeLayer100_[key] = new Float_t[100];
      t_->Branch("emeanTimeLayer100_"+key,  emeanTimeLayer100_[key], "emeanTimeLayer100_"+key+"[nlay]/F");

      emeanTimeLayer150_[key] = new Float_t[100];
      t_->Branch("emeanTimeLayer150_"+key,  emeanTimeLayer150_[key], "emeanTimeLayer150_"+key+"[nlay]/F");

      emeanTimeLayer200_[key] = new Float_t[100];
      t_->Branch("emeanTimeLayer200_"+key,  emeanTimeLayer200_[key], "emeanTimeLayer200_"+key+"[nlay]/F");

      maxTimeLayer_[key] = new Float_t[100];
      t_->Branch("maxTimeLayer_"+key,  maxTimeLayer_[key], "maxTimeLayer_"+key+"[nlay]/F");
      
      maxTimeEnergyLayer_[key] = new Float_t[100];
      t_->Branch("maxTimeEnergyLayer_"+key,  maxTimeEnergyLayer_[key], "maxTimeEnergyLayer_"+key+"[nlay]/F");

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

      widthep1_[key]    = new Float_t[100];
      t_->Branch("widthep1_"+key,     widthep1_[key],    "widthep1_"+key+"[nlay]/F");
      
      widthep2_[key]    = new Float_t[100];
      t_->Branch("widthep2_"+key,     widthep2_[key],    "widthep2_"+key+"[nlay]/F");

      width1_[key]    = new Float_t[100];
      t_->Branch("width1_"+key,     width1_[key],    "width1_"+key+"[nlay]/F");
      
      width2_[key]    = new Float_t[100];
      t_->Branch("width2_"+key,     width2_[key],    "width2_"+key+"[nlay]/F");

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
  edm::Handle<std::vector<SimTrack> > SimTk;
  iEvent.getByLabel(g4TracksSource_,SimTk);
  edm::Handle<std::vector<SimVertex> > SimVtx;
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
  std::vector<float> layerZ;
  for(size_t i=0; i<geometrySource_.size(); i++)  {
    iSetup.get<IdealGeometryRecord>().get(geometrySource_[i],geom[i]);
    const HGCalTopology &topo=geom[i]->topology();
    const HGCalDDDConstants &dddConst=topo.dddConstants(); 
    layerCtrOffset.push_back(nlay_);
    nlay_ += dddConst.layers(true);

    for(size_t ilay=0; ilay<dddConst.layers(true); ilay++)
      {
	uint32_t mySubDet(ForwardSubdetector::HGCEE);
	if(i==1) mySubDet=ForwardSubdetector::HGCHEF;
	else if(i==2) mySubDet=ForwardSubdetector::HGCHEB;
	uint32_t recoDetId = ( (i==0) ?
			       (uint32_t)HGCEEDetId(ForwardSubdetector(mySubDet),0,ilay+1,1,0,1):
			       (uint32_t)HGCHEDetId(ForwardSubdetector(mySubDet),0,ilay+1,1,0,1)
			       );
	
	//require to be on the same side of the generated particle
	const GlobalPoint pos( std::move( geom[i]->getPosition(recoDetId) ) );	  
	layerZ.push_back(fabs(pos.z()));
      }
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
      G4InteractionPositionInfo intInfo=getInteractionPosition(SimTk.product(),SimVtx.product(),genBarcodes->at(igen));
      math::XYZVectorD hitPos=intInfo.pos;
      genHitX_=hitPos.x();
      genHitY_=hitPos.y();
      genHitZ_=hitPos.z();
      hasInteractionBeforeHGC_=(fabs(hitPos.z())<317 && fabs(hitPos.z()) > 1e-3);
      hasChargedInteraction_=intInfo.info;

      // gen vertex positions
      genVertexX_ = p.vx();
      genVertexY_ = p.vy();
      genVertexZ_ = p.vz();

      //match nearest HGC layer in Z
      float dzMin(99999999.);
      for(size_t ilay=0; ilay<layerZ.size(); ilay++)
	{
	  float dz=fabs(fabs(layerZ[ilay])-fabs(genHitZ_));
	  if(dz>dzMin) continue;
	  dzMin=dz;
	  layerShowerStart_=ilay;
	}


      //hits
      std::unordered_map<uint32_t,uint32_t> recHitsIdMap; //needed to decode hits in clusters
      std::map<TString, std::map<uint32_t, std::pair<float,float> > > allEdeps, allCtrlEdeps;
      allEdeps["sim"]  = std::map<uint32_t,std::pair<float,float> >();
      allEdeps["rec"]  = std::map<uint32_t,std::pair<float,float> >();
      allCtrlEdeps["rec"]= std::map<uint32_t,std::pair<float,float> >();
      allEdeps["clus"] = std::map<uint32_t,std::pair<float,float> >();
      allEdeps["pf"]   = std::map<uint32_t,std::pair<float,float> >();
      std::map<TString, std::pair<uint32_t,std::pair<float,float> > > maxEdep;
      maxEdep["sim"]   = std::pair<uint32_t,std::pair<float,float> >(0,std::make_pair(0.f,0.f));
      maxEdep["rec"]   = std::pair<uint32_t,std::pair<float,float> >(0,std::make_pair(0.f,0.f));
      maxEdep["clus"]  = std::pair<uint32_t,std::pair<float,float> >(0,std::make_pair(0.f,0.f));     
      maxEdep["pf"]    = std::pair<uint32_t,std::pair<float,float> >(0,std::make_pair(0.f,0.f));
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
	      if(allEdeps[key].find(recoDetId)==allEdeps[key].end()) allEdeps[key][recoDetId].first = 0;
	      allEdeps[key][recoDetId].first += hitEnInMIPs;
	      
	      //check if maximum found
	      if(maxEdep[key].second.first>allEdeps[key][recoDetId].first) continue;
	      maxEdep[key].first=recoDetId;
	      maxEdep[key].second.first=allEdeps[key][recoDetId].first;
	    }

	  //RECO: save rec hits where sim hits exist or after shifting phi
	  key="rec";
	  uint32_t recHitCtr(0);
	  for(HGCRecHitCollection::const_iterator hit_it=recHits[i]->begin(); hit_it!=recHits[i]->end(); hit_it++,recHitCtr++)
	    {
	      //convert energy to keV
	      float hitEn(hit_it->energy()*1e6/mipEn_[i]);
              float hitTime(hit_it->time());
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
		      if(allCtrlEdeps[key].find(recoDetId)==allCtrlEdeps[key].end()) allCtrlEdeps[key][recoDetId].first = 0;
		      allCtrlEdeps[key][recoDetId].first += hitEn;
                      allCtrlEdeps[key][recoDetId].second = hitTime;
		    }
		  continue;
		}

	      //signal region
	      if(allEdeps[key].find(recoDetId)==allEdeps[key].end()) allEdeps[key][recoDetId].first = 0;
	      allEdeps[key][recoDetId].first += hitEn;
              allEdeps[key][recoDetId].second = hitTime;
	      
	      //check if maximum found
	      if(maxEdep[key].second.first>allEdeps[key][recoDetId].first) continue;
	      maxEdep[key].first=recoDetId;
	      maxEdep[key].second.first=allEdeps[key][recoDetId].first;
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
	    if(allEdeps[key].find(recoDetId)==allEdeps[key].end())  allEdeps[key][recoDetId].first = 0;
	    allEdeps[key][recoDetId].first += eclus;
            allEdeps[key][recoDetId].second = (*(recHits[subDetCtr]))[ recHitsIdMap[recoDetId] ].time();
	    
	    //check if maximum found
	    if(maxEdep[key].second.first>allEdeps[key][recoDetId].first) continue;
	    maxEdep[key].first=recoDetId;
	    maxEdep[key].second.first=allEdeps[key][recoDetId].first;
	  }
	}
       
      
      //RECHITS FROM SUPERCLUSTERS (e/gamma) or CLUSTERS IN  SELECTED PF CANDIDATE
      key="pf";
      if((abs(genId_)==11 || abs(genId_)==22) && pfClustersCollection_.find("pandora")==std::string::npos)
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
		      if(allEdeps[key].find(recoDetId)==allEdeps[key].end())  allEdeps[key][recoDetId].first = 0;
		      allEdeps[key][recoDetId].first += eclus;
                      allEdeps[key][recoDetId].second += (*(recHits[subDetCtr]))[ recHitsIdMap[recoDetId] ].time();
		      
		      //check if maximum found
		      if(maxEdep[key].second.first>allEdeps[key][recoDetId].first) continue;
		      maxEdep[key].first=recoDetId;
		      maxEdep[key].second.first=allEdeps[key][recoDetId].first;
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
		      if(allEdeps[key].find(recoDetId)==allEdeps[key].end())  allEdeps[key][recoDetId].first = 0;
		      allEdeps[key][recoDetId].first += eclus;
			  
		      //check if maximum found
		      if(maxEdep[key].second.first>allEdeps[key][recoDetId].first) continue;
		      maxEdep[key].first=recoDetId;
		      maxEdep[key].second.first=allEdeps[key][recoDetId].first;
		    }
		}
	}
	}
	
      //now compute the relevant variables for the regression
      int nFailed(0);
      for(std::map<TString, std::map<uint32_t,std::pair<float,float> > >::iterator stepIt=allEdeps.begin();
	  stepIt!=allEdeps.end();
	  stepIt++)
	{
	  TString key(stepIt->first);

	  //check if there is any reasonable energy for the max. hit
	  if(maxEdep[key].second.first > 0.5)
	    {
	      int subDetId((maxEdep[key].first >>25)&0x7);
	      int subDetCtr(0);
	      if(subDetId==ForwardSubdetector::HGCHEF) subDetCtr=1;
	      if(subDetId==ForwardSubdetector::HGCHEB) subDetCtr=2;
	      
	      hitMaxLayer_[key] = ((maxEdep[key].first >> 19) & 0x1f) + layerCtrOffset[subDetCtr]-1; 
	      hitMax_[key]      = maxEdep[key].second.first;
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
	  
          // useful constants
          //constexpr float cm_per_ns = 29.9792458;

	  //loop over the hits
	  int nhitsEE(0), nhitsHEF(0), nhitsHEB(0);
	  float totalEn(0), totalEnTDC(0),totalX0WgtEn(0), totalLambdaWgtEn(0),totalEnEE(0),totalEnHEF(0),totalEnHEB(0);
          
          //constexpr float tdc_reso = 0.080f;

	  for(std::map<uint32_t,std::pair<float,float> >::iterator detIt=stepIt->second.begin();
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
	       
	    
              constexpr std::array<float,6> tdc_resolutions = { { 0.020f, 0.050f, 0.080f, 0.100f, 0.150f, 0.200f} } ;

	      int layerIdx(hitLayer+layerCtrOffset[subDetCtr]-1); 
	      float en(detIt->second.first);
              float time(detIt->second.second);
              float time_smeared[6];// 20,50,80,100,150,200
              for( int ii = 0; ii < 6; ++ii ) {
                time_smeared[ii] = tdcReso_->fire(time,tdc_resolutions[ii]);
              }
              
	      totalEn                    += en;
              totalEnTDC                 += ( time > 0 ? en : 0.0 );
	      if(subDetCtr==0)           { totalEnEE += en; nhitsEE++; }
	      if(subDetCtr==1)           { totalEnHEF += en; nhitsHEF++; }
	      if(subDetCtr==2)           { totalEnHEB += en; nhitsHEB++; }
	      totalX0WgtEn               += en*getLayerWeight(layerIdx,true);
	      totalLambdaWgtEn           += en*getLayerWeight(layerIdx,false);
	      edeps_[key][layerIdx]         += en;
              
	      nhits_[key][layerIdx]      +=1;
	      nhits5mip_[key][layerIdx]  +=1*(en>5);
	      nhits10mip_[key][layerIdx] +=1*(en>10.);
	      emeanPhi_[key][layerIdx]      += en*hitPhi;
	      showerMeanPhi_[key]           += en*hitPhi;
              // correct the time back to the front face of the calorimeter
              if( time > 0 ) {
                /*
                std::cout << stepIt->first << ' ' << layerIdx << ' ' << time << ' ' << en << ' ' 
                          << (std::abs(hitZ) - std::abs(layerZ[0]))/cm_per_ns << ' ' << std::abs(hitZ) << ' ' 
                          << std::abs(layerZ[0]) << ' ' 
                          << time - (std::abs(hitZ) - std::abs(layerZ[0]))/cm_per_ns << std::endl;
                */
                emeanTimeLayer_[key][layerIdx] += en*( time ); // - (std::abs(hitZ) - std::abs(layerZ[0]))/cm_per_ns
                
                emeanTimeLayer20_[key][layerIdx] += en*( time_smeared[0] );
                emeanTimeLayer50_[key][layerIdx] += en*( time_smeared[1] );
                emeanTimeLayer80_[key][layerIdx] += en*( time_smeared[2] );
                emeanTimeLayer100_[key][layerIdx] += en*( time_smeared[3] );
                emeanTimeLayer150_[key][layerIdx] += en*( time_smeared[4] );
                emeanTimeLayer200_[key][layerIdx] += en*( time_smeared[5] );

                edepstdc_[key][layerIdx]      += en;

                if( time > maxTimeLayer_[key][layerIdx] ) {
                  maxTimeLayer_[key][layerIdx] = time;
                  maxTimeEnergyLayer_[key][layerIdx] = en;
                }

                ++nhitstdc_[key][layerIdx];
              }
              
	      emeanEta_[key][layerIdx]      += en*hitEta;
	      showerMeanEta_[key]           += en*hitEta;
	      emeanX_[key][layerIdx]        += en*hitX;
	      showerMeanX_[key]             += en*hitX;
	      emeanY_[key][layerIdx]        += en*hitY;
	      showerMeanY_[key]             += en*hitY;
	      showerMeanZ_[key]             += en*hitZ;
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
              float iTDCEn( edepstdc_[key][ilay] );
              if( iTDCEn > 0.001f ) {                
                emeanTime_[key] += emeanTimeLayer_[key][ilay];
                emeanTime20_[key] += emeanTimeLayer20_[key][ilay];
                emeanTime50_[key] += emeanTimeLayer50_[key][ilay];
                emeanTime80_[key] += emeanTimeLayer80_[key][ilay];
                emeanTime100_[key] += emeanTimeLayer100_[key][ilay];
                emeanTime150_[key] += emeanTimeLayer150_[key][ilay];
                emeanTime200_[key] += emeanTimeLayer200_[key][ilay];
                
                emeanTimeLayer_[key][ilay] /= iTDCEn; 
                emeanTimeLayer20_[key][ilay] /= iTDCEn; 
                emeanTimeLayer50_[key][ilay] /= iTDCEn; 
                emeanTimeLayer80_[key][ilay] /= iTDCEn; 
                emeanTimeLayer100_[key][ilay] /= iTDCEn; 
                emeanTimeLayer150_[key][ilay] /= iTDCEn; 
                emeanTimeLayer200_[key][ilay] /= iTDCEn; 
              } else {
                emeanTimeLayer_[key][ilay] = -1.f;
                emeanTimeLayer20_[key][ilay] = -1.f;
                emeanTimeLayer50_[key][ilay] = -1.f;
                emeanTimeLayer80_[key][ilay] = -1.f;
                emeanTimeLayer100_[key][ilay] = -1.f;
                emeanTimeLayer150_[key][ilay] = -1.f;
                emeanTimeLayer200_[key][ilay] = -1.f;
              }
	    }

          if( totalEnTDC > 0.001f) {
            emeanTime_[key] /= totalEnTDC;     
            emeanTime20_[key] /= totalEnTDC;
            emeanTime50_[key] /= totalEnTDC;
            emeanTime80_[key] /= totalEnTDC;
            emeanTime100_[key] /= totalEnTDC;
            emeanTime150_[key] /= totalEnTDC;
            emeanTime200_[key] /= totalEnTDC;
          } else {
            emeanTime_[key] = -1.f;
            emeanTime20_[key] = -1.f;
            emeanTime50_[key] = -1.f;
            emeanTime80_[key] = -1.f;
            emeanTime100_[key] = -1.f;
            emeanTime150_[key] = -1.f;
            emeanTime200_[key] = -1.f;
          }

	  //save energy sums
	  totalE_[key]=totalEn;
	  totalX0WgtE_[key]=totalX0WgtEn;
	  totalLambdaWgtE_[key]=totalLambdaWgtEn;
	  avgEPerHitEE_[key]  = nhitsEE> 0 ? totalEnEE/nhitsEE : 0.;
	  avgEPerHitHEF_[key] = nhitsHEF> 0 ? totalEnHEF/nhitsHEF : 0.;
	  avgEPerHitHEB_[key] = nhitsHEB> 0 ? totalEnHEB/nhitsHEB : 0.;

	  //loop once more over hits to determine distance to shower direction
	  std::vector<float> sihih(nlay_,0), sipip(nlay_,0), sihip(nlay_,0);
	  std::vector<float> sixix(nlay_,0), siyiy(nlay_,0), sixiy(nlay_,0);
	  std::vector<float> en2(nlay_,0);
	  for(std::map<uint32_t,std::pair<float,float> >::iterator detIt=stepIt->second.begin();
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

	      float en(detIt->second.first);
	      //int idx(abs((hitX-emeanX_[key][layerIdx])/cellSize));
	      //int idy(abs((hitY-emeanY_[key][layerIdx])/cellSize));
	      int idx(abs((hitX-emeanX_[key][layerIdx])/cellSize));
	      int idy(abs((hitY-emeanY_[key][layerIdx])/cellSize));
	      if(idx<=1 && idy<=1) edeps3x3_[key][layerIdx] += en;
	      if(idx<=3 && idy<=3) edeps5x5_[key][layerIdx] += en;
	      
	      float refRho(showerMeanEta_[key]!=0 ? fabs(hitZ/TMath::SinH(showerMeanEta_[key])) : 0. );
	      float refX( refRho*TMath::Cos( showerMeanPhi_[key] ) );
	      float refY( refRho*TMath::Sin( showerMeanPhi_[key] ) );
	      
	      float rho=sqrt(pow(hitX-refX,2)+pow(hitY-refY,2));
	      float avgToUse(avgEPerHitEE_[key]);
	      if(subDetCtr==1) avgToUse=avgEPerHitHEF_[key];
	      if(subDetCtr==2) avgToUse=avgEPerHitHEB_[key];
	      nhitsavg_[key][layerIdx] += (en>avgToUse);
	      edepdR_[key][layerIdx]   += en*rho;        //en*rho;
	      edepArea_[key][layerIdx] += cellSize;
	      sihih[layerIdx]          += pow(en*(hitEta-emeanEta_[key][layerIdx]),2);
	      sipip[layerIdx]          += pow(en*deltaPhi(hitPhi,emeanPhi_[key][layerIdx]),2);
	      sihip[layerIdx]          += -en*en*(hitEta-emeanEta_[key][layerIdx])*(deltaPhi(hitPhi,emeanPhi_[key][layerIdx]));
	      sixix[layerIdx]          += pow(en*(hitX-refX),2);
	      siyiy[layerIdx]          += pow(en*(hitY-refY),2);
	      sixiy[layerIdx]          += -en*en*(hitX-refX)*(hitY-refY);
	      en2[layerIdx]            += pow(en,2);
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

	      if(en2[ilay])
		{
		  sihih[ilay] /= en2[ilay];    sipip[ilay] /= en2[ilay];   sihip[ilay] /= en2[ilay];
		  double mvals[4]={sihih[ilay],sihip[ilay],sihip[ilay],sipip[ilay]};
		  TMatrixDSym m(2,mvals);
		  TMatrixDSymEigen me(m);
		  TVectorD eigenval = me.GetEigenValues();
		  widthep1_[key][ilay]=eigenval[0];
		  widthep2_[key][ilay]=eigenval[1];

		  sixix[ilay] /= en2[ilay];    siyiy[ilay] /= en2[ilay];   sixiy[ilay] /= en2[ilay];
		  double mvalsxy[4]={sixix[ilay],sixiy[ilay],sixiy[ilay],siyiy[ilay]};
		  TMatrixDSym mxy(2,mvalsxy);
		  TMatrixDSymEigen mexy(mxy);
		  TVectorD eigenvalxy = mexy.GetEigenValues();
		  width1_[key][ilay]=eigenvalxy[0];
		  width2_[key][ilay]=eigenvalxy[1];
		}
	      else
		{
		  widthep1_[key][ilay]=0; widthep2_[key][ilay]=0;
		  width1_[key][ilay]=0;   width2_[key][ilay]=0;
		}

	      float corrOverburden(getLayerWeight(ilay,false));
	      totalLength_[key] += corrOverburden;

	      //do not use EE or HEB for the volume
	      if(ilay<30)      totalVolumeEE_[key] += (edepArea_[key][ilay]>0 ? edepArea_[key][ilay] : 0 )*corrOverburden;
	      else if(ilay<42) totalVolumeHEF_[key] += (edepArea_[key][ilay]>0 ? edepArea_[key][ilay] : 0 )*corrOverburden;
	      else             totalVolumeHEB_[key] += (edepArea_[key][ilay]>0 ? edepArea_[key][ilay] : 0 )*corrOverburden;
	    }	 
	}

      //rec hits for control region
      for(std::map<uint32_t,std::pair<float,float> >::iterator ctrlIt=allCtrlEdeps["rec"].begin();
	  ctrlIt!= allCtrlEdeps["rec"].end();
	  ctrlIt++)
	{
	  int subDetId((ctrlIt->first>>25)&0x7);
	  int subDetCtr(0);
	  if(subDetId==ForwardSubdetector::HGCHEF) subDetCtr=1;
	  if(subDetId==ForwardSubdetector::HGCHEB) subDetCtr=2;
	  int hitLayer((ctrlIt->first >> 19) & 0x1f);
	  int layerIdx(hitLayer+layerCtrOffset[subDetCtr]-1); 
	  float en(ctrlIt->second.first);
	  ctrlnhits_["rec"][layerIdx]++;
	  ctrledeps_["rec"][layerIdx]+= en;
	}

      //log if failed
      if(nFailed) cout << "Failed to position " << nFailed << " det ids for " << key << " step" << endl;
           
      //that's all folks...
      t_->Fill();
    }
}


//define this as a plug-in
DEFINE_FWK_MODULE(HGCSimHitsAnalyzer);
