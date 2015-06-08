#include "UserCode/HGCanalysis/plugins/HGCROIAnalyzer.h"
#include "UserCode/HGCanalysis/interface/JetTools.h"
#include "UserCode/HGCanalysis/interface/HGCAnalysisTools.h"
#include "UserCode/HGCanalysis/interface/PCAShowerAnalysis.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFraction.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimG4CMS/Calo/interface/CaloHitID.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"

#include "DetectorDescription/OfflineDBLoader/interface/GeometryInfoDump.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/FCalGeometry/interface/HGCalGeometry.h"

#include <iostream>

using namespace std;

//
HGCROIAnalyzer::HGCROIAnalyzer( const edm::ParameterSet &iConfig ) : 
  slimmedRecHits_(new std::vector<SlimmedRecHit>),
  slimmedClusters_(new std::vector<SlimmedCluster>),
  slimmedROIs_(new std::vector<SlimmedROI>),
  slimmedVertices_(new std::vector<SlimmedVertex>),
  genVertex_(new TVector3),
  useSuperClustersAsROIs_(false)
{
  //configure analyzer
  g4TracksSource_   = iConfig.getUntrackedParameter< std::string >("g4TracksSource");
  g4VerticesSource_ = iConfig.getUntrackedParameter< std::string >("g4VerticesSource");
  genSource_        = iConfig.getUntrackedParameter< std::string >("genSource");
  genJetsSource_    = iConfig.getUntrackedParameter<std::string>("genJetsSource");
  recoVertexSource_ = iConfig.getUntrackedParameter<std::string>("recoVertexSource");
  useSuperClustersAsROIs_ = iConfig.getUntrackedParameter<bool>("useSuperClustersAsROIs");
  superClustersSource_ = iConfig.getUntrackedParameter<std::string>("superClustersSource");
  pfJetsSource_     = iConfig.getUntrackedParameter< std::string >("pfJetsSource");
  eeSimHitsSource_  = iConfig.getUntrackedParameter< std::string >("eeSimHitsSource");
  hefSimHitsSource_ = iConfig.getUntrackedParameter< std::string >("hefSimHitsSource");
  eeRecHitsSource_  = iConfig.getUntrackedParameter< std::string >("eeRecHitsSource");
  hefRecHitsSource_ = iConfig.getUntrackedParameter< std::string >("hefRecHitsSource");

  edm::Service<TFileService> fs;
  tree_=fs->make<TTree>("HGC","HGC");
  tree_->Branch("run",   &run_,   "run/I");
  tree_->Branch("event", &event_, "event/I");
  tree_->Branch("lumi",  &lumi_,  "lumi/I");
  tree_->Branch("RecHits",   "std::vector<SlimmedRecHit>",   &slimmedRecHits_);
  tree_->Branch("Clusters",  "std::vector<SlimmedCluster>",  &slimmedClusters_);
  tree_->Branch("ROIs",      "std::vector<SlimmedROI>",      &slimmedROIs_);
  tree_->Branch("Vertices",  "std::vector<SlimmedVertex>",   &slimmedVertices_);
  tree_->Branch("GenVertex", "TVector3",                     &genVertex_);
}

//
HGCROIAnalyzer::~HGCROIAnalyzer()
{
}

//
//store basic information on RecHits
//
void HGCROIAnalyzer::slimRecHits(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  slimmedRecHits_->clear();
  
  //EE hits
  edm::Handle<edm::PCaloHitContainer> eeSimHits;
  iEvent.getByLabel(edm::InputTag("g4SimHits",eeSimHitsSource_), eeSimHits);
  edm::Handle<HGCRecHitCollection> eeRecHits;
  iEvent.getByLabel(edm::InputTag("HGCalRecHit",eeRecHitsSource_),eeRecHits); 
  if(eeRecHits.isValid())
    {
      edm::ESHandle<HGCalGeometry> eeGeom;
      iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive",eeGeom);
      const HGCalTopology &topo=eeGeom->topology();
      const HGCalDDDConstants &dddConst=topo.dddConstants(); 
      
      float eeMipEn(55.1);
      for(HGCRecHitCollection::const_iterator hit_it=eeRecHits->begin(); 
	  hit_it!=eeRecHits->end(); 
	  hit_it++)
	{
	  uint32_t recoDetId(hit_it->id());
	  const GlobalPoint pos( std::move( eeGeom->getPosition(recoDetId) ) );	
	  slimmedRecHits_->push_back( SlimmedRecHit(recoDetId,
						    pos.x(),pos.y(),pos.z(),
						    hit_it->energy()*1e6/eeMipEn,
						    hit_it->time(),
						    dddConst.getFirstModule(true)->cellSize ) );
	}
      
      //add the simHits
      if(eeSimHits.isValid())
	{
	  for(edm::PCaloHitContainer::const_iterator hit_it = eeSimHits->begin(); 
	      hit_it != eeSimHits->end();
	      hit_it++)
	    {
	      //gang SIM->RECO cells to get final layer assignment  
	      HGCalDetId simId(hit_it->id());
	      int layer(simId.layer()),cell(simId.cell());
	      std::pair<int,int> recoLayerCell=dddConst.simToReco(cell,layer,topo.detectorType());
	      cell  = recoLayerCell.first;
	      layer = recoLayerCell.second;
	      if(layer<0) continue;
	      
	      uint32_t recoDetId( (uint32_t)HGCEEDetId(ForwardSubdetector(ForwardSubdetector::HGCEE),
						       simId.zside(),
						       layer,
						       simId.sector(),
						       simId.subsector(),
						       cell));
	      SlimmedRecHitCollection::iterator theHit=std::find(slimmedRecHits_->begin(),
								 slimmedRecHits_->end(),
								 SlimmedRecHit(recoDetId));
	      if(theHit == slimmedRecHits_->end()) continue;

	      float dist2center( sqrt( theHit->x_*theHit->x_+theHit->y_*theHit->y_+theHit->z_*theHit->z_) );
	      float tof(hit_it->time()-dist2center/(0.1*CLHEP::c_light)+1.0);
	      float emf(hit_it->energyEM()/hit_it->energy());
	      theHit->addSimHit( hit_it->energy()*1e6/eeMipEn,tof,emf );
	    }
	}
    }
  
  //HEF hits
  edm::Handle<edm::PCaloHitContainer> hefSimHits;
  iEvent.getByLabel(edm::InputTag("g4SimHits",hefSimHitsSource_), hefSimHits);
  edm::Handle<HGCRecHitCollection> hefRecHits;
  iEvent.getByLabel(edm::InputTag("HGCalRecHit",hefRecHitsSource_),hefRecHits); 
  if(hefRecHits.isValid())
    {
      edm::ESHandle<HGCalGeometry> hefGeom;
      iSetup.get<IdealGeometryRecord>().get("HGCalHESiliconSensitive",hefGeom);
      const HGCalTopology &topo=hefGeom->topology();
      const HGCalDDDConstants &dddConst=topo.dddConstants(); 
      
      float hefMipEn(85.0);
      for(HGCRecHitCollection::const_iterator hit_it=hefRecHits->begin(); 
	  hit_it!=hefRecHits->end(); 
	  hit_it++)
	{
	  uint32_t recoDetId(hit_it->id());
	  const GlobalPoint pos( std::move( hefGeom->getPosition(recoDetId) ) );	
	  slimmedRecHits_->push_back( SlimmedRecHit(recoDetId,
						    pos.x(),pos.y(),pos.z(),
						    hit_it->energy()*1e6/hefMipEn,
						    hit_it->time(),
						    dddConst.getFirstModule(true)->cellSize ) );
	}
      
      //add the simHits
      if(hefSimHits.isValid())
	{
	  for(edm::PCaloHitContainer::const_iterator hit_it = hefSimHits->begin(); 
	      hit_it != hefSimHits->end();
	      hit_it++)
	    {
	      //gang SIM->RECO cells to get final layer assignment  
	      HGCalDetId simId(hit_it->id());
	      int layer(simId.layer()),cell(simId.cell());
	      std::pair<int,int> recoLayerCell=dddConst.simToReco(cell,layer,topo.detectorType());
	      cell  = recoLayerCell.first;
	      layer = recoLayerCell.second;
	      if(layer<0) continue;
	      
	      uint32_t recoDetId( (uint32_t)HGCHEDetId(ForwardSubdetector(ForwardSubdetector::HGCHEF),
						       simId.zside(),
						       layer,
						       simId.sector(),
						       simId.subsector(),
						       cell));
	      SlimmedRecHitCollection::iterator theHit=std::find(slimmedRecHits_->begin(),
								 slimmedRecHits_->end(),
								 SlimmedRecHit(recoDetId));
	      if(theHit == slimmedRecHits_->end()) continue;
	      
	      float dist2center( sqrt( theHit->x_*theHit->x_+theHit->y_*theHit->y_+theHit->z_*theHit->z_) );
	      float tof(hit_it->time()-dist2center/(0.1*CLHEP::c_light)+1.0);
	      float emf(hit_it->energyEM()/hit_it->energy());
	      theHit->addSimHit( hit_it->energy()*1e6/hefMipEn,tof,emf );
	    }
	}
    }
}

//
void HGCROIAnalyzer::doMCJetMatching(edm::Handle<std::vector<reco::PFJet> > &pfJets,
				     edm::Handle<reco::GenJetCollection> &genJets,
				     edm::Handle<edm::View<reco::Candidate> > &genParticles,
				     std::unordered_map<uint32_t,uint32_t> &reco2genJet,
				     std::unordered_map<uint32_t,uint32_t> &genJet2Parton,
				     std::unordered_map<uint32_t,uint32_t> &genJet2Stable)
{
  
  //
  // match gen jets
  // gen particles
  // - find all status2 parton (after showering) within R=0.4
  // - minimize in pT
  // reco matching
  // based on http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2013_125_v3.pdf
  // - give preference to higher pT gen jets first (collections are already ordered)
  // - find reco jet which minimizes deltaR within 0.2 cone
  // - remove matched reco jet from next matches
  // - iterate until all gen jets are matched
  //
  for(size_t j=0; j<genJets->size(); j++)
    {
      const reco::GenJet& genjet=genJets->at(j);
      float pt=genjet.pt();
      float abseta=fabs(genjet.eta());
      if(abseta<1.5 || abseta>3.0) continue;
     
      //gen particle matching
      float minDPt2Stable(99999.),minDpt2Parton(9999.);
      for(size_t i = 0; i < genParticles->size(); ++ i)
	{
	  const reco::GenParticle & p = dynamic_cast<const reco::GenParticle &>( (*genParticles)[i] );
	  float dR(deltaR(p,genjet));
	  if(dR>0.4) continue;
	  float dPt( fabs(p.pt()-pt));

	  if(p.status()==1)
	    {
	      if(dPt>minDPt2Stable) continue;
	      minDPt2Stable=dPt;
	      genJet2Stable[j]=i;
	    }
	  else if(p.status()==2 && ( (abs(p.pdgId())==21 || abs(p.pdgId())<6) ) )
	    {
	      if(dPt>minDpt2Parton) continue;
	      minDpt2Parton=dPt;
	      genJet2Parton[j]=i;
	    }
	}
  
      //reco matching
      float minDR(0.2);
      for(size_t i=0; i<pfJets->size(); i++)
	{
	  const reco::PFJet &jet=pfJets->at(i);
	  float dR=deltaR(jet,genjet);
	  if(dR>minDR) continue;
	  minDR=dR;
	  if(reco2genJet.find(i)!=reco2genJet.end()) continue;
	  reco2genJet[i]=j;
	}
    }
}

void HGCROIAnalyzer::doMCJetMatching(edm::Handle<reco::SuperClusterCollection> &superClusters,
				     edm::Handle<reco::GenJetCollection> &genJets,
				     edm::Handle<edm::View<reco::Candidate> > &genParticles,
				     std::unordered_map<uint32_t,uint32_t> &reco2genJet,
				     std::unordered_map<uint32_t,uint32_t> &genJet2Parton,
				     std::unordered_map<uint32_t,uint32_t> &genJet2Stable)
{
  
  //
  // match gen jets
  // gen particles
  // - find all status2 parton (after showering) within R=0.4
  // - minimize in pT
  // reco matching
  // based on http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2013_125_v3.pdf
  // - give preference to higher pT gen jets first (collections are already ordered)
  // - find reco jet which minimizes deltaR within 0.2 cone
  // - remove matched reco jet from next matches
  // - iterate until all gen jets are matched
  //
  for(size_t j=0; j<genJets->size(); j++)
    {
      const reco::GenJet& genjet=genJets->at(j);
      float pt=genjet.pt();
      float abseta=fabs(genjet.eta());
      if(abseta<1.5 || abseta>3.0) continue;
     
      //gen particle matching
      float minDPt2Stable(99999.),minDpt2Parton(9999.);
      for(size_t i = 0; i < genParticles->size(); ++ i)
	{
	  const reco::GenParticle & p = dynamic_cast<const reco::GenParticle &>( (*genParticles)[i] );
	  float dR(deltaR(p,genjet));
	  if(dR>0.4) continue;
	  float dPt( fabs(p.pt()-pt));

	  if(p.status()==1)
	    {
	      if(dPt>minDPt2Stable) continue;
	      minDPt2Stable=dPt;
	      genJet2Stable[j]=i;
	    }
	  else if(p.status()==2 && ( (abs(p.pdgId())==21 || abs(p.pdgId())<6) ) )
	    {
	      if(dPt>minDpt2Parton) continue;
	      minDpt2Parton=dPt;
	      genJet2Parton[j]=i;
	    }
	}
  
      //reco matching
      float minDR(0.2);
      int i(0);
      for(reco::SuperClusterCollection::const_iterator c_it=superClusters->begin();
	  c_it!=superClusters->end();
	  c_it++,i++)
	{
	  float dR=deltaR(c_it->position(),genjet);
	  if(dR>minDR) continue;
	  minDR=dR;
	  if(reco2genJet.find(i)!=reco2genJet.end()) continue;
	  reco2genJet[i]=j;
	}
    }
}



//
void HGCROIAnalyzer::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{

  //event header
  run_   = iEvent.id().run();
  event_ = iEvent.id().event();
  lumi_  = iEvent.luminosityBlock();

  //parse rec hits
  slimRecHits(iEvent,iSetup);

  //Geant4 collections
  edm::Handle<std::vector<SimTrack> > SimTk;
  iEvent.getByLabel(g4TracksSource_,SimTk);
  edm::Handle<std::vector<SimVertex> > SimVtx;
  iEvent.getByLabel(g4VerticesSource_,SimVtx); 
  edm::Handle<std::vector<int> > genBarcodes;
  iEvent.getByLabel("genParticles",genBarcodes);  

  //PV collection 
  edm::Handle<reco::VertexCollection> vtxH;
  iEvent.getByLabel(recoVertexSource_, vtxH);
  slimmedVertices_->clear();
  std::vector<size_t> selVtx;
  for(size_t iv=0; iv<vtxH->size(); iv++)
    {
      const reco::Vertex &vtx=vtxH->at(iv);
      if(!vtx.isValid()) continue;
      if(vtx.isFake()) continue;
      selVtx.push_back(iv);
      slimmedVertices_->push_back(SlimmedVertex(vtx.nTracks(),vtx.x(),vtx.y(),vtx.z(),vtx.p4().pt(),vtx.normalizedChi2()) );
    }

  //hard process vertex
  edm::Handle<edm::View<reco::Candidate> > genParticles;
  iEvent.getByLabel(edm::InputTag(genSource_), genParticles);
  for(size_t i = 0; i < genParticles->size(); ++ i)
    {
       const reco::GenParticle & p = dynamic_cast<const reco::GenParticle &>( (*genParticles)[i] );
       if(p.status()!=3) continue;
       genVertex_->SetXYZ(p.vx(),p.vy(),p.vz());
       break;
    }
  
  //jet analysis
  edm::Handle<std::vector<reco::PFJet> > pfJets;
  iEvent.getByLabel(edm::InputTag(pfJetsSource_),pfJets);
  edm::Handle<reco::SuperClusterCollection> superClusters;
  iEvent.getByLabel(edm::InputTag(superClustersSource_,""),superClusters);
  edm::Handle<reco::GenJetCollection> genJets;
  iEvent.getByLabel(edm::InputTag(genJetsSource_), genJets);
  std::unordered_map<uint32_t,uint32_t> reco2genJet,genJet2Parton,genJet2Stable;
  if(useSuperClustersAsROIs_) doMCJetMatching(superClusters,genJets,genParticles,reco2genJet,genJet2Parton,genJet2Stable);
  else                        doMCJetMatching(pfJets,genJets,genParticles,reco2genJet,genJet2Parton,genJet2Stable);

  //
  // Analyze reco jets fiducial in HGC
  // 
  slimmedROIs_->clear();
  slimmedClusters_->clear();

  if(useSuperClustersAsROIs_)
    {
      int j(0);
      for(reco::SuperClusterCollection::const_iterator c_it=superClusters->begin();
	  c_it!=superClusters->end();
	  c_it++,j++)
	{
	  	  
	  if(c_it->energy()<10 || fabs(c_it->eta())<1.5 || fabs(c_it->eta())>3.0) continue;
	  
	  SlimmedROI slimSuperCluster(c_it->energy()/TMath::CosH(c_it->eta()),c_it->eta(),c_it->phi(),0.,0.);
	  for(size_t isv=0; isv<selVtx.size(); isv++) slimSuperCluster.addBetaStar(0);
	  
	  slimSuperCluster.setPFEnFractions(0.,1.,0.);
	  slimSuperCluster.setPFMultiplicities(0.,c_it->clustersSize(),0.);
	  
	  if( reco2genJet.find(j) != reco2genJet.end())
	    {
	      uint32_t genJetIdx=reco2genJet[j];
	      const reco::GenJet& genjet=genJets->at( genJetIdx );
	      slimSuperCluster.setGenJet(genjet.pt(),genjet.eta(),genjet.phi(),genjet.mass(),genjet.jetArea());
	      
	      if(genJet2Parton.find( genJetIdx )!=genJet2Parton.end())
		{
		  const reco::GenParticle & p = dynamic_cast<const reco::GenParticle &>( (*genParticles)[ genJet2Parton[genJetIdx] ] );
		  slimSuperCluster.setParton(p.pt(),p.eta(),p.phi(),p.pdgId());
		}
	      
	      if(genJet2Stable.find(genJetIdx)!=genJet2Stable.end())
		{
		  const reco::GenParticle & p = dynamic_cast<const reco::GenParticle &>( (*genParticles)[ genJet2Stable[genJetIdx] ] );
		  G4InteractionPositionInfo intInfo=getInteractionPosition(SimTk.product(),SimVtx.product(),genBarcodes->at(genJet2Stable[genJetIdx]));
		  math::XYZVectorD hitPos=intInfo.pos; 
		  slimSuperCluster.setStable(p.pdgId(),
					     p.pt(),      p.eta(),     p.phi(),  
					     hitPos.x(),  hitPos.y(),  hitPos.z());
		}
	    }

	  //iterate over clusters
	  for(reco::CaloCluster_iterator cIt = c_it->clustersBegin(); cIt!=c_it->clustersEnd(); cIt++)
	    {	  
	      const reco::CaloCluster *cl = cIt->get(); 
	      SlimmedCluster slimCluster( cl->energy(),cl->eta(),cl->phi(),(int)cl->hitsAndFractions().size() );
	      slimCluster.roiidx_=slimmedROIs_->size();	  
	      
	      //run pca analysis
	      PCAShowerAnalysis pca;
	      PCAShowerAnalysis::PCASummary_t pcaSummary=pca.computeShowerParameters( *cl, *slimmedRecHits_);
	      slimCluster.center_x_ = pcaSummary.center_x;
	      slimCluster.center_y_ = pcaSummary.center_y;
	      slimCluster.center_z_ = pcaSummary.center_z;
	      slimCluster.axis_x_   = pcaSummary.axis_x;
	      slimCluster.axis_y_   = pcaSummary.axis_y;
	      slimCluster.axis_z_   = pcaSummary.axis_z;
	      slimCluster.ev_1_     = pcaSummary.ev_1;
	      slimCluster.ev_2_     = pcaSummary.ev_2;
	      slimCluster.ev_3_     = pcaSummary.ev_3;
	      slimCluster.sigma_1_  = pcaSummary.sigma_1;
	      slimCluster.sigma_2_  = pcaSummary.sigma_2;
	      slimCluster.sigma_3_  = pcaSummary.sigma_3;
	      
	      GlobalPoint pcaShowerPos(pcaSummary.center_x,pcaSummary.center_y,pcaSummary.center_z);
	      GlobalVector pcaShowerDir(pcaSummary.axis_x,pcaSummary.axis_y,pcaSummary.axis_z);
	      for (unsigned int ih=0;ih<cl->hitsAndFractions().size();++ih) 
		{
		  uint32_t id = (cl->hitsAndFractions())[ih].first.rawId();
		  SlimmedRecHitCollection::iterator theHit=std::find(slimmedRecHits_->begin(),
								     slimmedRecHits_->end(),
								     SlimmedRecHit(id));
		  if(theHit==slimmedRecHits_->end()) continue;
		  
		  theHit->clustId_=slimmedClusters_->size();
		  
		  GlobalPoint cellPos(theHit->x_,theHit->y_,theHit->z_);
		  float cellSize = theHit->cellSize_;
		  float lambda = (cellPos.z()-pcaShowerPos.z())/pcaShowerDir.z();
		  GlobalPoint interceptPos = pcaShowerPos + lambda*pcaShowerDir;
		  float absdx=std::fabs(cellPos.x()-interceptPos.x());
		  float absdy=std::fabs(cellPos.y()-interceptPos.y());
		  
		  theHit->isIn3x3_ = (absdx<cellSize*3./2. && absdy<cellSize*3./2.);
		  theHit->isIn5x5_ = (absdx<cellSize*5./2. && absdy<cellSize*5./2.);
		  theHit->isIn7x7_ = (absdx<cellSize*7./2. && absdy<cellSize*7./2.);
		}
	      
	      slimmedClusters_->push_back(slimCluster);
	    }

	  //all done with this jet
	  slimmedROIs_->push_back(slimSuperCluster);
	}
    }
  else
    {
      for(size_t j=0; j<pfJets->size(); j++)
	{
	  const reco::PFJet &jet=pfJets->at(j);
	  
	  if(jet.pt()<10 || fabs(jet.eta())<1.5 || fabs(jet.eta())>3.0) continue;
	  
	  SlimmedROI slimJet(jet.pt(),jet.eta(),jet.phi(),jet.mass(),jet.jetArea());
	  
	  for(size_t isv=0; isv<selVtx.size(); isv++)
	    {
	      size_t iv=selVtx[isv];
	      std::pair<float,float> beta=betaVariables( &jet, &(vtxH->at(iv)), *vtxH);
	      slimJet.addBetaStar(beta.second);
	    }
	  
	  slimJet.setPFEnFractions(jet.neutralHadronEnergyFraction(),
				   jet.photonEnergyFraction(),
				   jet.chargedHadronEnergyFraction());
	  slimJet.setPFMultiplicities(jet.neutralHadronMultiplicity(),
				      jet.photonMultiplicity(),
				      jet.chargedHadronMultiplicity());
	  
	  if( reco2genJet.find(j) != reco2genJet.end())
	    {
	      uint32_t genJetIdx=reco2genJet[j];
	      
	      const reco::GenJet& genjet=genJets->at( genJetIdx );
	      slimJet.setGenJet(genjet.pt(),genjet.eta(),genjet.phi(),genjet.mass(),genjet.jetArea());
	      
	      if(genJet2Parton.find( genJetIdx )!=genJet2Parton.end())
		{
		  const reco::GenParticle & p = dynamic_cast<const reco::GenParticle &>( (*genParticles)[ genJet2Parton[genJetIdx] ] );
		  slimJet.setParton(p.pt(),p.eta(),p.phi(),p.pdgId());
		  
		}
	      
	      if(genJet2Stable.find(genJetIdx)!=genJet2Stable.end())
		{
		  const reco::GenParticle & p = dynamic_cast<const reco::GenParticle &>( (*genParticles)[ genJet2Stable[genJetIdx] ] );
		  G4InteractionPositionInfo intInfo=getInteractionPosition(SimTk.product(),SimVtx.product(),genBarcodes->at(genJet2Stable[genJetIdx]));
		  math::XYZVectorD hitPos=intInfo.pos; 
		  slimJet.setStable(p.pdgId(),
				    p.pt(),      p.eta(),     p.phi(),  
				    hitPos.x(),  hitPos.y(),  hitPos.z());
		}
	    }
	  	  
	  //first find all pf clusters used
	  std::set<const reco::PFBlockElementCluster *> pfClusters;
	  std::vector<reco::PFCandidatePtr> jetConst(jet.getPFConstituents());
	  for(std::vector<reco::PFCandidatePtr>::iterator cIt=jetConst.begin();
	      cIt!=jetConst.end(); 
	      cIt++)
	    {
	      const reco::PFCandidate::ElementsInBlocks&einb=(*cIt)->elementsInBlocks();
	      for(size_t ieleinb=0; ieleinb<einb.size(); ieleinb++)
		{
		  const reco::PFBlockRef blockRef = einb[ieleinb].first;
		  
		  const edm::OwnVector< reco::PFBlockElement > &eleList=blockRef->elements();
		  for(unsigned int iEle=0; iEle<eleList.size(); iEle++)
		    {
		      reco::PFBlockElement::Type eletype = eleList[iEle].type();
		      if(eletype!=reco::PFBlockElement::HGC_ECAL && eletype!=reco::PFBlockElement::HGC_HCALF && eletype!=reco::PFBlockElement::HGC_HCALB) continue;
		      pfClusters.insert( dynamic_cast<const reco::PFBlockElementCluster*>(&(eleList[iEle])) );
		    }
		}
	    }
	  
	  //iterate of the clusters
	  for(std::set<const reco::PFBlockElementCluster *>::iterator cIt=pfClusters.begin();
	      cIt!=pfClusters.end();
	      cIt++)
	    {
	      const reco::PFClusterRef &cl=(*cIt)->clusterRef();
	      SlimmedCluster slimCluster( cl->energy(),cl->eta(),cl->phi(),(int)cl->hitsAndFractions().size() );
	      slimCluster.roiidx_=slimmedROIs_->size();	  
	      
	      //run pca analysis
	      PCAShowerAnalysis pca;
	      const reco::CaloCluster *caloCl=dynamic_cast<const reco::CaloCluster *>(cl.get());
	      PCAShowerAnalysis::PCASummary_t pcaSummary=pca.computeShowerParameters( *caloCl, *slimmedRecHits_);
	      slimCluster.center_x_ = pcaSummary.center_x;
	      slimCluster.center_y_ = pcaSummary.center_y;
	      slimCluster.center_z_ = pcaSummary.center_z;
	      slimCluster.axis_x_   = pcaSummary.axis_x;
	      slimCluster.axis_y_   = pcaSummary.axis_y;
	      slimCluster.axis_z_   = pcaSummary.axis_z;
	      slimCluster.ev_1_     = pcaSummary.ev_1;
	      slimCluster.ev_2_     = pcaSummary.ev_2;
	      slimCluster.ev_3_     = pcaSummary.ev_3;
	      slimCluster.sigma_1_  = pcaSummary.sigma_1;
	      slimCluster.sigma_2_  = pcaSummary.sigma_2;
	      slimCluster.sigma_3_  = pcaSummary.sigma_3;
	      
	      GlobalPoint pcaShowerPos(pcaSummary.center_x,pcaSummary.center_y,pcaSummary.center_z);
	      GlobalVector pcaShowerDir(pcaSummary.axis_x,pcaSummary.axis_y,pcaSummary.axis_z);
	      for (unsigned int ih=0;ih<cl->hitsAndFractions().size();++ih) 
		{
		  uint32_t id = (cl->hitsAndFractions())[ih].first.rawId();
		  SlimmedRecHitCollection::iterator theHit=std::find(slimmedRecHits_->begin(),
								     slimmedRecHits_->end(),
								     SlimmedRecHit(id));
		  if(theHit==slimmedRecHits_->end()) continue;
		  
		  theHit->clustId_=slimmedClusters_->size();
		  
		  GlobalPoint cellPos(theHit->x_,theHit->y_,theHit->z_);
		  float cellSize = theHit->cellSize_;
		  float lambda = (cellPos.z()-pcaShowerPos.z())/pcaShowerDir.z();
		  GlobalPoint interceptPos = pcaShowerPos + lambda*pcaShowerDir;
		  float absdx=std::fabs(cellPos.x()-interceptPos.x());
		  float absdy=std::fabs(cellPos.y()-interceptPos.y());
		  
		  theHit->isIn3x3_ = (absdx<cellSize*3./2. && absdy<cellSize*3./2.);
		  theHit->isIn5x5_ = (absdx<cellSize*5./2. && absdy<cellSize*5./2.);
		  theHit->isIn7x7_ = (absdx<cellSize*7./2. && absdy<cellSize*7./2.);
		}
	      
	      slimmedClusters_->push_back(slimCluster);
	    }

	  //all done with this jet
	  slimmedROIs_->push_back(slimJet);
	}
    }
  
  //remove unclustered vertices
  std::cout << slimmedRecHits_->size() << "->";
  slimmedRecHits_->erase(remove_if(slimmedRecHits_->begin(), slimmedRecHits_->end(), SlimmedRecHit::IsNotClustered),
			 slimmedRecHits_->end());
  std::cout << slimmedRecHits_->size() << std::endl;

  //all done, fill tree
  if(slimmedROIs_->size())  tree_->Fill();
}

//
void HGCROIAnalyzer::endJob() 
{ 
}


//define this as a plug-in
DEFINE_FWK_MODULE(HGCROIAnalyzer);
