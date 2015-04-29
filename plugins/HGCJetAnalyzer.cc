#include "UserCode/HGCanalysis/plugins/HGCJetAnalyzer.h"
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
#include <iostream>
#include <unordered_map>

using namespace std;

//
HGCJetAnalyzer::HGCJetAnalyzer( const edm::ParameterSet &iConfig )
{
  //configure analyzer
  genSource_        = iConfig.getUntrackedParameter< std::string >("genSource");
  genJetsSource_    = iConfig.getUntrackedParameter<std::string>("genJetsSource");
  pfJetsSource_     = iConfig.getUntrackedParameter< std::string >("pfJetsSource");
  eeRecHitsSource_  = iConfig.getUntrackedParameter< std::string >("eeRecHitsSource");
  hefRecHitsSource_  = iConfig.getUntrackedParameter< std::string >("hefRecHitsSource");

  edm::Service<TFileService> fs;

  const Double_t PTBINS[]={0,5,10,20,40,60,80,100,120,150,200,250,500,1000};
  const Int_t NPTBINS=sizeof(PTBINS)/sizeof(Double_t)-1;

  histMap_["gen"]           = fs->make<TH2F>("gen",        ";Transverse momentum [GeV];Pseudo-rapidity;Jets/bin width",              NPTBINS,PTBINS, 15, 1.5, 3.0);
  histMap_["gen_parton"]    = fs->make<TH2F>("gen_parton", ";Transverse momentum [GeV];Pseudo-rapidity;Parton-matching efficiency",  NPTBINS,PTBINS, 15, 1.5, 3.0);
  histMap_["gen_reco"]      = fs->make<TH2F>("gen_reco",   ";Transverse momentum [GeV];Pseudo-rapidity;Reco-matching efficiency",    NPTBINS,PTBINS, 15, 1.5, 3.0);
  for(std::map<TString,TH2F *>::iterator it=histMap_.begin(); it!=histMap_.end(); it++) it->second->Sumw2();

  jetTree_=fs->make<TTree>("HGCJets","HGCJets");
  jetTree_->Branch("jpt",&jpt_,"jpt/F");
  jetTree_->Branch("jeta",&jeta_,"jeta/F");
  jetTree_->Branch("jphi",&jphi_,"jphi/F");
  jetTree_->Branch("jmass",&jmass_,"jmass/F");
  jetTree_->Branch("jnhf",&jnhf_,"jnhf/F");
  jetTree_->Branch("jnhe",&jnhe_,"jnhe/F");
  jetTree_->Branch("jnhm",&jnhm_,"jnhm_/F");
  jetTree_->Branch("jgf",&jgf_,"jgf/F");
  jetTree_->Branch("jge",&jge_,"jge/F");
  jetTree_->Branch("jgm",&jgm_,"jgm/F");
  jetTree_->Branch("jchf",&jchf_,"jchf/F");
  jetTree_->Branch("jche",&jche_,"jche/F");
  jetTree_->Branch("jchm",&jchm_,"jchm/F");
  jetTree_->Branch("gjpt",&gjpt_,"gjpt/F");
  jetTree_->Branch("gjeta",&gjeta_,"gjeta/F");
  jetTree_->Branch("gjphi",&gjphi_,"gjphi/F");
  jetTree_->Branch("gjmass",&gjmass_,"gjmass/F");
  jetTree_->Branch("ppt",&ppt_,"ppt/F");
  jetTree_->Branch("peta",&peta_,"peta/F");
  jetTree_->Branch("pphi",&pphi_,"pphi/F");
  jetTree_->Branch("pmass",&pmass_,"pmass/F");
  jetTree_->Branch("pid",&pid_,"pid/F");
  jetTree_->Branch("nTDCHits",&nTDCHits_,"nTDCHits/I");
}

//
HGCJetAnalyzer::~HGCJetAnalyzer()
{
}

//
void HGCJetAnalyzer::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
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
  edm::Handle<reco::GenJetCollection> genJets;
  iEvent.getByLabel(edm::InputTag(genJetsSource_), genJets);
  edm::Handle<edm::View<reco::Candidate> > genParticles;
  iEvent.getByLabel(edm::InputTag(genSource_), genParticles);
  edm::Handle<std::vector<reco::PFJet> > pfJets;
  iEvent.getByLabel(edm::InputTag(pfJetsSource_),pfJets);
  std::unordered_map<uint32_t,uint32_t> reco2genJet,genJet2Parton;
  for(size_t j=0; j<genJets->size(); j++)
    {
      const reco::GenJet& genjet=genJets->at(j);
      float pt=genjet.pt();
      float abseta=fabs(genjet.eta());
      if(abseta<1.5 || abseta>3.0) continue;
     
      //gen particle matching
      bool genMatched(false);
      float minDPt(99999.);
      for(size_t i = 0; i < genParticles->size(); ++ i)
	{
	  const reco::GenParticle & p = dynamic_cast<const reco::GenParticle &>( (*genParticles)[i] );
	  if(p.status()!=2) continue;
	  if( !(abs(p.pdgId())==21 || abs(p.pdgId())<6) ) continue;
	  float dR(deltaR(p,genjet));
	  if(dR>0.4) continue;
	  float dPt( fabs(p.pt()-pt));
	  if(dPt>minDPt) continue;
	  minDPt=dPt;
	  genJet2Parton[j]=i;
	  genMatched=true;
	}
      
      //reco matching
      bool recoMatched(false);
      float minDR=0.2;
      for(size_t i=0; i<pfJets->size(); i++)
	{
	  const reco::PFJet &jet=pfJets->at(i);
	  float dR=deltaR(jet,genjet);
	  if(dR>minDR) continue;
	  minDR=dR;
	  if(reco2genJet.find(i)!=reco2genJet.end()) continue;
	  reco2genJet[i]=j;
	  recoMatched=true;
	}

      //kinematics and matching efficiencies
      Int_t ptbin=histMap_["gen"]->GetXaxis()->FindBin(pt);
      Float_t binwidth=histMap_["gen"]->GetXaxis()->GetBinWidth(ptbin);
      histMap_["gen"]->Fill(pt,abseta,1./binwidth);
      if(genMatched)  histMap_["gen_parton"]->Fill(pt,abseta,1./binwidth);
      if(recoMatched) histMap_["gen_reco"]->Fill(pt,abseta,1./binwidth);
    }


  //
  //map RECHits by DetId (to be used later)
  //
  std::unordered_map<uint32_t,uint32_t> eeRecHitsIdMap;
  edm::Handle<HGCRecHitCollection> eeRecHits;
  iEvent.getByLabel(edm::InputTag("HGCalRecHit",eeRecHitsSource_),eeRecHits); 
  if(eeRecHits.isValid())
    {
      uint32_t recHitCtr=0;
      for(HGCRecHitCollection::const_iterator hit_it=eeRecHits->begin(); hit_it!=eeRecHits->end(); hit_it++,recHitCtr++)
	{
	  uint32_t recoDetId(hit_it->id());
	  eeRecHitsIdMap[recoDetId]=recHitCtr;
	}
    }
  std::unordered_map<uint32_t,uint32_t> hefRecHitsIdMap;
  edm::Handle<HGCRecHitCollection> hefRecHits;
  iEvent.getByLabel(edm::InputTag("HGCalRecHit",hefRecHitsSource_),hefRecHits); 
  if(hefRecHits.isValid())
    {
      uint32_t recHitCtr=0;
      for(HGCRecHitCollection::const_iterator hit_it=hefRecHits->begin(); hit_it!=hefRecHits->end(); hit_it++,recHitCtr++)
	{
	  uint32_t recoDetId(hit_it->id());
	  hefRecHitsIdMap[recoDetId]=recHitCtr;
	}
    }
  
  //
  // Analyze matched jets
  // 
  for(std::unordered_map<uint32_t,uint32_t>::iterator recoIt=reco2genJet.begin();
      recoIt!=reco2genJet.end();
      recoIt++)
    {
      const reco::PFJet &jet=pfJets->at(recoIt->first);
      jpt_  = jet.pt();
      jeta_ = jet.eta();
      jphi_ = jet.phi();
      jmass_= jet.mass();
      jnhf_ = jet.neutralHadronEnergyFraction();
      jnhe_ = jet.neutralHadronEnergy();
      jnhm_ = jet.neutralHadronMultiplicity();
      jgf_  = jet.photonEnergyFraction();
      jge_  = jet.photonEnergy();
      jgm_  = jet.photonMultiplicity();
      jchf_ = jet.chargedHadronEnergyFraction();
      jche_ = jet.chargedHadronEnergy();
      jchm_ = jet.chargedHadronMultiplicity();

      const reco::GenJet& genjet=genJets->at(recoIt->second);
      gjpt_ = genjet.pt();
      gjeta_ = genjet.eta();
      gjphi_ = genjet.phi();
      gjmass_ = genjet.mass();

      ppt_=0; peta_=0; pphi_=0; pmass_=0; pid_=0;
      if(genJet2Parton.find(recoIt->second)!=genJet2Parton.end())
	{
	  const reco::GenParticle & p = dynamic_cast<const reco::GenParticle &>( (*genParticles)[ genJet2Parton[recoIt->second] ] );
	  ppt_   = p.pt();
	  peta_  = p.eta();
	  pphi_  = p.phi();
	  pmass_ = p.mass();
	  pid_   = p.pdgId();
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
	  
      //analyze rec hits which have been clustered
      nTDCHits_=0;
      for(std::set<const reco::PFBlockElementCluster *>::iterator sc=pfClusters.begin();
	  sc!=pfClusters.end();
	  sc++)
	{
	  for( const auto& rhf : (*sc)->clusterRef()->hitsAndFractions() )
	    {
	      uint32_t recoDetId( rhf.first.rawId() );
	      
	      const HGCRecHit *hitPtr=0;
	      std::unordered_map<uint32_t,uint32_t>::iterator ptoHit=eeRecHitsIdMap.find(recoDetId);
	      if(ptoHit==eeRecHitsIdMap.end()) {
		ptoHit=hefRecHitsIdMap.find(recoDetId);
		if(ptoHit==hefRecHitsIdMap.end()) continue;
		hitPtr = & ( (*(hefRecHits))[ ptoHit->second] );
	      }
	      else{
		hitPtr = & ( (*(eeRecHits))[ ptoHit->second] );
	      }
	      //uint32_t layIdx(((recoDetId>>19)&0x1f));
	      float toa=hitPtr->time();
	      //float frac=rhf.second;
	      //float en=hitPtr->energy()*1.0e6;		
	      
	      nTDCHits_ += (toa>0);
	      //if(toa>0)
	      //	std::cout << "Jet #"<< recoIt->first
	      //		  << " eta=" << jeta_
	      //		  << " pt/GeV=" << jpt_
	      //		  << " gen pt/GeV=" << gjpt_ 
	      //		  << " DetId 0x" << hex << recoDetId << dec
	      //		  << " layer #" << layIdx << dec
	      //		  << " E/keV=" << en
	      //		  << " frac=" << frac
	      //		  << " t/ns=" << toa 
	      //		  << std::endl;
	    }
	}
      cout << nTDCHits_ << " " << recoIt->first << " " << recoIt->second 
	   << " " << gjpt_ << " " << gjeta_ << " " << gjphi_ <<  endl;

      //fill tree with information
      jetTree_->Fill();
    }
  std::cout << "-----" << std::endl;

}

//
void HGCJetAnalyzer::endJob() 
{ 
  histMap_["gen_parton"]->Divide(histMap_["gen"]);
  histMap_["gen_reco"]->Divide(histMap_["gen"]);
}


//define this as a plug-in
DEFINE_FWK_MODULE(HGCJetAnalyzer);
