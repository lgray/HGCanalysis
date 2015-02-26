#include "UserCode/HGCanalysis/plugins/HGCJetAnalyzer.h"

#include "DataFormats/Math/interface/deltaR.h"

#include <iostream>

using namespace std;

//
HGCJetAnalyzer::HGCJetAnalyzer( const edm::ParameterSet &iConfig )
{
  //configure analyzer
  genJetsSource_    = iConfig.getUntrackedParameter<std::string>("genJetsSource");
  pfJetsSource_     = iConfig.getUntrackedParameter< std::string >("pfJetsSource");

  edm::Service<TFileService> fs;
  histMap_["kin_gen"]         = fs->make<TH2F>("kin_gen",         ";Transverse momentum [GeV];Pseudo-rapidity;Jets",    40, 0, 1000, 15, 1.5, 3.0);
  histMap_["kin_gen_matched"] = fs->make<TH2F>("kin_gen_matched", ";Transverse momentum [GeV];Pseudo-rapidity;Jets",    40, 0, 1000, 15, 1.5, 3.0);
  histMap_["kin_rec"]         = fs->make<TH2F>("kin_rec",         ";Transverse momentum [GeV];Pseudo-rapidity;Jets",    40, 0, 1000, 15, 1.5, 3.0);
  histMap_["kin_rec_unm"]     = fs->make<TH2F>("kin_rec_unm",     ";Transverse momentum [GeV];Pseudo-rapidity;Jets",    40, 0, 1000, 15, 1.5, 3.0);
  histMap_["ptresp_pt"]       = fs->make<TH2F>("ptresp_pt",       ";p_{T}(gen) [GeV];p_{T}(reco)/p_{T}(gen);Jets",      40, 0,   1000, 100,0,4);
  histMap_["ptresp_eta"]      = fs->make<TH2F>("ptresp_eta",      ";Pseudo-rapidity (gen);p_{T}(reco)/p_{T}(gen);Jets", 15, 1.5, 3.0,  100,0,4);
  TString comp[]={"gamma","chf","nhf"};
  TString jetType[]={"","_unm"};
  for(size_t i=0; i<sizeof(comp)/sizeof(TString); i++)
    {
      for(size_t j=0; j<sizeof(jetType)/sizeof(TString); j++)
	{
	  histMap_[comp[i]+"_pt_enfrac"+jetType[j] ] = fs->make<TH2F>(comp[i]+"_pt_enfrac"+jetType[j],     ";Transverse momentum [GeV];"+comp[i]+" fraction;Jets",             40, 0, 1000, 50,0,1);
	  histMap_[comp[i]+"_pt_en"+jetType[j] ] = fs->make<TH2F>(comp[i]+"_pt_en"+jetType[j],         ";Transverse momentum [GeV];"+comp[i]+" total energy [GeV];Jets",         40, 0, 1000, 50,0,500);
	  histMap_[comp[i]+"_pt_inden"+jetType[j]] = fs->make<TH2F>(comp[i]+"_pt_inden"+jetType[j],      ";Transverse momentum [GeV];"+comp[i]+" energy [GeV];Jets", 40, 0, 1000, 50,0,250);
	  histMap_[comp[i]+"_pt_mult"+jetType[j]] = fs->make<TH2F>(comp[i]+"_pt_mult"+jetType[j],       ";Transverse momentum [GeV];"+comp[i]+" multiplicity;Jets",         40, 0, 1000, 80,0,80);
	  histMap_[comp[i]+"_eta_enfrac"+jetType[j]] = fs->make<TH2F>(comp[i]+"_eta_enfrac"+jetType[j],     ";Pseudo-rapidity;"+comp[i]+" fraction;Jets",             15, 1.5, 3.0, 50,0,1);
	  histMap_[comp[i]+"_eta_en"+jetType[j]] = fs->make<TH2F>(comp[i]+"_eta_en"+jetType[j],         ";Pseudo-rapidity;"+comp[i]+" total energy [GeV];Jets",         15, 1.5, 3.0, 50,0,250);
	  histMap_[comp[i]+"_eta_inden"+jetType[j]] = fs->make<TH2F>(comp[i]+"_eta_inden"+jetType[j],      ";Pseudo-rapidity;"+comp[i]+" energy [GeV];Jets", 15, 1.5, 3.0, 50,0,250);
	  histMap_[comp[i]+"_eta_mult"+jetType[j]] = fs->make<TH2F>(comp[i]+"_eta_mult"+jetType[j],       ";Pseudo-rapidity;"+comp[i]+" multiplicity;Jets",         15, 1.5, 3.0, 50,0,50);
	}
    }
  
 
  for(std::map<TString,TH2F *>::iterator it=histMap_.begin(); it!=histMap_.end(); it++) it->second->Sumw2();
}

//
HGCJetAnalyzer::~HGCJetAnalyzer()
{
}

//
void HGCJetAnalyzer::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  //generator level jets
  edm::Handle<reco::GenJetCollection> genJets;
  iEvent.getByLabel(edm::InputTag(genJetsSource_), genJets);

  for(size_t j=0; j<genJets->size(); j++)
    {
      const reco::GenJet& genjet=genJets->at(j);
      if(fabs(genjet.eta())<1.7 || fabs(genjet.eta())>2.8) continue;
      histMap_["kin_gen"]->Fill(genjet.pt(),fabs(genjet.eta()));
    }


  //PF jets
  edm::Handle<std::vector<reco::PFJet> > pfJets;
  iEvent.getByLabel(edm::InputTag(pfJetsSource_),pfJets);

  //loop over jets
  std::map<TString, std::vector<float> > compInfo;
  compInfo["gamma"] = std::vector<float>(3,0);
  compInfo["chf"]   = std::vector<float>(3,0);
  compInfo["nhf"]   = std::vector<float>(3,0);
  for(size_t i=0; i<pfJets->size(); i++)
    {
      const reco::PFJet &jet=pfJets->at(i);

      float pt  = jet.pt();
      float eta = fabs(jet.eta());
      if(eta<1.7 || eta>2.8) continue;

      compInfo["nhf"][0]   = jet.neutralHadronEnergyFraction();
      compInfo["nhf"][1]   = jet.neutralHadronEnergy();
      compInfo["nhf"][2]   = jet.neutralHadronMultiplicity();
      compInfo["gamma"][0] = jet.photonEnergyFraction();
      compInfo["gamma"][1] = jet.photonEnergy();
      compInfo["gamma"][2] = jet.photonMultiplicity();
      compInfo["chf"][0]   = jet.chargedHadronEnergyFraction();
      compInfo["chf"][1]   = jet.chargedHadronEnergy();
      compInfo["chf"][2]   = jet.chargedHadronMultiplicity();

      //match by DR
      float genpt(0),geneta(0),minDR(0.4),minDeltaPt(99999.);
      for(size_t j=0; j<genJets->size(); j++)
	{
	  const reco::GenJet& genjet=genJets->at(j);
	  if(fabs(genjet.eta())<1.7 || fabs(genjet.eta())>2.8) continue;

	  float dR=deltaR(genjet,jet);
	  if(dR>minDR) continue;
	  float dPt=fabs(pt-genjet.pt());
	  if(dPt>minDeltaPt) continue;
	  minDR=dR;
	  minDeltaPt=dPt;
	  genpt=genjet.pt();
	  geneta=fabs(genjet.eta());
	}

      //fill histos
      TString jetType("");
      if(genpt>0) 
	{
	  histMap_["kin_gen_matched"]->Fill(genpt,geneta);
	  histMap_["kin_rec"]->Fill(pt,eta);
	  histMap_["ptresp_pt"]->Fill(genpt,pt/genpt);
	  histMap_["ptresp_eta"]->Fill(geneta,pt/genpt);
	}
      else
	{
	  jetType="_unm";
	  histMap_["kin_rec_unm"]->Fill(pt,eta);
	}

      //jet components
      float ptToUse(genpt>0 ? genpt : pt);
      float etaToUse(genpt>0 ? geneta : eta);
      for(std::map<TString, std::vector<float> >::iterator it=compInfo.begin();
	  it!=compInfo.end();
	  it++)
	{
	  if(pt<20) continue;

	  histMap_[ it->first + "_pt_enfrac"+jetType ] ->Fill( ptToUse, it->second[0] );
	  histMap_[ it->first + "_pt_en"+jetType ]     ->Fill( ptToUse, it->second[1] );
	  histMap_[ it->first + "_pt_mult"+jetType ]   ->Fill( ptToUse, it->second[2] );

	  histMap_[ it->first + "_eta_enfrac"+jetType ] ->Fill( etaToUse, it->second[0] );
	  histMap_[ it->first + "_eta_en"+jetType ]     ->Fill( etaToUse, it->second[1] );
	  histMap_[ it->first + "_eta_mult"+jetType ]   ->Fill( etaToUse, it->second[2] );
	}
      
      //loop over constituents
      std::vector<reco::PFCandidatePtr> jetConst(jet.getPFConstituents());
      for(std::vector<reco::PFCandidatePtr>::iterator cIt=jetConst.begin();cIt!=jetConst.end(); cIt++)
	{
	  TString compType("");
	  if( (*cIt)->particleId() == reco::PFCandidate::gamma )     compType="gamma";
	  else if( (*cIt)->particleId() == reco::PFCandidate::h )    compType="chf";
	  else if( (*cIt)->particleId() == reco::PFCandidate::h0 )   compType="nhf";
	  else continue;
	  histMap_[ compType + "_pt_inden"  + jetType ]->Fill( ptToUse, (*cIt)->energy() );
	  histMap_[ compType + "_eta_inden" + jetType ]->Fill( etaToUse, (*cIt)->energy() );
 	}
    }
}

//
void HGCJetAnalyzer::endJob() { }


//define this as a plug-in
DEFINE_FWK_MODULE(HGCJetAnalyzer);
