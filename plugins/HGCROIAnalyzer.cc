#include "UserCode/HGCanalysis/plugins/HGCROIAnalyzer.h"

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
#include "TMath.h"
#include "Math/DistFunc.h"

#include <iostream>

using namespace std;

//
HGCROIAnalyzer::HGCROIAnalyzer( const edm::ParameterSet &iConfig )
{
  //configure analyzer
  genSource_            = iConfig.getUntrackedParameter<std::string>("genSource");
  genJetsSource_        = iConfig.getUntrackedParameter<std::string>("genJetsSource");
  geometrySource_       = iConfig.getUntrackedParameter< std::vector<std::string> >("geometrySource");
  hitCollections_       = iConfig.getUntrackedParameter< std::vector<std::string> >("hitCollections");
  mipEn_                = iConfig.getUntrackedParameter< std::vector<double> >("mipEn");
  taggingMode_          = iConfig.getUntrackedParameter<bool>("taggingMode");
  saveHitTree_          = iConfig.getUntrackedParameter<bool>("saveHitTree");
  roipuParamFile_       = iConfig.getUntrackedParameter<edm::FileInPath>("roipuParamFile");
  

  edm::Service<TFileService> fs;

  //init regions of interest
  nLayerBins_=54;
  nEtaBins_=15;
  regsH_=fs->make<TH2F>("regs",";Layer;Pseudo-rapidity;",nLayerBins_,1,nLayerBins_+1,nEtaBins_,1.5,3.0);
  for(Int_t xbin=1; xbin<=regsH_->GetYaxis()->GetNbins(); xbin++)
    {
      roi_.add(regsH_->GetYaxis()->GetBinCenter(xbin),0.);
      roi_.rotateInPhi(xbin-1);
    }

  //energy flux distribution wrt to ROI centers
  ndRbins_=12;   drMin_=0;    drMax_=0.6;
  nCsiBins_=56;  csiMin_=-3;  csiMax_=10;
  csidrH_=fs->make<TH2F>("csidr",";#Delta R;log #xi;",
			 ndRbins_*nLayerBins_, drMin_,  drMax_+(drMax_-drMin_)*(nLayerBins_-1),
			 nCsiBins_*nEtaBins_,  csiMin_, csiMax_+(csiMax_-csiMin_)*(nEtaBins_-1));
  csitdrH_=fs->make<TH2F>("csitdr",";#Delta R;log #xi_{T};",
			  ndRbins_*nLayerBins_, drMin_,  drMax_+(drMax_-drMin_)*(nLayerBins_-1),
			  nCsiBins_*nEtaBins_,  csiMin_, csiMax_+(csiMax_-csiMin_)*(nEtaBins_-1));
  
  //read cell counting from file
  TFile *fIn=TFile::Open(roipuParamFile_.fullPath().c_str());
  medianPU_csiH_  = (TH2F *)fIn->Get("medianPU_csi");  medianPU_csiH_->SetDirectory(0);
  widthPU_csiH_   = (TH2F *)fIn->Get("widthPU_csi");   widthPU_csiH_->SetDirectory(0);
  medianPU_csitH_ = (TH2F *)fIn->Get("medianPU_csit");  medianPU_csitH_->SetDirectory(0);
  widthPU_csitH_  = (TH2F *)fIn->Get("widthPU_csit");   widthPU_csitH_->SetDirectory(0);
  fIn->Close();

  //init tree only if not in tagging mode
  if(!taggingMode_)
    {
      if(saveHitTree_)
	{
	  hitT_=fs->make<TTree>("HGC","Event Summary");
	  initHGCSimulationEventTree(hitT_,simEvt_);
	}
      roiT_=fs->make<TTree>("HGCROI","ROI Summary");
      initHGCROITree(roiT_, roiEvt_);
    }
  
  //start counter
  evtCtr_=0;
}

//
HGCROIAnalyzer::~HGCROIAnalyzer()
{
}

//
void HGCROIAnalyzer::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  evtCtr_++;

  //event header
  simEvt_.run    = iEvent.id().run();          roiEvt_.run   =iEvent.id().run();
  simEvt_.lumi   = iEvent.luminosityBlock();   roiEvt_.lumi  = iEvent.luminosityBlock(); 
  simEvt_.event  = evtCtr_;                    roiEvt_.event = iEvent.id().event();

  //generator level particles
  simEvt_.ngen=0;
  edm::Handle<edm::View<reco::Candidate> > genParticles;
  iEvent.getByLabel(edm::InputTag(genSource_), genParticles);
  if(genParticles.isValid())
    {
      for(size_t i = 0; i < genParticles->size(); ++ i)
	{
	  const reco::GenParticle & p = dynamic_cast<const reco::GenParticle &>( (*genParticles)[i] );
	  
	  //hard process
	  bool isOutgoingNeutrino( p.status()==1 && (abs(p.pdgId())==12||abs(p.pdgId())==14||abs(p.pdgId())==16) );
	  //bool isStatus2qg( p.status()==2 && (abs(p.pdgId())<6 || abs(p.pdgId())==21) );
	  bool isStatus3q( p.status()==3 && abs(p.pdgId())<6);
	  bool isStatus2l( p.status()==2 && (abs(p.pdgId())==11 || abs(p.pdgId())==13 || abs(p.pdgId())==14));
	  //if(!isOutgoingNeutrino && !isStatus2qg && !isStatus2l) continue;
	  if(!isOutgoingNeutrino && !isStatus3q && !isStatus2l) continue;
	  
	  simEvt_.gen_id[simEvt_.ngen]=p.pdgId();
	  simEvt_.gen_pt[simEvt_.ngen]=p.pt();
	  simEvt_.gen_eta[simEvt_.ngen]=p.eta();
	  simEvt_.gen_phi[simEvt_.ngen]=p.phi();
	  simEvt_.gen_en[simEvt_.ngen]=p.energy();
	  simEvt_.ngen++;

	  if(simEvt_.ngen>=MAXGENPEREVENT) break;
	}
    }

  //reset or shuffle RoI
  if(!taggingMode_) { roi_.reset(); }
  else              { for(Int_t xbin=1; xbin<=regsH_->GetYaxis()->GetNbins(); xbin++) roi_.rotateInPhi(xbin-1); }

  //analyze gen jets
  simEvt_.njgen=0;
  edm::Handle<reco::GenJetCollection> genJets;
  iEvent.getByLabel(edm::InputTag(genJetsSource_), genJets);
  std::vector< std::pair<const reco::GenJet *,int> > matchedGenJets;
  std::vector< int > tagJetIdx;
  if(genJets.isValid() && !taggingMode_)
    {
      //pre-select jets
      for(size_t k=0; k<genJets->size(); ++k)
	{
	  const reco::GenJet & j =(*genJets)[k];
	  if(j.pt()<15 || fabs(j.eta())>4.7) continue;
	  
	  //require not to be matched to lepton (charged or neutrino)
	  int matchedId(0),matchedIdx(-1);
	  for(int igen=0; igen<simEvt_.ngen; igen++)
	    {
	      float dR=deltaR(simEvt_.gen_eta[igen],simEvt_.gen_phi[igen],j.eta(),j.phi());
	      if(dR>0.4) continue;
	      if(matchedIdx>=0)
		if(fabs(simEvt_.gen_en[igen]-j.energy())>fabs(simEvt_.gen_en[matchedIdx]-j.energy())) continue;
	      matchedId=simEvt_.gen_id[igen];
	      matchedIdx=igen;
	    }
	  if(abs(matchedId)>=11 && abs(matchedId)<=16) continue;

	  //store information
	  matchedGenJets.push_back( std::pair<const reco::GenJet *,int>(&j,matchedId) );
	  if(simEvt_.njgen<MAXGENPEREVENT)
	    {
	      simEvt_.genj_en[simEvt_.njgen]=j.energy();
	      simEvt_.genj_eta[simEvt_.njgen]=j.eta();
	      simEvt_.genj_phi[simEvt_.njgen]=j.phi();
	      simEvt_.genj_pt[simEvt_.njgen]=j.pt();
	      simEvt_.genj_emfrac[simEvt_.njgen]=j.emEnergy()/j.energy();
	      simEvt_.genj_hadfrac[simEvt_.njgen]=j.hadEnergy()/j.energy();
	      simEvt_.genj_invfrac[simEvt_.njgen]=j.invisibleEnergy()/j.energy();
	      simEvt_.njgen++;
	    }
	}

      //now do a VBF pre-selection mode
      if(matchedGenJets.size()>=2)
	{
	  float dijetMass=(matchedGenJets[0].first->p4()+matchedGenJets[1].first->p4()).mass();
	  float deta=fabs(matchedGenJets[0].first->eta()-matchedGenJets[1].first->eta());
	  int idx1(0), idx2(1);
	  for(size_t ij=0; ij<matchedGenJets.size(); ij++)
	    for(size_t kj=ij+1; kj<matchedGenJets.size(); kj++)
	      {
		float newDijetMass=(matchedGenJets[ij].first->p4()+matchedGenJets[kj].first->p4()).mass();
		float newDeta=fabs(matchedGenJets[ij].first->eta()-matchedGenJets[kj].first->eta());
		if(newDijetMass<dijetMass) continue;
		dijetMass=newDijetMass;
		deta=newDeta;
		idx1=ij; idx2=kj;
	      }

	  //save if good
	  if(dijetMass>250 && deta>2.5)
	    {
	      if(fabs(matchedGenJets[idx1].first->eta())>1.5 && fabs(matchedGenJets[idx1].first->eta())<3.0)
		{
		  roi_.add(matchedGenJets[idx1].first->eta(),matchedGenJets[idx1].first->phi());
		  tagJetIdx.push_back(idx1);
		}
	      if(fabs(matchedGenJets[idx2].first->eta())>1.5 && fabs(matchedGenJets[idx2].first->eta())<3.0)
		{
		  roi_.add(matchedGenJets[idx2].first->eta(),matchedGenJets[idx2].first->phi());
		  tagJetIdx.push_back(idx2);
		}
	    }
	}
    }
  
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
	  
	  //check if hits are in ROI
	  int roiIdx=roi_.findROI(refPos.eta(),refPos.phi());	      
	  if(roiIdx<0) continue;

	  //compute weight to assign based on layer and distance to center
	  Int_t layerBin     = regsH_->GetXaxis()->FindBin(layer);
	  Int_t etaBin       = regsH_->GetYaxis()->FindBin(TMath::Abs(roi_.roiInfo_[roiIdx].center_eta));
	  
	  float dR=deltaR(roi_.roiInfo_[roiIdx].center_eta,roi_.roiInfo_[roiIdx].center_phi,refPos.eta(),refPos.phi());
	  float cell_csit = dR>0.01 ? TMath::Log( hitEn/TMath::CosH(fabs(refPos.eta())) ) : -9999;
	  float cell_csi  = dR>0.01 ? TMath::Log( hitEn ) : -9999;
	  if(dR>0.01)
	    {
	      if(cell_csi>csiMax_)  cell_csi=csiMax_;
	      if(cell_csi<csiMin_)  cell_csi=csiMin_;
	      if(cell_csit>csiMax_) cell_csit=csiMax_;
	      if(cell_csit<csiMin_) cell_csit=csiMin_;

	      regsH_->Fill( layer, TMath::Abs(roi_.roiInfo_[roiIdx].center_eta) );

	      float dR_ext       = dR+(layerBin-1)*(drMax_-drMin_);

	      float cell_csi_ext = cell_csi+(etaBin-1)*(csiMax_-csiMin_);
	      csidrH_->Fill(dR_ext,cell_csi_ext);
	      
	      float cell_csit_ext = cell_csit+(etaBin-1)*(csiMax_-csiMin_);
	      csitdrH_->Fill(dR_ext,cell_csit_ext);
	    }

	  //if not characterizing, save hits to tree
	  if(taggingMode_) continue;
	  
	  float hitWeight(1.0);
	  if(medianPU_csiH_ && widthPU_csiH_)
	    {
	      Int_t dRbin=medianPU_csiH_->GetYaxis()->FindBin(dR);
	      Int_t weightsYbin( dRbin + (etaBin-1)*ndRbins_ );

	      float medianVal=medianPU_csiH_->GetBinContent(layerBin,weightsYbin);
	      float sigmaVal=widthPU_csiH_->GetBinContent(layerBin,weightsYbin);
	      float chi2=sigmaVal>0 ? (cell_csi>medianVal)*TMath::Power((cell_csi-medianVal)/sigmaVal,2) : 0;
	      hitWeight=chi2 > 0 ? ROOT::Math::chisquared_cdf(chi2,1) : 0;
	    }
	  
	  float hitWeightT(1.0);
	  if(medianPU_csitH_ && widthPU_csitH_)
	    {
	      Int_t dRbin=medianPU_csitH_->GetYaxis()->FindBin(dR);
	      Int_t weightsYbin( dRbin + (etaBin-1)*ndRbins_ );

	      float medianVal=medianPU_csitH_->GetBinContent(layerBin,weightsYbin);
	      float sigmaVal=widthPU_csitH_->GetBinContent(layerBin,weightsYbin);
	      float chi2=sigmaVal>0 ? (cell_csi>medianVal)*TMath::Power((cell_csi-medianVal)/sigmaVal,2) : 0;
	      hitWeightT=chi2 > 0 ? ROOT::Math::chisquared_cdf(chi2,1) : 0;
	    }

	  //update for hit tree
	  simEvt_.hit_layer[simEvt_.nhits]   = layer;
	  simEvt_.hit_x[simEvt_.nhits]       = refPos.x();
	  simEvt_.hit_y[simEvt_.nhits]       = refPos.y();
	  simEvt_.hit_z[simEvt_.nhits]       = refPos.z();
	  simEvt_.hit_eta[simEvt_.nhits]     = refPos.eta();
	  simEvt_.hit_phi[simEvt_.nhits]     = refPos.phi();
	  simEvt_.hit_edep[simEvt_.nhits]    = hitEn;
	  simEvt_.hit_wgt[simEvt_.nhits]     = hitWeight;
	  simEvt_.hit_wgt_t[simEvt_.nhits]   = hitWeightT;
	  simEvt_.nhits++;

	  //update ROI info
	  roi_.roiInfo_[roiIdx].add(hitEn,refPos.x(), refPos.y(), refPos.z(), refPos.eta(), refPos.phi(), hitWeight,layer);
	}
      
      const HGCalTopology &topo=geom->topology();
      const HGCalDDDConstants &dddConst=topo.dddConstants();
      layerCtrOffset +=dddConst.layers(true);
    }
  
  //fill trees
  if(!taggingMode_)
    {
      if(saveHitTree_) hitT_->Fill();
      
      //
      for(size_t i=0; i<roi_.roiInfo_.size(); i++)
	{
	  roi_.roiInfo_[i].finalize();

	  const reco::GenJet *j=matchedGenJets[ tagJetIdx[i] ].first;
	  roiEvt_.gen_id=matchedGenJets[tagJetIdx[i] ].second;
	  roiEvt_.gen_pt=j->pt();
	  roiEvt_.gen_en=j->energy();
	  roiEvt_.gen_eta=j->eta();
	  roiEvt_.gen_phi=j->phi();
	  roiEvt_.gen_emfrac=j->emEnergy()/j->energy();
	  roiEvt_.gen_hadfrac=j->hadEnergy()/j->energy();
	  roiEvt_.gen_invfrac=j->invisibleEnergy()/j->energy();

	  for(size_t j=0; j<4; j++)
	    {
	      //per sub detector
	      for(size_t k=0; k<3; k++)
		{
		  roiEvt_.en[j][k]          = roi_.roiInfo_[i].en[j][k][0];
		  roiEvt_.eta[j][k]         = roi_.roiInfo_[i].eta[j][k][0];
		  roiEvt_.phi[j][k]         = roi_.roiInfo_[i].phi[j][k][0];
		  roiEvt_.width[j][k]       = roi_.roiInfo_[i].width[j][k][0];
		  roiEvt_.totalVolume[j][k] = roi_.roiInfo_[i].totalVolume[j][k][0];
		  roiEvt_.nhits[j][k]       = roi_.roiInfo_[i].nhits[j][k][0];

		  roiEvt_.wgt_en[j][k]          = roi_.roiInfo_[i].en[j][k][1];
		  roiEvt_.wgt_eta[j][k]         = roi_.roiInfo_[i].eta[j][k][1];
		  roiEvt_.wgt_phi[j][k]         = roi_.roiInfo_[i].phi[j][k][1];
		  roiEvt_.wgt_width[j][k]       = roi_.roiInfo_[i].width[j][k][1];
		  roiEvt_.wgt_totalVolume[j][k] = roi_.roiInfo_[i].totalVolume[j][k][1];
		  roiEvt_.wgt_nhits[j][k]       = roi_.roiInfo_[i].nhits[j][k][1];
		}

	      //per layer
	      for(size_t k=0; k<54; k++)
		{
		  roiEvt_.x[j][k]     = roi_.roiInfo_[i].x[j][k][0];
		  roiEvt_.y[j][k]     = roi_.roiInfo_[i].y[j][k][0];
		  roiEvt_.z[k]        = roi_.roiInfo_[i].z[k];
		  roiEvt_.wgt_x[j][k] = roi_.roiInfo_[i].x[j][k][1];
		  roiEvt_.wgt_y[j][k] = roi_.roiInfo_[i].y[j][k][1];
		}
	    }


	  roiT_->Fill();
	}
    }
}

//
void HGCROIAnalyzer::endJob()
{
  cout << "Analysed " << evtCtr_ << " events" << endl;
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCROIAnalyzer);
