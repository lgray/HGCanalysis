#include "UserCode/HGCanalysis/plugins/HGCROIAnalyzer.h"
#include "UserCode/HGCanalysis/interface/HGCAnalysisTools.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/JetReco/interface/TrackJet.h"
#include "DataFormats/JetReco/interface/TrackJetCollection.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "SimG4CMS/Calo/interface/CaloHitID.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"

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
#include <map>

using namespace std;

//
HGCROIAnalyzer::HGCROIAnalyzer( const edm::ParameterSet &iConfig )
{
  //configure analyzer
  genSource_              = iConfig.getUntrackedParameter<std::string>("genSource");
  genJetsSource_          = iConfig.getUntrackedParameter<std::string>("genJetsSource");
  geometrySource_         = iConfig.getUntrackedParameter< std::vector<std::string> >("geometrySource");
  hitCollections_         = iConfig.getUntrackedParameter< std::vector<std::string> >("hitCollections");
  mipEn_                  = iConfig.getUntrackedParameter< std::vector<double> >("mipEn");
  taggingMode_            = iConfig.getUntrackedParameter<bool>("taggingMode");
  saveHitTree_            = iConfig.getUntrackedParameter<bool>("saveHitTree");
  roipuParamFile_         = iConfig.getUntrackedParameter<edm::FileInPath>("roipuParamFile");
  trackJetCollection_     = iConfig.getUntrackedParameter< std::string >("trackJetCollection");
  vtxCollection_          = iConfig.getUntrackedParameter< std::string >("vtxCollection");

  edm::Service<TFileService> fs;

  //init regions of interest
  nLayerBins_=54;
  nEtaBins_=15;
  regsH_=fs->make<TH2F>("regs",";Layer;Pseudo-rapidity;",nLayerBins_,1,nLayerBins_+1,nEtaBins_,1.5,3.0);

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
  medianPU_csiH_  = (TH2F *)fIn->Get("medianPU_csi");    medianPU_csiH_->SetDirectory(0);
  widthPU_csiH_   = (TH2F *)fIn->Get("widthPU_csi");     widthPU_csiH_->SetDirectory(0);
  sigma1PU_csiH_  = (TH2F *)fIn->Get("sigma1PU_csi");    sigma1PU_csiH_->SetDirectory(0);
  sigma2PU_csiH_  = (TH2F *)fIn->Get("sigma2PU_csi");    sigma2PU_csiH_->SetDirectory(0);
  medianPU_csitH_ = (TH2F *)fIn->Get("medianPU_csit");   medianPU_csitH_->SetDirectory(0);
  widthPU_csitH_  = (TH2F *)fIn->Get("widthPU_csit");    widthPU_csitH_->SetDirectory(0);
  fIn->Close();
  
  //init tree only if not in tagging mode
  if(!taggingMode_)
    {
      roiT_=fs->make<TTree>("HGCROI","ROI Summary");
      initHGCROITree(roiT_, roiEvt_);
      //      if(saveHitTree_)
    }
  
}

//
HGCROIAnalyzer::~HGCROIAnalyzer()
{
}

//
void HGCROIAnalyzer::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  if(taggingMode_) tagEvent(iEvent,iSetup);
  else             buildROI(iEvent,iSetup);
}

//
void HGCROIAnalyzer::tagEvent( const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  //generate random phi directions
  std::vector<float> randPhi;
  for(Int_t ybin=1; ybin<=regsH_->GetYaxis()->GetNbins(); ybin++)
    randPhi.push_back( rand_.Uniform(-TMath::Pi(),TMath::Pi()) );
    
  //chacterize energy flow
  int layerCtrOffset=1;
  for(size_t i=0; i<geometrySource_.size(); i++)
    {
      edm::Handle<HGCRecHitCollection> recHits;
      iEvent.getByLabel(edm::InputTag("HGCalRecHit",hitCollections_[i]),recHits);
      
      edm::ESHandle<HGCalGeometry> geomH;
      iSetup.get<IdealGeometryRecord>().get(geometrySource_[i],geomH);
      const HGCalGeometry *geom=geomH.product();
      
      for(HGCRecHitCollection::const_iterator hit_it=recHits->begin(); hit_it!=recHits->end(); hit_it++)
	{
	  //convert energy to keV
	  float hitEn(hit_it->energy()*1e6/mipEn_[i]);
	  if(hitEn<0.5) continue;

	  uint32_t recoDetId(hit_it->id());

	  //decode position
	  const GlobalPoint refPos( std::move( geom->getPosition(recoDetId) ) );
	  int layer( ((recoDetId >> 19) & 0x1f) + layerCtrOffset-1 );
	  Int_t layerBin     = regsH_->GetXaxis()->FindBin(layer);

	  for(Int_t etaBin=1; etaBin<=regsH_->GetYaxis()->GetNbins(); etaBin++)
	    {
	      float eta=regsH_->GetYaxis()->GetBinCenter(etaBin);
	      float phi=randPhi[etaBin-1];
	      float dR=deltaR(eta,phi,refPos.eta(),refPos.phi());
	      if(dR>0.6 || dR<0.001) continue;
	      float cell_csit =  TMath::Log( hitEn/TMath::CosH(fabs(refPos.eta())) );
	      float cell_csi  =  TMath::Log( hitEn );
	      if(cell_csi>csiMax_)  cell_csi=csiMax_;
	      if(cell_csi<csiMin_)  cell_csi=csiMin_;
	      if(cell_csit>csiMax_) cell_csit=csiMax_;
	      if(cell_csit<csiMin_) cell_csit=csiMin_;
	      regsH_->Fill( layer, eta);
		  
	      float dR_ext       = dR+(layerBin-1)*(drMax_-drMin_);
	      float cell_csi_ext = cell_csi+(etaBin-1)*(csiMax_-csiMin_);
	      csidrH_->Fill(dR_ext,cell_csi_ext);
	      float cell_csit_ext = cell_csit+(etaBin-1)*(csiMax_-csiMin_);
	      csitdrH_->Fill(dR_ext,cell_csit_ext);
	    }
	}
      const HGCalTopology &topo=geom->topology();
      const HGCalDDDConstants &dddConst=topo.dddConstants();
      layerCtrOffset +=dddConst.layers(true);
    }
}


//
void HGCROIAnalyzer::buildROI( const edm::Event &iEvent, const edm::EventSetup &iSetup)
{

  //event header
  roiEvt_.newEvent(iEvent.id().run(),iEvent.luminosityBlock(),iEvent.id().event());

  //
  // GENERATOR LEVEL
  //

  //generator level particles
  edm::Handle<edm::View<reco::Candidate> > genParticles;
  iEvent.getByLabel(edm::InputTag(genSource_), genParticles);
  std::vector<int> selGenParticles;
  if(genParticles.isValid())
    {
      for(size_t i = 0; i < genParticles->size(); ++ i)
	{
	  const reco::GenParticle & p = dynamic_cast<const reco::GenParticle &>( (*genParticles)[i] );
	  if(p.status()==3)  
	    {
	      //save the generated vertex will update to the latest parton in the list
	      if(abs(p.pdgId())==21 || abs(p.pdgId())<6)
		 roiEvt_.addGenVertex( p.vertex().x(), p.vertex().y(), p.vertex().z() );

	      //add the VBF quarks (both precede the Higgs)
	      if(abs(p.pdgId())==25)
		{
		  selGenParticles.push_back(i-2);
		  selGenParticles.push_back(i-1);
		}
	      
	      continue;
	    }

	  if(p.status()!=2) continue;     //unstable partons
	  if(! ( abs(p.pdgId())==21 || abs(p.pdgId())<6) ) continue; //only care about quarks and gluons
	  selGenParticles.push_back(i);
	}
    }

  //generator level jets
  edm::Handle<reco::GenJetCollection> genJets;
  iEvent.getByLabel(edm::InputTag(genJetsSource_), genJets);
  std::vector<const reco::GenJet *> selGenJets;
  if(genJets.isValid() && genParticles.isValid())
    {
      //pre-select jets
      for(size_t k=0; k<genJets->size(); ++k)
	{
	  const reco::GenJet & j =(*genJets)[k];
	  if(j.pt()<10 || fabs(j.eta())>4.7) continue;
	  
	  //require to be matched to a quark or gluon
	  int genid(0),genstatus(0);
	  float minDR(0.4);
	  for(size_t igen=0; igen<selGenParticles.size(); igen++)
	    {
	      const reco::GenParticle & p = dynamic_cast<const reco::GenParticle &>( (*genParticles)[ selGenParticles[igen] ] );
	      float dr=deltaR(j,p);
	      if(dr>minDR) continue;

	      //if already matched to a status 3 don't care more
	      if(genstatus==3) continue;

	      minDR=dr;
	      genid=p.pdgId();
	      genstatus=p.status();
	    }
	  if(genid==0) continue;

	  selGenJets.push_back( &j );

	  //if fiducial in HGC save
	  if(abs(j.eta())<1.5 || abs(j.eta())>3.0) continue;
	  roiEvt_.addGenJet(genid,genstatus,j.pt(),j.eta(),j.phi(),j.emEnergy()/j.energy(),j.hadEnergy()/j.energy(),j.invisibleEnergy()/j.energy());
	}

      //use this to select in the fiducial VBF selection region later
      roiEvt_.gen_mjj=0;
      roiEvt_.gen_detajj=0;
      if(selGenJets.size()>=2)
	{
	  for(size_t ij=0; ij<selGenJets.size(); ij++)
	    for(size_t kj=ij+1; kj<selGenJets.size(); kj++)
	      {
		float newDijetMass=(selGenJets[ij]->p4()+selGenJets[kj]->p4()).mass();
		float newDeta=fabs(selGenJets[ij]->eta()-selGenJets[kj]->eta());
		if(newDijetMass<roiEvt_.gen_mjj) continue;
		roiEvt_.gen_mjj=newDijetMass;
		roiEvt_.gen_detajj=newDeta;
	      }
	}
    }

  //
  // RECO LEVEL
  //

  //PV collection (select all vertices within 10 mm of the generated vertex)
  edm::Handle<reco::VertexCollection> vtxH;
  iEvent.getByLabel(vtxCollection_, vtxH);
  for(size_t iv=0; iv<vtxH->size(); iv++)
    {
      const reco::Vertex &vtx=vtxH->at(iv);
      if(!vtx.isValid()) continue;
      if(vtx.isFake()) continue;
      float dzToSimVtx( fabs(vtx.position().z()-roiEvt_.gen_z) );
      if(dzToSimVtx>0.1) continue; 
      roiEvt_.addRecoVertex(iv,vtx.nTracks(),vtx.x(),vtx.y(),vtx.z(),vtx.p4().pt(),vtx.normalizedChi2());
    }

  //track jets
  edm::Handle<std::vector<reco::TrackJet> > trackJets;
  iEvent.getByLabel(trackJetCollection_,trackJets);
  for(size_t i=0; i<trackJets->size(); i++)
    {
      //pre-select the jet
      const reco::TrackJet *jet=&(trackJets->at(i));
      if(fabs(jet->eta())<1.5 || fabs(jet->eta())>3.0) continue;
      if(jet->pt()<3) continue;
      const reco::VertexRef vtx = jet->primaryVertex();
      if(vtx.isNull()) continue;
      if(!vtx->isValid()) continue;
      if(vtx->isFake()) continue;

      //check index of PV (require to be pre-selected)
      int pvIdx=-1;
      for(Int_t iv=0; iv<roiEvt_.nvtx; iv++)
	{
	  math::XYZPoint vtxDist( vtx->position().x()-roiEvt_.vtx_x[iv], 
				  vtx->position().y()-roiEvt_.vtx_y[iv], 
				  vtx->position().z()-roiEvt_.vtx_z[iv]);
	  if( vtxDist.R()>0 ) continue;
	  pvIdx=iv;
	  break;
	}
      if(pvIdx<0) continue;

      //check generated jet index
      int genIdx(-1);
      float minDR(0.4);
      for(Int_t igen=0; igen<roiEvt_.ngen; igen++)
	{
	  float dr=deltaR(jet->eta(),jet->phi(),roiEvt_.gen_eta[igen],roiEvt_.gen_phi[igen]);
	  if(dr>minDR) continue;
	  minDR=dr;
	  genIdx=igen;
	}
      
      //add information
      roiEvt_.addTrackJet(pvIdx,genIdx,jet->pt(), jet->eta(), jet->phi(), jet->numberOfTracks() );
      roiEvt_.roi_info->push_back( ROIInfo(jet->eta(),jet->phi()) );
    }

  //
  // ROI analysis
  //

  //re-center RoI based on RecHit energy sum (EE hits >25MIP)
  std::vector<float> etaSum(roiEvt_.roi_info->size(),0),phiSum(roiEvt_.roi_info->size(),0),enSum(roiEvt_.roi_info->size(),0);
  for(size_t i=0; i<1; i++)
    {
  
      edm::Handle<HGCRecHitCollection> recHits;
      iEvent.getByLabel(edm::InputTag("HGCalRecHit",hitCollections_[i]),recHits);
      
      edm::ESHandle<HGCalGeometry> geomH;
      iSetup.get<IdealGeometryRecord>().get(geometrySource_[i],geomH);
      const HGCalGeometry *geom=geomH.product();
      
      for(HGCRecHitCollection::const_iterator hit_it=recHits->begin(); hit_it!=recHits->end(); hit_it++)
	{
	  //convert energy to MIP
	  float hitEn(hit_it->energy()*1e6/mipEn_[i]);
	  if(hitEn<25) continue;
	  
	  uint32_t recoDetId(hit_it->id());
	  
	  //decode position
	  const GlobalPoint refPos( std::move( geom->getPosition(recoDetId) ) );
	  
	  for(size_t iroi=0; iroi<roiEvt_.roi_info->size(); iroi++)
	    {
	      float dr=deltaR(refPos.eta(),refPos.phi(),(*(roiEvt_.roi_info))[iroi].center_eta,(*(roiEvt_.roi_info))[iroi].center_phi);
	      if(dr>0.4) continue;
	      etaSum[iroi] += hitEn*refPos.eta();
	      phiSum[iroi] += hitEn*TVector2::Phi_mpi_pi(refPos.phi());
	      enSum[iroi]  += hitEn;
	    }
	}
    }
  
  //update RoI wrt to initial track jet estimate
  for(size_t iroi=0; iroi<roiEvt_.roi_info->size(); iroi++)
    {
      if(enSum[iroi]==0) continue;
      float newEta=etaSum[iroi]/enSum[iroi];
      float newPhi=TVector2::Phi_mpi_pi(phiSum[iroi]/enSum[iroi]);
      (*(roiEvt_.roi_info))[iroi].center_eta=newEta;
      (*(roiEvt_.roi_info))[iroi].center_phi=newPhi;
    }

  //store hits in re-centered RoIs
  Int_t layerCtrOffset(1);
  for(size_t i=0; i<geometrySource_.size(); i++)
    {
      edm::Handle<HGCRecHitCollection> recHits;
      iEvent.getByLabel(edm::InputTag("HGCalRecHit",hitCollections_[i]),recHits);
      
      edm::ESHandle<HGCalGeometry> geomH;
      iSetup.get<IdealGeometryRecord>().get(geometrySource_[i],geomH);
      const HGCalGeometry *geom=geomH.product();
      
      int hitCtr(0);
      for(HGCRecHitCollection::const_iterator hit_it=recHits->begin(); hit_it!=recHits->end(); hit_it++, hitCtr++)
	{
	  //convert energy to MIP
	  float hitEn(hit_it->energy()*1e6/mipEn_[i]);
	  if(hitEn<0.5) continue;

	  uint32_t recoDetId(hit_it->id());

	  //decode position
	  const GlobalPoint refPos( std::move( geom->getPosition(recoDetId) ) );
	  int layer( ((recoDetId >> 19) & 0x1f) + layerCtrOffset-1 );

	  //add to RoIs
	  for(size_t iroi=0; iroi<roiEvt_.roi_info->size(); iroi++)
	    {
	      //check if hit is to be considered
	      float center_eta=(*(roiEvt_.roi_info))[iroi].center_eta;
	      float center_phi=(*(roiEvt_.roi_info))[iroi].center_phi;
	      float dR=deltaR(refPos.eta(),refPos.phi(),center_eta,center_phi);
	      if(dR>0.6) continue;

	      //compute weight to assign based on layer and distance to center
	      Int_t layerBin     = regsH_->GetXaxis()->FindBin(layer);
	      Int_t etaBin       = regsH_->GetYaxis()->FindBin(TMath::Abs(center_eta));

	      //compute weights
	      float cell_csi  = dR>0.01 ? TMath::Log( hitEn ) : -9999;
	      float hitWeight(1.0),hitWeight2(1.0); 
	      if(medianPU_csiH_ && sigma1PU_csiH_ && sigma2PU_csiH_)
		{
		  Int_t dRbin=medianPU_csiH_->GetYaxis()->FindBin(dR);
		  Int_t weightsYbin( dRbin + (etaBin-1)*ndRbins_ );

		  float medianVal=medianPU_csiH_->GetBinContent(layerBin,weightsYbin);
		  float sigmaVal=sigma1PU_csiH_->GetBinContent(layerBin,weightsYbin);
		  float chi2=sigmaVal>0 ? (cell_csi>medianVal)*TMath::Power((cell_csi-medianVal)/sigmaVal,2) : 0;
		  hitWeight=chi2 > 0 ? ROOT::Math::chisquared_cdf(chi2,1) : 0;
		  
		  sigmaVal=sigma2PU_csiH_->GetBinContent(layerBin,weightsYbin);
		  chi2=sigmaVal>0 ? (cell_csi>medianVal)*TMath::Power((cell_csi-medianVal)/sigmaVal,2) : 0;
		  hitWeight2=chi2 > 0 ? ROOT::Math::chisquared_cdf(chi2,1) : 0; 
		}
	      
	      //update ROI info
	      (*(roiEvt_.roi_info))[iroi].add(hitEn,refPos.x(), refPos.y(), refPos.z(), refPos.eta(), refPos.phi(), hitWeight,hitWeight2, layer);
	    }
	}
      
      const HGCalTopology &topo=geom->topology();
      const HGCalDDDConstants &dddConst=topo.dddConstants();
      layerCtrOffset +=dddConst.layers(true);
    }
  
  
  //finalize and fill
  for(size_t iroi=0; iroi<roiEvt_.roi_info->size(); iroi++) (*(roiEvt_.roi_info))[iroi].finalize();
  roiT_->Fill();



  
  /*
  
  //reset or shuffle RoI
  if(!taggingMode_) { roi_.reset(); }
  else              { for(Int_t xbin=1; xbin<=regsH_->GetYaxis()->GetNbins(); xbin++) roi_.rotateInPhi(xbin-1); }






  //store hits in re-centered ROI
  layerCtrOffset=1;
  simEvt_.nhits=0;
  for(size_t i=0; i<geometrySource_.size(); i++)
    {
      edm::Handle<HGCRecHitCollection> recHits;
      iEvent.getByLabel(edm::InputTag("HGCalRecHit",hitCollections_[i]),recHits);
      
      edm::ESHandle<HGCalGeometry> geomH;
      iSetup.get<IdealGeometryRecord>().get(geometrySource_[i],geomH);
      const HGCalGeometry *geom=geomH.product();
      
      int hitCtr(0);
      for(HGCRecHitCollection::const_iterator hit_it=recHits->begin(); hit_it!=recHits->end(); hit_it++, hitCtr++)
	{
	  //convert energy to keV
	  float hitEn(hit_it->energy()*1e6/mipEn_[i]);
	  if(hitEn<0.5) continue;

	  uint32_t recoDetId(hit_it->id());

	  //decode position
	  const GlobalPoint refPos( std::move( geom->getPosition(recoDetId) ) );
	  int layer( ((recoDetId >> 19) & 0x1f) + layerCtrOffset-1 );

       	  //check if hits are in ROI
	  int roiIdx=roi_.findROI(refPos.eta(),refPos.phi());	

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
	  
	  float hitWeight(1.0),hitWeightT(1.0); 
	  if(medianPU_csiH_ && widthPU_csiH_)
	    {
	      Int_t dRbin=medianPU_csiH_->GetYaxis()->FindBin(dR);
	      Int_t weightsYbin( dRbin + (etaBin-1)*ndRbins_ );

	      float medianVal=medianPU_csiH_->GetBinContent(layerBin,weightsYbin);
	      //float sigmaVal=widthPU_csiH_->GetBinContent(layerBin,weightsYbin);
	      float sigmaVal=sigma1PU_csiH_->GetBinContent(layerBin,weightsYbin);
	      float chi2=sigmaVal>0 ? (cell_csi>medianVal)*TMath::Power((cell_csi-medianVal)/sigmaVal,2) : 0;
	      hitWeight=chi2 > 0 ? ROOT::Math::chisquared_cdf(chi2,1) : 0;

	      sigmaVal=sigma2PU_csiH_->GetBinContent(layerBin,weightsYbin);
	      chi2=sigmaVal>0 ? (cell_csi>medianVal)*TMath::Power((cell_csi-medianVal)/sigmaVal,2) : 0;
	      hitWeightT=chi2 > 0 ? ROOT::Math::chisquared_cdf(chi2,1) : 0; 
	    }
	  
	  // float hitWeightT(1.0);
	  // if(medianPU_csitH_ && widthPU_csitH_)
	  //   {
	  //     Int_t dRbin=medianPU_csitH_->GetYaxis()->FindBin(dR);
	  //     Int_t weightsYbin( dRbin + (etaBin-1)*ndRbins_ );

	  //     float medianVal=medianPU_csitH_->GetBinContent(layerBin,weightsYbin);
	  //     float sigmaVal=widthPU_csitH_->GetBinContent(layerBin,weightsYbin);
	  //     float chi2=sigmaVal>0 ? (cell_csi>medianVal)*TMath::Power((cell_csi-medianVal)/sigmaVal,2) : 0;
	  //     hitWeightT=chi2 > 0 ? ROOT::Math::chisquared_cdf(chi2,1) : 0;
	  //   }

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
	  roi_.roiInfo_[roiIdx].add(hitEn,refPos.x(), refPos.y(), refPos.z(), refPos.eta(), refPos.phi(), hitWeight,hitWeightT, layer);
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

	  roiEvt_.nvtx = vtxH->size();

	  roiEvt_.ncandsc=selTrackJets.size();

	  roiEvt_.roi_eta=roi_.roiInfo_[i].center_eta;
	  roiEvt_.roi_phi=roi_.roiInfo_[i].center_phi;
	  roiEvt_.roi_nsc=nMatchedTkJets[i];

	  
	  roiEvt_.tkj_pt  = roiTrackJets[i] ? roiTrackJets[i]->pt() : -1;
	  roiEvt_.tkj_eta = roiTrackJets[i] ? roiTrackJets[i]->eta() : 0;
	  roiEvt_.tkj_phi = roiTrackJets[i] ? roiTrackJets[i]->phi() : 0;
	  roiEvt_.tkj_ntk = roiTrackJets[i] ? roiTrackJets[i]->numberOfTracks() : -1; 
	  roiEvt_.tkj_vtx = roiTrackJets[i] ? roiTrackJetsVtxIndex[i] : -1;

	  const reco::GenJet *j=matchedGenJets[ tagJetIdx[i] ].first;
	  roiEvt_.gen_id=matchedGenJets[tagJetIdx[i] ].second;
	  roiEvt_.gen_pt=j->pt();
	  roiEvt_.gen_en=j->energy();
	  roiEvt_.gen_eta=j->eta();
	  roiEvt_.gen_phi=j->phi();
	  roiEvt_.gen_emfrac=j->emEnergy()/j->energy();
	  roiEvt_.gen_hadfrac=j->hadEnergy()/j->energy();
	  roiEvt_.gen_invfrac=j->invisibleEnergy()/j->energy();

	  for(size_t j=0; j<8; j++)
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

		  roiEvt_.wgt2_en[j][k]          = roi_.roiInfo_[i].en[j][k][2];
		  roiEvt_.wgt2_eta[j][k]         = roi_.roiInfo_[i].eta[j][k][2];
		  roiEvt_.wgt2_phi[j][k]         = roi_.roiInfo_[i].phi[j][k][2];
		  roiEvt_.wgt2_width[j][k]       = roi_.roiInfo_[i].width[j][k][2];
		  roiEvt_.wgt2_totalVolume[j][k] = roi_.roiInfo_[i].totalVolume[j][k][2];
		  roiEvt_.wgt2_nhits[j][k]       = roi_.roiInfo_[i].nhits[j][k][2];
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

*/
}

//
void HGCROIAnalyzer::endJob() { }


//define this as a plug-in
DEFINE_FWK_MODULE(HGCROIAnalyzer);
