#include "UserCode/HGCanalysis/plugins/HGCSimpleHitAnalyzer.h"
#include "SimDataFormats/CaloTest/interface/HcalTestNumbering.h"


using namespace std;

//
HGCSimpleHitAnalyzer::HGCSimpleHitAnalyzer(const edm::ParameterSet& iConfig) : sdH_(0)
{
  hitCollections_ = iConfig.getParameter<std::vector<std::string> >("hitCollections");
  geometrySource_ = iConfig.getParameter< std::vector<std::string> >("geometrySource");
}

//
HGCSimpleHitAnalyzer::~HGCSimpleHitAnalyzer() 
{
}

//
void HGCSimpleHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel("genParticles", genParticles);

  bool sdHneedsLabels(sdH_==0);
  if(sdHneedsLabels) 
    {
      edm::Service<TFileService> fs;
      sdH_=fs->make<TH2F>("sd",";Layer;Sub-detector;",54,0,54,hitCollections_.size(),0.,hitCollections_.size());
    }

  int nlay(0);
  for(size_t i=0; i<hitCollections_.size(); i++)
    {
      std::string hitsName(hitCollections_[i]);
      edm::Handle<edm::PCaloHitContainer> caloHits;
      iEvent.getByLabel("g4SimHits", hitsName, caloHits);

      if(sdHneedsLabels) sdH_->GetYaxis()->SetBinLabel(i+1,hitsName.c_str());

      std::string geoName(geometrySource_[i]);
      if(hitsName.find("HGCHits")!=std::string::npos)
	{
	  edm::ESHandle<HGCalGeometry> geom;
	  iSetup.get<IdealGeometryRecord>().get(geoName,geom);
	  const HGCalTopology &topo=geom->topology();
	  const HGCalDDDConstants &dddConst=topo.dddConstants();

	  uint32_t mySubDet(ForwardSubdetector::HGCEE); 
	  if(i==1) mySubDet=ForwardSubdetector::HGCHEF;  
	  else if(i==2) mySubDet=ForwardSubdetector::HGCHEB; 

	  for(edm::PCaloHitContainer::const_iterator hit_it = caloHits->begin(); hit_it != caloHits->end(); ++hit_it)
	    {
	      //use only layers which will be used at RECO level
	      HGCalDetId simId(hit_it->id());
	      float energy = hit_it->energy()*1e6;
	      float time   = hit_it->time();
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
	      for(size_t igen=0; igen<genParticles->size(); igen++)
		{
		  const reco::GenParticle & p = (*genParticles)[igen];
		  if(pos.eta()*p.eta()<0) continue;
		  countHit(layer+nlay,energy,time,pos.eta()-p.eta(),(p.charge())*deltaPhi(pos.phi(),p.phi()));
		  sdH_->Fill(layer+nlay,i);
		}
	    }

	  //increment baseline of layer counting
	  nlay += dddConst.layers(true);
	}
      else{

	for(edm::PCaloHitContainer::const_iterator hit_it = caloHits->begin(); hit_it != caloHits->end(); ++hit_it)
	  {
	    //convert sim to reco id
	    uint32_t simId( hit_it->id() );
	    int subdet, z, depth0, eta0, phi0, lay;
	    HcalTestNumbering::unpackHcalIndex(simId, subdet, z, depth0, eta0, phi0, lay);
	    int sign = (z==0) ? (-1):(1);
	    if(subdet!=int(HcalEndcap)) continue;
	    HcalDDDRecConstants::HcalID id = hcalDDD_->getHCID(subdet, eta0, phi0, lay, depth0);
	    HcalDetId detId(HcalEndcap,sign*id.eta,id.phi,id.depth);
	    float energy  = hit_it->energy()*1e6;
	    float time    = hit_it->time();

	    //get global position
	    int ieta   = detId.ietaAbs();
	    int iphi   = detId.iphi();
	    int layer  = detId.depth();
	    std::pair<double,double> etaphi = hcalDDD_->getEtaPhi(subdet,ieta,iphi);
	    double rz = hcalDDD_->getRZ(subdet,ieta,iphi);
	    HepGeom::Point3D<float> pos(rz*cos(etaphi.second)/cosh(etaphi.first),rz*sin(etaphi.second)/cosh(etaphi.first),rz*tanh(etaphi.first));
	    
	    for(size_t igen=0; igen<genParticles->size(); igen++)
	      {
		const reco::GenParticle & p = (*genParticles)[igen];
		if(pos.eta()*p.eta()<0) continue;
		countHit(layer+nlay,energy,time,pos.eta()-p.eta(),(p.charge())*deltaPhi(pos.phi(),p.phi()));
		sdH_->Fill(layer+nlay,i);
	      }
	  }
      }
    }
}

//
void HGCSimpleHitAnalyzer::countHit(int layer,float en, float time, float deta, float dphi)
{
  if(histos_["en"].find(layer)==histos_["en"].end())
    {
      edm::Service<TFileService> fs;
  
      TString pf("layer"); pf+=layer;
      histos_["en"][layer]   = fs->make<TH1F>("en_"+pf,  ";Energy [keV];Hits",100,0,1000);
      histos_["time"][layer] = fs->make<TH1F>("time_"+pf,";Time-of-arrival [ns];Hits",150,-25,50);
      histos_["deta"][layer] = fs->make<TH1F>("deta_"+pf,";#Delta#eta;Hits",100,-1,1);
      histos_["dphi"][layer] = fs->make<TH1F>("dphi_"+pf,";charge x #Delta#phi [rad];Hits",200,-2,2);
    }
  histos_["en"][layer]->Fill(en);
  histos_["time"][layer]->Fill(time);
  histos_["deta"][layer]->Fill(deta);
  histos_["dphi"][layer]->Fill(dphi);
}

// ------------ method called once each job just before starting event loop  ------------
void HGCSimpleHitAnalyzer::beginJob() 
{
}

// ------------ method called once each job just after ending the event loop  ------------
void HGCSimpleHitAnalyzer::endJob() { }

// ------------ method called when starting to processes a run  ------------
void HGCSimpleHitAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) 
{
  edm::ESHandle<HcalDDDRecConstants> pHRNDC; 
  iSetup.get<HcalRecNumberingRecord>().get( pHRNDC ); 
  hcalDDD_   = &(*pHRNDC);  
}

// ------------ method called when ending the processing of a run  ------------
void HGCSimpleHitAnalyzer::endRun(edm::Run const&, edm::EventSetup const&) 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HGCSimpleHitAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
//define this as a plug-in
DEFINE_FWK_MODULE(HGCSimpleHitAnalyzer);
