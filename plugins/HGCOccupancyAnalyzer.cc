#include "UserCode/HGCanalysis/plugins/HGCOccupancyAnalyzer.h"

#include "DetectorDescription/OfflineDBLoader/interface/GeometryInfoDump.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "SimG4CMS/Calo/interface/CaloHitID.h"

#include "DetectorDescription/Core/interface/DDFilter.h"
#include "DetectorDescription/Core/interface/DDFilteredView.h"
#include "DetectorDescription/Core/interface/DDSolid.h"

#include "DataFormats/GeometryVector/interface/Basic3DVector.h"

#include "CLHEP/Units/GlobalSystemOfUnits.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Plane3D.h"
#include "CLHEP/Geometry/Vector3D.h"

#include "TVector2.h"

#include <iostream>

using namespace std;

//
HGCOccupancyAnalyzer::HGCOccupancyAnalyzer( const edm::ParameterSet &iConfig ) : isInit_(false), evtSizeH_(0)
{
  //configure analyzer
  geometrySource_   = iConfig.getUntrackedParameter< std::vector<std::string> >("geometrySource");
  digiCollections_  = iConfig.getUntrackedParameter< std::vector<std::string> >("digiCollections");

  edm::Service<TFileService> fs;
  evtSizeH_ = fs->make<TH2F>("evtsize",";Event size (log_{10} byte);Sub-detector",100,0,10,geometrySource_.size()+1,0,geometrySource_.size()+1);
  evtSizeH_->GetYaxis()->SetBinLabel(1,"HGC");
  for(size_t i=0; i<geometrySource_.size(); i++)
    {
      TString subDetName(geometrySource_[i]);
      subDetName.ReplaceAll("Sensitive","");
      evtSizeH_->GetYaxis()->SetBinLabel(i+2,subDetName);
    }

}

//
HGCOccupancyAnalyzer::~HGCOccupancyAnalyzer()
{
}

//
void HGCOccupancyAnalyzer::endJob()
{
  //normalize to the number of events analyzed
  int nEvents_(0);
  for(size_t i=0; i<occHistos_.size(); i++) nEvents_=occHistos_[i].finalize();

  if(nEvents_==0) return;
  std::cout << "[endJob] analyzed " << nEvents_ << " events" << std::endl;
  evtSizeH_->Scale(1./nEvents_);
}

//
void HGCOccupancyAnalyzer::prepareAnalysis(std::map<int,const HGCalGeometry *> &hgcGeometries)
{
  if(isInit_) return;
  isInit_=true;

  //init histograms
  edm::Service<TFileService> fs;
  for(std::map<int,const HGCalGeometry *>::iterator it = hgcGeometries.begin(); it!= hgcGeometries.end(); it++)
    {
      const std::vector<DetId> &detIds=it->second->getValidDetIds();
      std::cout << "\t HGC subdet: " << it->first << " has " << detIds.size() << " valid detIds" <<  std::endl;

      //init new histogram class
      occHistos_.push_back( SubdetectorOccupancyHisto(it->first,&fs) );
      for(size_t i=0; i<detIds.size(); i++)
	{
	  int layer((detIds[i].rawId() >>19)&0x1F);
	  float eta(it->second->getPosition( detIds[i] ).eta());
	  occHistos_[it->first].initLayer(layer);
	  occHistos_[it->first].normHistos_[layer]->Fill(eta);
	}
    }
}

  
//
void HGCOccupancyAnalyzer::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  //read geometry from event setup
  std::map<int,const HGCalGeometry *> hgcGeometries;
  for(size_t i=0; i<geometrySource_.size(); i++)
    {
      edm::ESHandle<HGCalGeometry> hgcGeo;
      iSetup.get<IdealGeometryRecord>().get(geometrySource_[i],hgcGeo);
      hgcGeometries[i]=hgcGeo.product();
    }
  prepareAnalysis(hgcGeometries);

  float totalEvtSize(0);
  for(size_t i=0; i<digiCollections_.size(); i++)
    {
      if(digiCollections_[i].find("HE") != std::string::npos)
	{
	  edm::Handle<HGCHEDigiCollection> heDigis;
	  iEvent.getByLabel(edm::InputTag("mix",digiCollections_[i]),heDigis);
	  float heEvtSize=analyzeHEDigis(i,heDigis,hgcGeometries[i]);
	  totalEvtSize+=heEvtSize;
	  evtSizeH_->Fill(TMath::Log10(heEvtSize/8.),i+1);
	}
      else
	{
	  edm::Handle<HGCEEDigiCollection> eeDigis;
	  iEvent.getByLabel(edm::InputTag("mix",digiCollections_[i]),eeDigis);
	  float eeEvtSize=analyzeEEDigis(i,eeDigis,hgcGeometries[i]);
	  totalEvtSize+=eeEvtSize;
	  evtSizeH_->Fill(TMath::Log10(eeEvtSize/8.),i+1);
	}
    }
  evtSizeH_->Fill(TMath::Log10(totalEvtSize/8.),0.);
}


//
float HGCOccupancyAnalyzer::analyzeHEDigis(size_t isd,edm::Handle<HGCHEDigiCollection> &heDigis, const HGCalGeometry *geom)
{
  //check inputs
  if(!heDigis.isValid()) return 0;
  
  //analyze hits
  for(HGCHEDigiCollection::const_iterator hit_it = heDigis->begin(); hit_it != heDigis->end(); ++hit_it) 
    {
      if(hit_it->size()==0) continue;
      int adc=hit_it->sample(0).raw();      
      HGCHEDetId detId(hit_it->id());
      int layer=detId.layer();
      float eta( geom->getPosition(detId).eta());
      occHistos_[isd].count(layer,eta,adc);
    }
  
  return occHistos_[isd].endEvent();
}

//
float HGCOccupancyAnalyzer::analyzeEEDigis(size_t isd,edm::Handle<HGCEEDigiCollection> &eeDigis, const HGCalGeometry *geom)
{
  //check inputs
  if(!eeDigis.isValid()) return 0.;
  
  //analyze hits
  for(HGCEEDigiCollection::const_iterator hit_it = eeDigis->begin(); hit_it != eeDigis->end(); ++hit_it) 
    {
      if(hit_it->size()==0) continue;
      int adc=hit_it->sample(0).raw();      
      HGCEEDetId id(hit_it->id());
      int layer=id.layer();
      float eta( geom->getPosition(id).eta() );
      occHistos_[isd].count(layer,eta,adc);
    }

  return occHistos_[isd].endEvent();
}


//define this as a plug-in
DEFINE_FWK_MODULE(HGCOccupancyAnalyzer);
