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
HGCOccupancyAnalyzer::HGCOccupancyAnalyzer( const edm::ParameterSet &iConfig ) : isInit_(false)
{
  //configure analyzer
  geometrySource_   = iConfig.getUntrackedParameter< std::vector<std::string> >("geometrySource");
  digiCollections_  = iConfig.getUntrackedParameter< std::vector<std::string> >("digiCollections");
}

//
HGCOccupancyAnalyzer::~HGCOccupancyAnalyzer()
{
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

  for(size_t i=0; i<digiCollections_.size(); i++)
    {
      if(digiCollections_[i].find("HE") != std::string::npos)
	{
	  edm::Handle<HGCHEDigiCollection> heDigis;
	  iEvent.getByLabel(edm::InputTag("mix",digiCollections_[i]),heDigis);
	  analyzeHEDigis(i,heDigis,hgcGeometries[i]);
	}
      else
	{
	  edm::Handle<HGCEEDigiCollection> eeDigis;
	  iEvent.getByLabel(edm::InputTag("mix",digiCollections_[i]),eeDigis);
	  analyzeEEDigis(i,eeDigis,hgcGeometries[i]);
	}
    }
}


//
void HGCOccupancyAnalyzer::analyzeHEDigis(size_t isd,edm::Handle<HGCHEDigiCollection> &heDigis, const HGCalGeometry *geom)
{
  //check inputs
  if(!heDigis.isValid()) return;

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
  
  occHistos_[isd].finalize();
}

//
void HGCOccupancyAnalyzer::analyzeEEDigis(size_t isd,edm::Handle<HGCEEDigiCollection> &eeDigis, const HGCalGeometry *geom)
{
  //check inputs
  if(!eeDigis.isValid()) return;
  
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

  occHistos_[isd].finalize();
}


//define this as a plug-in
DEFINE_FWK_MODULE(HGCOccupancyAnalyzer);
