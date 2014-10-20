#ifndef _HGCHitsAnalyzer_h_
#define _HGCHitsAnalyzer_h_

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DetectorDescription/Core/interface/DDCompactView.h"

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimG4CMS/Calo/interface/HGCNumberingScheme.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/HGCDigi/interface/HGCDigiCollections.h"

#include "SimCalorimetry/HGCSimProducers/interface/HGCDigitizerBase.h"  

#include "UserCode/HGCanalysis/interface/HGCSimulationEvent.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "Geometry/FCalGeometry/interface/HGCalGeometry.h"

#include "TH1F.h"
#include "TTree.h"

#include <string>

/**
   @class HGCHitsAnalyzer
   @author P. Silva (CERN)
*/

class HGCHitsAnalyzer : public edm::EDAnalyzer 
{
  
 public:
  
  explicit HGCHitsAnalyzer( const edm::ParameterSet& );
  ~HGCHitsAnalyzer();
  virtual void analyze( const edm::Event&, const edm::EventSetup& );

 private:

  /**
     @short loops over genparticles and saves a summary of stable (status=1) particles incoming to the detector
   */
  void analyzeGenParticles(edm::Handle<edm::View<reco::Candidate> > &genParticles);


  /**
     @short loops over SimTracks and SimVertices and stores the most relevant information
     cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideDataFormatSimG4Core
  */
  void analyzeG4information(edm::Handle<edm::SimTrackContainer> &SimTk,edm::Handle<edm::SimVertexContainer> &SimVtx);


  /** 
      @short analyzes tracks point to the endcap region and propagates them to the different layers of HGCal
      cf. https://cmssdt.cern.ch/SDT/lxr/source/RecoEgamma/EgammaPhotonProducers/src/ConversionProducer.cc?v=CMSSW_6_2_0_SLHC10
   */
  void analyzeTrackingInformation(edm::Handle<reco::TrackCollection> &tracks, edm::ESHandle<TrackerGeometry> &tkGeom, edm::ESHandle<MagneticField> &bField,std::map<int,const HGCalGeometry *> &hgcGeometries);
  
  /**
     @short accumulate sim hits
   */
  void analyzeHits(size_t isd,edm::Handle<edm::PCaloHitContainer> &caloHits,const HGCalGeometry *geom);

  /**
     @short save digis
   */
  void analyzeHEDigis(size_t isd,edm::Handle<HGCHEDigiCollection> &heDigis);
  void analyzeEEDigis(size_t isd,edm::Handle<HGCEEDigiCollection> &eeDigis);

  /**
     @short reset accumulators
   */
  inline void resetAccumulators()
    {
      for(size_t isd=0; isd<simHitData_.size(); isd++)
	for( HGCSimHitDataAccumulator::iterator it = simHitData_[isd].begin(); it!=simHitData_[isd].end(); it++)
	  std::fill(it->second.begin(), it->second.end(),0.);
    }

  //tree and summary ntuple
  TTree *t_;
  HGCSimEvent_t simEvt_;
  
  //gen level
  bool saveGenParticles_;
  std::string genSource_;
  
  //geant4
  bool saveG4_;
  std::string g4TracksSource_,g4VerticesSource_;

  //tracking
  bool saveTkExtrapol_;
  std::string trackSource_;
  PropagatorWithMaterial *piTkPropagator_;
  std::map<int, std::vector<ReferenceCountingPointer<BoundDisk> > > plusSurface_, minusSurface_;

  //hgcal
  std::vector<std::string> geometrySource_;
  std::vector<std::string> hitCollections_;
  std::vector<std::string> digiCollections_;
  std::vector<HGCSimHitDataAccumulator> simHitData_;
};
 

#endif
