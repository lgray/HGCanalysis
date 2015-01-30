#ifndef _HGCTrackerInteractionsFilter_h_
#define _HGCTrackerInteractionsFilter_h_

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

/**
   @class HGCTrackerInteractionsFilter
   @author P. Silva (CERN)
*/

class HGCTrackerInteractionsFilter : public edm::EDFilter
{
  
 public:
  
  explicit HGCTrackerInteractionsFilter( const edm::ParameterSet& );
  virtual ~HGCTrackerInteractionsFilter();
  virtual bool filter(edm::Event &, const edm::EventSetup &);

 private:

  //Geant4
  std::string g4TracksSource_, g4VerticesSource_;
};
 

#endif
