#include "UserCode/HGCanalysis/plugins/HGCTrackerInteractionsFilter.h"
#include "UserCode/HGCanalysis/interface/HGCAnalysisTools.h"

using namespace std;

//
HGCTrackerInteractionsFilter::HGCTrackerInteractionsFilter( const edm::ParameterSet &iConfig )
{
  //configure analyzer
  g4TracksSource_           = iConfig.getUntrackedParameter<std::string>("g4TracksSource");
  g4VerticesSource_         = iConfig.getUntrackedParameter<std::string>("g4VerticesSource");
}

//
HGCTrackerInteractionsFilter::~HGCTrackerInteractionsFilter()
{
}

//
bool HGCTrackerInteractionsFilter::filter(edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  //generator level particles
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel("genParticles", genParticles);
  size_t maxGenParts(genParticles->size());
  if(maxGenParts>2) maxGenParts=2;
  if(genParticles->size()>maxGenParts) std::cout << "[Warning] found more than " << maxGenParts << " gen particles, will save only first " << maxGenParts << std::endl;

  //Geant4 collections
  edm::Handle<edm::SimTrackContainer> SimTk;
  iEvent.getByLabel(g4TracksSource_,SimTk);
  edm::Handle<edm::SimVertexContainer> SimVtx;
  iEvent.getByLabel(g4VerticesSource_,SimVtx); 
  edm::Handle<std::vector<int> > genBarcodes;
  iEvent.getByLabel("genParticles",genBarcodes);  
  

  //ready to roll and loop over generator level particles
  size_t nHitsBeforeHGC(0);
  for(size_t igen=0; igen<maxGenParts; igen++)
    {
      //mc truth
      //const reco::GenParticle & p = (*genParticles)[igen];

      //sim tracks and vertices      
      math::XYZVectorD hitPos=getInteractionPosition(SimTk.product(),SimVtx.product(),genBarcodes->at(igen)).pos;
      nHitsBeforeHGC=(fabs(hitPos.z())<317);
    }
  
  bool accept(nHitsBeforeHGC<maxGenParts);
  cout << "For " << maxGenParts << " analyzed found " << nHitsBeforeHGC << " interacting in tracker => decision=" << accept << endl;
  return accept;
}


//define this as a plug-in
DEFINE_FWK_MODULE(HGCTrackerInteractionsFilter);
