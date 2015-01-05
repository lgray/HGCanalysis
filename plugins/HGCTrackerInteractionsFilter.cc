#include "UserCode/HGCanalysis/plugins/HGCTrackerInteractionsFilter.h"

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
      const reco::GenParticle & p = (*genParticles)[igen];

      //sim tracks and vertices
      math::XYZVectorD hitPos=getInteractionPosition(p,SimTk,SimVtx,genBarcodes->at(igen));
      nHitsBeforeHGC=(fabs(hitPos.z())<317);
    }
  
  bool accept(nHitsBeforeHGC<maxGenParts);
  cout << "For " << maxGenParts << " analyzed found " << nHitsBeforeHGC << " interacting in tracker => decision=" << accept << endl;
  return accept;
}


//
math::XYZVectorD HGCTrackerInteractionsFilter::getInteractionPosition(const reco::GenParticle & genp,
							     edm::Handle<edm::SimTrackContainer> &SimTk,
							     edm::Handle<edm::SimVertexContainer> &SimVtx,
							     int barcode)
{
  //loop over vertices
  for (const SimVertex &simVtx : *SimVtx) 
    {
      //require the parent to be the given barcode
      bool noParent( simVtx.noParent() );
      if(noParent) continue;
      int pIdx( simVtx.parentIndex() );
      if( pIdx!=barcode) continue;

      int vtxIdx(simVtx.vertexId());
      int rawTkMult(0),eTkMult(0),gTkMult(0),nTkMult(0),nucleiTkMult(0),pTkMult(0);
      for (const SimTrack &vtxTk : *SimTk)
	{
	  int tkVtxIdx( vtxTk.vertIndex() ); 
	  if(tkVtxIdx!=vtxIdx) continue;
	  
	  int tkType=vtxTk.type();
	  rawTkMult++;
	  eTkMult      += (abs(tkType)==11);
	  gTkMult      += (abs(tkType)==22);
	  nTkMult      += (abs(tkType)==2112 || abs(tkType)==2212);
	  nucleiTkMult += (abs(tkType)>1000000000);
	  pTkMult      += (abs(tkType)==211 || abs(tkType)==111);
	}
      
      if(rawTkMult<2) continue;
      return math::XYZVectorD(simVtx.position());
    }

  return math::XYZVectorD(0,0,0);
}


//define this as a plug-in
DEFINE_FWK_MODULE(HGCTrackerInteractionsFilter);
