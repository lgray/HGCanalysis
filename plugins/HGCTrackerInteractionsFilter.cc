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
      nHitsBeforeHGC=(hitPos.z()<317);
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
  for (const SimTrack &simtrack : *SimTk) 
    {
      if (simtrack.genpartIndex()!=barcode) continue;
      int simid = simtrack.trackId();
      for (const SimVertex &simvertex : *SimVtx) 
	{
	  //for neutrals only one vertex
	  if(genp.charge()==0)
	    {
	      if (simvertex.parentIndex()!=simid) continue;
	      return math::XYZVectorD(simvertex.position());
	    }
	  else
	    {
	      uint32_t tkMult(0);
	      for (const SimTrack &dausimtrack : *SimTk)
		{
		  int dausimid=dausimtrack.trackId();
		  if(dausimid==simid) continue;
		  unsigned int vtxIdx=dausimtrack.vertIndex(); 
		  if(vtxIdx!=simvertex.vertexId()) continue;
		  int tkType=abs(dausimtrack.type());
		  
		  //neglect ionization products
		  if(tkType==11) continue;

		  //check for nucleons or nuclei, neutral pions or gammas
		  if(tkType==2112 || tkType==2212 || tkType>1000000000 || tkType==111 || tkType==22) tkMult++;
		}
	      if(tkMult>2) return  math::XYZVectorD(simvertex.position());
	    }
	}
    }

  return math::XYZVectorD(0,0,0);
}


//define this as a plug-in
DEFINE_FWK_MODULE(HGCTrackerInteractionsFilter);
