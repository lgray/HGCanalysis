#ifndef _HGCOccupancyAnalyzer_h_
#define _HGCOccupancyAnalyzer_h_

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
#include "TH2F.h"
#include "TString.h"

#include <string>

/**
   @class SubdetectorOccupancyHisto
   @author P. Silva (CERN)
*/
class SubdetectorOccupancyHisto
{
public:
  
  /**
     @short CTOR
   */
  SubdetectorOccupancyHisto(int sdcode,edm::Service<TFileService> *fs) : sdcode_(sdcode),fs_(fs) 
  { 
    thr_.push_back(2);
    thr_.push_back(4);
    thr_.push_back(20);
    thr_.push_back(40);
  }
  
  /**
     @short check if layer exists already
   */
  bool hasLayer(int layer) 
  { 
    return normHistos_.find(layer)!=normHistos_.end();
  }
  
  /**
     @short initiate montoring histos for a layer
   */
  void initLayer(int layer)
  {
    if(hasLayer(layer)) return;
    TString name("sd_"); name += sdcode_; name += "_layer"; name += layer;
    normHistos_[layer] = (*fs_)->make<TH1F>(name+"_nch",";Pseudo-rapidity;Number of channels",10,1.5,3.0);
    mipHistos_[layer] = (*fs_)->make<TH2F>(name+"_mip",";# MIPs;Pseudo-rapidity",250,0,250,10,1.5,3.0);
    mipHistos_[layer]->Sumw2();
    std::map<int,TH1F *> layerCountHistos;   
    std::map<int,TH2F *> layerOccHistos;
    for(size_t ithr=0; ithr<thr_.size(); ithr++)
      {
	TString thrName(name); thrName += "_thr"; thrName += thr_[ithr];
	layerCountHistos[ thr_[ithr] ] = new TH1F(thrName+"_count",";Counts;Pseudo-rapidity",10,1.5,3.0);
	layerCountHistos[ thr_[ithr] ]->SetDirectory(0);

	layerOccHistos[ thr_[ithr] ] = (*fs_)->make<TH2F>(thrName+"_occ",";Occupancy;Pseudo-rapidity",100,0,1,10,1.5,3.0);
	layerOccHistos[ thr_[ithr] ]->Sumw2();
      }
    
    countHistos_[layer]=layerCountHistos;
    occHistos_[layer]=layerOccHistos;
  }

  /**
     @short accumulate for a new hit
  */
  void count(int layer,float eta,int adc)
  {
    mipHistos_[layer]->Fill(adc*0.25,eta);
    for(size_t ithr=0; ithr<thr_.size(); ithr++)
      {
	if(adc<thr_[ithr]) continue;
	countHistos_[layer][ thr_[ithr] ]->Fill(eta);
      }
  }

  /**
     @short to be called at the end of an event
  */
  void finalize(bool reset=true)
  {
    //iterate over layers
    for(std::map< int, std::map<int,TH2F *> >::iterator it=occHistos_.begin();
	it!=occHistos_.end();
	it++)
      {
	TH1F *normH=normHistos_[it->first];

	//iterate over thresholds
	for(std::map<int,TH2F *>::iterator jt=it->second.begin();
	    jt!=it->second.end();
	    jt++)
	  {
	    TH2F *occH=jt->second;
	    TH1F *countH=countHistos_[it->first][jt->first];

	    //fill occupancy histos
	    for(int xbin=1; xbin<=countH->GetXaxis()->GetNbins(); xbin++)
	      {
		float eta=countH->GetXaxis()->GetBinCenter(xbin);
		float cts=countH->GetBinContent(xbin);
		float norm=normH->GetBinContent(xbin);
		float occ(norm>0 ? cts/norm : 0.);
		occH->Fill(occ,eta);
	      }

	    //reset counts
	    if(reset) countH->Reset("ICE");
	  }
      }
  }

  /**
     @short DTOR
   */
  ~SubdetectorOccupancyHisto() {}

  //all histos are public and can be manipulated...
  std::map< int, TH1F *>                normHistos_;
  std::map< int, TH2F *>                mipHistos_;
  std::map< int, std::map<int,TH1F *> > countHistos_;
  std::map< int, std::map<int,TH2F *> > occHistos_;
 
 private:
  
  int sdcode_;
  edm::Service<TFileService> *fs_;
  std::vector<int> thr_;
};


/**
   @class HGCOccupancyAnalyzer
   @author P. Silva (CERN)
*/

class HGCOccupancyAnalyzer : public edm::EDAnalyzer 
{
  
 public:
  
  explicit HGCOccupancyAnalyzer( const edm::ParameterSet& );
  ~HGCOccupancyAnalyzer();
  virtual void analyze( const edm::Event&, const edm::EventSetup& );

 private:
  /**
   */
  void prepareAnalysis(std::map<int,const HGCalGeometry *> &hgcGeometries);
  
  /**
     @short digi analyzers
   */
  void analyzeHEDigis(size_t isd,edm::Handle<HGCHEDigiCollection> &heDigis,const HGCalGeometry *geom);
  void analyzeEEDigis(size_t isd,edm::Handle<HGCEEDigiCollection> &eeDigis,const HGCalGeometry *geom);
  
  //gen level
  bool isInit_;

  //hgcal
  std::vector<std::string> geometrySource_;
  std::vector<std::string> digiCollections_;
  std::vector<SubdetectorOccupancyHisto> occHistos_;
};
 

#endif
