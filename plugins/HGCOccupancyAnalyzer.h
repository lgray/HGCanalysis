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
 SubdetectorOccupancyHisto(int sdcode,edm::Service<TFileService> *fs) : sdcode_(sdcode),fs_(fs), nEvents_(0)
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
    mipHistos_[layer] = (*fs_)->make<TH2F>(name+"_mip",";# MIPs;Pseudo-rapidity",100,0,50,10,1.5,3.0);
    mipHistos_[layer]->Sumw2();

    dataCountHistos_[layer] = new TH1F(name+"_datacount",";Data volume [bit];Pseudo-rapidity",10,1.5,3.0);
    dataCountHistos_[layer]->SetDirectory(0);
    dataVolHistos_[layer] = (*fs_)->make<TH2F>(name+"_datavol",";Data volume [bit];Pseudo-rapidity",60,0,12,10,1.5,3.0);
    dataVolHistos_[layer]->Sumw2();

    trigCountHistos_[layer] = new TH1F(name+"_trigcount",";Trigger volume [bit];Pseudo-rapidity",10,1.5,3.0);
    trigCountHistos_[layer]->SetDirectory(0);
    trigVolHistos_[layer] = (*fs_)->make<TH2F>(name+"_trigvol",";Trigger volume [bit];Pseudo-rapidity",60,0,6,10,1.5,3.0);
    trigVolHistos_[layer]->Sumw2();

    std::map<int,TH1F *> layerCountHistos;   
    std::map<int,TH2F *> layerOccHistos;

    Double_t occbins[100];
    for(size_t ibin=0; ibin<10; ibin++)   occbins[ibin]=ibin*1e-4;
    for(size_t ibin=10; ibin<19; ibin++)  occbins[ibin]=(ibin-9.)*1e-3;
    for(size_t ibin=19; ibin<29; ibin++)  occbins[ibin]=(ibin-18)*1e-2;
    for(size_t ibin=29; ibin<100; ibin++) occbins[ibin]=0.1+(1-0.1)/(100-29)*(ibin-28);
    for(size_t ithr=0; ithr<thr_.size(); ithr++)
      {
	TString thrName(name); thrName += "_thr"; thrName += thr_[ithr];
	layerCountHistos[ thr_[ithr] ] = new TH1F(thrName+"_count",";Counts;Pseudo-rapidity",10,1.5,3.0);
	layerCountHistos[ thr_[ithr] ]->SetDirectory(0);
	
	//layerOccHistos[ thr_[ithr] ] = (*fs_)->make<TH2F>(thrName+"_occ",";Occupancy;Pseudo-rapidity",100,0,1,10,1.5,3.0);
	layerOccHistos[ thr_[ithr] ] = (*fs_)->make<TH2F>(thrName+"_occ",";Occupancy;Pseudo-rapidity",99,occbins,10,1.5,3.0);
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
    dataCountHistos_[layer]->Fill(eta,computeDataVol(adc,eta));
    trigCountHistos_[layer]->Fill(eta,computeTriggerVol(adc,eta));
    for(size_t ithr=0; ithr<thr_.size(); ithr++)
      {
	if(adc<thr_[ithr]) continue;
	countHistos_[layer][ thr_[ithr] ]->Fill(eta);
      }
  }

  /**
     @short computes number of bits to send
   */
  int computeDataVol(float adc,float eta)
  {
    float nmips=0.25*adc;
    int nbits(1);
    if(nmips>6.4)      nbits=12;
    else if(nmips>0.4) nbits=8;
    return nbits;
  }
  int computeTriggerVol(float adc,float eta)
  {
    float nmips=0.25*adc;
    int nbits(1);
    if(fabs(eta)<2)
      {
	if(nmips>96)      nbits=10;
        else if(nmips>10) nbits=6;
      }
    else if(fabs(eta)<2.5)
      {
	if(nmips>96)      nbits=10;
	else if(nmips>10) nbits=6;
      }
    else
      {
	if(nmips>192)     nbits=10;
	else if(nmips>25) nbits=6;
      }
    return nbits;
  }

  /**
     @short to be called at the end of an event
  */
  float endEvent(bool reset=true)
  {
    nEvents_++;

    //iterate over layers
    float totalEvtSize(0);
    for(std::map< int, std::map<int,TH2F *> >::iterator it=occHistos_.begin();
	it!=occHistos_.end();
	it++)
      {

	//occupancies
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
	
	//data volumes
	TH1F *dataH=dataCountHistos_[it->first];
	TH1F *trigH=trigCountHistos_[it->first];
	for(int xbin=1; xbin<=dataH->GetXaxis()->GetNbins(); xbin++)
	  {
	    float ncells=normH->GetBinContent(xbin);
	    float eta=dataH->GetXaxis()->GetBinCenter(xbin);
	    float data=dataH->GetBinContent(xbin);
	    totalEvtSize += data;
	    dataVolHistos_[it->first]->Fill(ncells>0 ? data/ncells : 1. ,eta);
	    float trig=trigH->GetBinContent(xbin);
	    trigVolHistos_[it->first]->Fill(ncells>0 ? trig/ncells : 1. ,eta);
	  }
	dataH->Reset("ICE");
	trigH->Reset("ICE");
      }

    return  totalEvtSize;
  }

  /**
     @short normalize according to the number of events analyzed
   */
  int finalize()
  {
    if(nEvents_==0) return 0;
    for(std::map< int, TH2F *>::iterator it=mipHistos_.begin();
	it!=mipHistos_.end();
	it++)
      {
	//for the energy spectrum divide by the number of cells analyzed
	TH1F *normH=normHistos_[it->first];
	for(int ybin=1; ybin<=it->second->GetYaxis()->GetNbins(); ybin++)
	  {
	    float ncells=normH->GetBinContent(ybin);
	    if(ncells==0) continue;
	    for(int xbin=0; xbin<=it->second->GetXaxis()->GetNbins()+1; xbin++)
	      {
		float counts=it->second->GetBinContent(xbin,ybin);
		float countsErr=it->second->GetBinError(xbin,ybin);
		it->second->SetBinContent(xbin,ybin,counts/(nEvents_*ncells));
		it->second->SetBinError(xbin,ybin,countsErr/(nEvents_*ncells));
	      }
	  }

	//normalized by the number of events (already averaged per cell)
	dataVolHistos_[it->first]->Scale(1./nEvents_);
	trigVolHistos_[it->first]->Scale(1./nEvents_);
	for(std::map<int,TH2F *>::iterator jt=occHistos_[it->first].begin();
	    jt!=occHistos_[it->first].end();
	    jt++)
	  {
	    jt->second->Scale(1./nEvents_);
	  }
      }
    return nEvents_;
  }

  /**
     @short DTOR
   */
  ~SubdetectorOccupancyHisto() {}

  //all histos are public and can be manipulated...
  std::map< int, TH1F *>                normHistos_, dataCountHistos_, trigCountHistos_;
  std::map< int, TH2F *>                mipHistos_,  dataVolHistos_,   trigVolHistos_;
  std::map< int, std::map<int,TH1F *> > countHistos_;
  std::map< int, std::map<int,TH2F *> > occHistos_;
 
 private:
  
  int sdcode_;
  edm::Service<TFileService> *fs_;
  std::vector<int> thr_;
  int nEvents_;
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
  void endJob();

 private:
  /**
   */
  void prepareAnalysis(std::map<int,const HGCalGeometry *> &hgcGeometries);
  
  /**
     @short digi analyzers
   */
  float analyzeHEDigis(size_t isd,edm::Handle<HGCHEDigiCollection> &heDigis,const HGCalGeometry *geom);
  float analyzeEEDigis(size_t isd,edm::Handle<HGCEEDigiCollection> &eeDigis,const HGCalGeometry *geom);
  
  //gen level
  bool isInit_;

  //hgcal
  std::vector<std::string> geometrySource_;
  std::vector<std::string> digiCollections_;
  std::vector<SubdetectorOccupancyHisto> occHistos_;

  //event size 
  TH2F *evtSizeH_;
};
 

#endif
