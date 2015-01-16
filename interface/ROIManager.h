#ifndef _roimanager_h_
#define _roimanager_h_

#include "UserCode/HGCanalysis/interface/ROIInfo.h"

#include <vector>

#include "TRandom.h"


/**
   @short handles ROI for an event
 */
class ROIManager
{
 public:

  /**
     @short CTOR
   */
 ROIManager(float acceptCone=0.6) : acceptCone_(acceptCone) { }

  /**
     @short DTOR
   */
  ~ROIManager() {}
  
  /**
     @short resets information on the ROIs
  */
  void reset() { roiInfo_.clear(); }
  
  /**
     @short add new ROI 
  */
  void add(float eta,float phi)  { roiInfo_.push_back( ROIInfo(eta,phi) ); }
 
  /**
     @short useful for testing purposes
  */
  void rotateInPhi(int roiIdx)  { roiInfo_[roiIdx].center_phi=rand_.Uniform(0,2*3.1415); }
  
  /**
     @short checks if position can be accepted by a ROI
  */
  int findROI(float eta,float phi) 
  {
    int roiIdx(-1);
    float minDR(acceptCone_);
    for(size_t i=0; i<roiInfo_.size(); i++)
      {
	float dR=deltaR(roiInfo_[i].center_eta,roiInfo_[i].center_phi,eta,phi);
	if(dR>minDR) continue;
	minDR=dR;
	roiIdx=i;
      }
    return roiIdx;
  }
 
  //ROI information is public
  std::vector<ROIInfo> roiInfo_;
 
 private:
  TRandom rand_;
  float acceptCone_;
};

#endif
	
