#ifndef _hgc_roi_h_
#define _hgc_roi_h_

#include <vector>
#include <iostream>

#include "DataFormats/Math/interface/deltaR.h"

#include "TRandom.h"

/**
   @short define ROI in an event
 */
class HGCROI
{
 public:
  HGCROI(float roiCone=0.6) : roiCone_(roiCone) {};
  ~HGCROI() {};
  void reset() { roiCenters_.clear(); }
  void add(float eta,float phi) { roiCenters_.push_back( std::pair<float,float>(eta,phi)); }
  void rotateInPhi(int roiIdx)  { roiCenters_[roiIdx].second=rand_.Uniform(0,2*3.1415);  }
  int findROI(float eta,float phi) 
  {
    int roiIdx(-1);
    float minDR(roiCone_);
    for(size_t i=0; i<roiCenters_.size(); i++)
      {
	float dR=fabs(getDeltaR2ROI(i,eta,phi));
	if(dR>minDR) continue;
	minDR=dR;
	roiIdx=i;
      }
    return roiIdx;
  }
  float getDeltaR2ROI(int idx,float eta,float phi)
  {
    if(idx<0) return 99999.;
    float dR=deltaR(roiCenters_[idx].first,roiCenters_[idx].second,eta,phi);
    if(fabs(eta)>fabs(roiCenters_[idx].first)) dR *=-1;
    return dR;
  }
  std::pair<float,float> getROIcenter(int idx) { return roiCenters_[idx]; }
  void print()
  {
    std::cout << "Pseudo-rapidity\tPhi" << std::endl
	      << "==============================================" << std::endl;
    for(size_t i=0; i<roiCenters_.size(); i++)
      {
	std::cout << roiCenters_[i].first << "\t" << roiCenters_[i].second << std::endl;
      }
  }
 private:
  TRandom rand_;
  float roiCone_;
  std::vector<std::pair<float,float> > roiCenters_;
};

#endif
