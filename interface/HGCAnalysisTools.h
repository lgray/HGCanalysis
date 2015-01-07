#ifndef _hgc_analysis_tools_h_
#define _hgc_analysis_tools_h_

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"

/**
   @short tries to find the interaction position based on G4 information
   @return interaction position and a flag for the interaction type
 */
struct G4InteractionPositionInfo
{
  math::XYZVectorD pos;
  int info;
};
G4InteractionPositionInfo getInteractionPosition(const std::vector<SimTrack> *SimTk, 
						 const std::vector<SimVertex> *SimVtx, 
						 int barcode);

#endif
