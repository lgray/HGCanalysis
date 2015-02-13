#include "UserCode/HGCanalysis/interface/HuffmanAlgo.h"
#include "UserCode/HGCanalysis/interface/HGCAnalysisTools.h"
#include "UserCode/HGCanalysis/interface/ROIInfo.h"
#include "UserCode/HGCanalysis/interface/HGCROISummary.h"


#ifdef __CINT__

#pragma link off all class; 
#pragma link off all function; 
#pragma link off all global; 
#pragma link off all typedef;

#pragma link C++ class HuffmanTreeNode;
#pragma link C++ struct HuffmanTreeNodeCmp;
#pragma link C++ class HuffmanCode;
#pragma link C++ typedef HuffmanCodeMap;
#pragma link C++ function BuildHuffmanTree;
#pragma link C++ function GenerateHuffmanCodes;
#pragma link C++ function getHuffmanCodesFrom;
#pragma link C++ function getTriggerBits;
#pragma link C++ function getReadoutBits;
#pragma link C++ function testCompressionAlgos;
#pragma link C++ struct G4InteractionPositionInfo;
#pragma link C++ function getInteractionPosition;
#pragma link C++ function getEffSigma;
#pragma link C++ class ROIInfo;
#pragma link C++ class std::vector<ROIInfo>;
#pragma link C++ getLambdaForHGCLayer;
#pragma link C++ class HGCROISummary;
#pragma link C++ function initHGCROITree;
#pragma link C++ function attachHGCROITree;


#endif

// Local Variables:
// mode: c++
// mode: sensitive
// c-basic-offset: 8
// End:

