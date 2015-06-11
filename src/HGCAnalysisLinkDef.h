#include "UserCode/HGCanalysis/interface/HuffmanAlgo.h"
#include "UserCode/HGCanalysis/interface/HGCAnalysisTools.h"
#include "UserCode/HGCanalysis/interface/ROOTTools.h"


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
#pragma link C++ function getLambdaForHGCLayer;
#pragma link C++ struct CircleParams_t;
#pragma link C++ function fitCircleTo;
#pragma link C++ function circle_fcn;

#endif

// Local Variables:
// mode: c++
// mode: sensitive
// c-basic-offset: 8
// End:

