#include "UserCode/HGCanalysis/interface/SlimmedRecHit.h"
#include "UserCode/HGCanalysis/interface/SlimmedCluster.h"
#include "UserCode/HGCanalysis/interface/SlimmedVertex.h"
#include "UserCode/HGCanalysis/interface/SlimmedROI.h"

#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/Common/interface/AssociationMap.h"

SlimmedRecHit srh;
std::vector<SlimmedRecHit> vsrh;
edm::Wrapper<std::vector<SlimmedRecHit> > wvsrh;

SlimmedCluster sc;
std::vector<SlimmedCluster> vsc;
edm::Wrapper<std::vector<SlimmedCluster> > wvsc;

SlimmedVertex sv;
std::vector<SlimmedVertex> vsv;
edm::Wrapper<std::vector<SlimmedVertex> > wvsv;

SlimmedROI sj;
std::vector<SlimmedROI> vsj;
edm::Wrapper<std::vector<SlimmedROI> > wvsj;
