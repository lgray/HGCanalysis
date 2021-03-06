import FWCore.ParameterSet.Config as cms

import sys

process = cms.Process("Analysis")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source( "PoolSource",
                             fileNames = cms.untracked.vstring(),
                             skipEvents = cms.untracked.uint32(0)
                             )

from UserCode.HGCanalysis.storeTools_cff import fillFromStore
#process.source.fileNames=fillFromStore('/store/cmst3/group/hgcal/CMSSW/SingleElectron_v15/')
process.source.fileNames=cms.untracked.vstring('/store/cmst3/group/hgcal/CMSSW/MinBias_GEN_1.root')
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

process.HepMCconvert = cms.EDAnalyzer( "MyEDMtoHEPMCAnalyzer" )
process.p = cms.Path( process.HepMCconvert )
