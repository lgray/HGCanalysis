import FWCore.ParameterSet.Config as cms

process = cms.Process("ROIAnalysis")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')    
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuon_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False),
                                        SkipEvent = cms.untracked.vstring('ProductNotFound')
                                        ) 

# configure from command line
# cmsRun test/runHGCSimHitsAnalyzer_cfg.py tag
# where tag can be any sub-directory under /store/cmst3/group/hgcal/CMSSW
#           or any upgrade relval sample (may need tweaking for new releases...)
import os,sys
if(len(sys.argv)<2):
    print '\ncmsRun runHGCROIAnalyzer_cfg.py input [outputfile]\n'
    sys.exit()
input2process=sys.argv[2]
outputName='HGCROIAnalyzer.root'
if(len(sys.argv)>3):
    outputName=sys.argv[3]

print '[runHGCROIAnalyzer] processing from %s and output name is %s'%(input2process,outputName)

#configure the source (list all files in directory within range [ffile,ffile+step[
from UserCode.HGCanalysis.storeTools_cff import fillFromStore
process.source = cms.Source("PoolSource",                            
                            fileNames=cms.untracked.vstring()
                            )
if input2process.find('file')>=0:
    process.source.fileNames=cms.untracked.vstring(input2process)
else :
    process.source.fileNames=fillFromStore(input2process)
print process.source.fileNames
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#prepare alternative jet collections 
#from RecoJets.Configuration.RecoPFJets_cff import *
#process.ak3PFJetsPandora = ak4PFJets.clone( src=cms.InputTag('pandorapfanew'), rParam = 0.3 )
#process.load('RecoJets.Configuration.GenJetParticles_cff')
#from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets
#process.ak3GenJets = ak5GenJets.clone(rParam = 0.3)

#load the analyzer
process.TFileService = cms.Service("TFileService", fileName = cms.string(outputName))
process.load('UserCode.HGCanalysis.hgcROIAnalyzer_cfi')

#run it
process.p = cms.Path(process.analysis)

