import FWCore.ParameterSet.Config as cms

process = cms.Process("HGCSimHitsAnalysis")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')    
process.load('FWCore.MessageService.MessageLogger_cfi')
#v5 geometry
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuon_cff')
#v4 geometry
#process.load('Configuration.Geometry.GeometryExtended2023HGCalV4MuonReco_cff')
#process.load('Configuration.Geometry.GeometryExtended2023HGCalV4Muon_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False),
                                        SkipEvent = cms.untracked.vstring('ProductNotFound')
                                        ) 

process.RandomNumberGeneratorService.analysis = cms.PSet( initialSeed = cms.untracked.uint32(12345),
                                                          engineName = cms.untracked.string('HepJamesRandom') )

# configure from command line
# cmsRun test/runHGCSimHitsAnalyzer_cfg.py tag
# where tag can be any sub-directory under /store/cmst3/group/hgcal/CMSSW
#           or any upgrade relval sample (may need tweaking for new releases...)
ffile=0
step=-1
preFix='Single13_CMSSW_6_2_0_SLHC18'
import os,sys
if(len(sys.argv)<3):
    print '\ncmsRun runHGCSimHitsAnalyzer_cfg.py doFullAnalysis tag first_file step\n'
    print '\ttag - process tag'
    print '\tfirst_file - first file to process'
    print '\tstep - number of files to process\n'
    sys.exit()

preFix=sys.argv[2]
if(len(sys.argv)>3):
    if(sys.argv[3].isdigit()) : ffile=int(sys.argv[3])
if(len(sys.argv)>4):
    if(sys.argv[4].isdigit()) : step=int(sys.argv[4])
print '[runHGCSimHitsAnalyzer] processing %d files of %s, starting from %d'%(step,preFix,ffile)

#configure the source (list all files in directory within range [ffile,ffile+step[
from UserCode.HGCanalysis.storeTools_cff import fillFromStore
process.source = cms.Source("PoolSource",                            
                            fileNames=cms.untracked.vstring()
                            )
if preFix.find('/store/')>=0 :
    process.source.fileNames=fillFromStore(preFix,ffile,step)
elif preFix.find('file')>=0:
    process.source.fileNames=cms.untracked.vstring(preFix)
elif preFix.find('lpc:')>=0:
    preFix=preFix.split(':')[1]
    process.source.fileNames=fillFromStore('srm://cmseos.fnal.gov:8443/srm/v2/server?SFN=/eos/uscms/store/user/lpchgcal/HGCAL_Samples/%s'%preFix,ffile,step)
else :
    process.source.fileNames=fillFromStore('/store/cmst3/group/hgcal/CMSSW/%s'%preFix,ffile,step)
print process.source.fileNames
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#load the analyzer
import getpass
whoami=getpass.getuser()
outputTag=preFix.replace('/','_')
process.TFileService = cms.Service("TFileService", fileName = cms.string('/tmp/%s/%s_SimHits_%d.root'%(whoami,outputTag,ffile)))
process.load('UserCode.HGCanalysis.hgcSimHitsAnalyzer_cfi')
process.load('UserCode.HGCanalysis.hgcTrackerInteractionsFilter_cfi')

#run it
process.p = cms.Path(process.analysis)
#process.p = cms.Path(process.analysis*process.trackerIntFilter)

