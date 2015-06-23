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

process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) )

inputs = """
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_1.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_10.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_13.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_15.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_19.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_2.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_22.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_26.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_27.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_28.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_30.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_33.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_36.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_38.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_39.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_4.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_41.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_43.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_44.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_46.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_5.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_6.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_8.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_9.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_1.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_10.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_12.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_14.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_15.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_17.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_19.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_2.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_20.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_21.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_22.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_23.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_24.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_25.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_26.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_27.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_28.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_29.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_3.root
""".split()

process.source.fileNames=cms.untracked.vstring(inputs)

print process.source.fileNames

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


