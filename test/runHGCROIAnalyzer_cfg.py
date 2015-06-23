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

process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )

inputs = """
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
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_24.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_25.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_26.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_27.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_28.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_29.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_3.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_30.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_31.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_33.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_34.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_35.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_36.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_37.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_38.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_39.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_4.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_40.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_41.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_42.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_44.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_45.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_46.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_47.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_48.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_49.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_5.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_50.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_6.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_10_9.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_50_1.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_50_11.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_50_18.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_50_19.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_50_2.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_50_23.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_50_25.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_50_26.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_50_28.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_50_29.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_50_3.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_50_36.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_50_37.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_50_38.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_50_40.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_50_43.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_50_44.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_50_45.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_50_49.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_50_5.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_50_50.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_50_6.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_50_7.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_50_8.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_1.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_10.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_13.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_15.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_19.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_2.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_26.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_28.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_30.root
/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/Events_211_100_33.root
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
""".split()


process.source.fileNames = cms.untracked.vstring(inputs)
print process.source.fileNames

#prepare alternative jet collections 
#from RecoJets.Configuration.RecoPFJets_cff import *
#process.ak3PFJetsPandora = ak4PFJets.clone( src=cms.InputTag('pandorapfanew'), rParam = 0.3 )
#process.load('RecoJets.Configuration.GenJetParticles_cff')
#from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets
#process.ak3GenJets = ak5GenJets.clone(rParam = 0.3)

#load the analyzer
process.TFileService = cms.Service("TFileService", fileName = cms.string(outputName))
process.load('UserCode.HGCanalysis.hgcROIAnalyzer_cfi')

#use superclusters(True) or jets (False)?
process.analysis.useSuperClustersAsROIs = cms.untracked.bool(False)

#run it
process.p = cms.Path(process.analysis)

