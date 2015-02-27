import FWCore.ParameterSet.Config as cms

from RecoLocalCalo.CaloTowersCreator.calotowermaker_cfi import *

process = cms.Process("pandoraRECO")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("RecoLocalCalo.CaloTowersCreator.calotowermaker_cfi")

# No data of type "HGCalGeometry" with label "HGCalEESensitive" in record "IdealGeometryRecord"
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')

#The below three lines were added to solve an error Exception Message:
#No "CaloGeometryRecord" record found in the EventSetup.
### process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi");
### process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
### process.load("Geometry.CaloEventSetup.CaloTopology_cfi");

#Add the next two lines to solve an error Exception Message:
#No "IdealMagneticFieldRecord" record found in the EventSetup.
### process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

#Add the next three lines to solve an error Exception Message:
#No "TransientTrackRecord" record found in the EventSetup.
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')
# process.GlobalTag.globaltag = 'START70_V1::All'

process.TrackAssociatorRecord = cms.ESSource("EmptyESSource",
        recordName = cms.string('TrackAssociatorRecord'),
        iovIsRunNotTime = cms.bool(True),
        firstValid = cms.vuint32(1)
)
process.load('SimTracker.TrackAssociation.quickTrackAssociatorByHits_cfi')

import os,sys
if(len(sys.argv)<5):
    print '\ncmsRun digitizeAndMix_cfg.py SampleName First_File Step_File [outputDir]'
    sys.exit()
preFix = sys.argv[2]
ffile  = int(sys.argv[3])
step   = int(sys.argv[4])
outputDir='./'
if len(sys.argv)>5 : outputDir=sys.argv[5]
skipThese=0
if len(sys.argv)>6 : skipThese=int(sys.argv[6])
processThese=-1
if len(sys.argv)>7 : processThese=int(sys.argv[7])

#'root://eoscms//eos/cms/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC21/RECO-PU0/Events_211_50_8.root'
#/tmp/lgray/Events_211_10_34
# Input source
from UserCode.HGCanalysis.storeTools_cff import fillFromStore
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring()
)
if preFix.find('/store')>=0:
    process.source.fileNames=fillFromStore(preFix,ffile,step)
else:
    process.source.fileNames=fillFromStore('/store/cmst3/group/hgcal/CMSSW/%s'%preFix,ffile,step)

process.options = cms.untracked.PSet()

process.source.skipEvents = cms.untracked.uint32(skipThese)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(processThese) )

print 'Will process %d events from '%processThese
print process.source.fileNames
print ' ....starting at %d'%skipThese

process.load('SimGeneral.MixingModule.mixNoPU_cfi')   
process.load('Configuration.EventContent.EventContent_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("RecoParticleFlow/PFClusterProducer/particleFlowRecHitHGCEE_cfi")

process.load('HGCal.PandoraTranslator.HGCalTrackCollection_cfi')

process.load('HGCal.PandoraTranslator.runPandora_cfi')

#process.load("UserCode/HGCanalysis/hgcTrackerInteractionsFilter_cfi")

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    fileName = cms.untracked.string('file:%s/Events_%d_%d_%d.root'%(outputDir,ffile,skipThese,processThese)),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RECO')
    )
)
process.FEVTDEBUGHLToutput.outputCommands.append('keep *')

process.kt4PFJetsPandora = process.kt4PFJets.clone( src=cms.InputTag('pandorapfanew') )
process.ak4PFJetsPandora = process.ak4PFJets.clone( src=cms.InputTag('pandorapfanew') )
process.ak5PFJetsPandora = process.ak5PFJets.clone( src=cms.InputTag('pandorapfanew') )



process.reconstruction_step = cms.Path(#process.trackerIntFilter*
                                       process.particleFlowRecHitHGCEE*
                                       process.pfTrack*
                                       process.HGCalTrackCollection*
                                       #process.trackingParticleRecoTrackAsssociation*
                                       process.pandorapfanew*
                                       process.kt4PFJetsPandora*
                                       process.ak4PFJetsPandora*
                                       process.ak5PFJetsPandora)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

process.schedule = cms.Schedule(process.reconstruction_step,
                                process.FEVTDEBUGHLToutput_step)
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023HGCalMuon
process = cust_2023HGCalMuon(process)
