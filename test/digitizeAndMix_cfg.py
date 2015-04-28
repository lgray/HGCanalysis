import FWCore.ParameterSet.Config as cms

procName='DIGI'
process = cms.Process(procName)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mix_POISSON_average_cfi')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
#process.load('Configuration.StandardSequences.RawToDigi_cff')
#process.load('Configuration.StandardSequences.L1Reco_cff')
#process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

#configure from command line
import os,sys
if(len(sys.argv)<7):
    print '\ncmsRun digitizeAndMix_cfg.py SampleName First_File Step_File MinBiasName AvgPU [OutputDir]\n'
    sys.exit() 
preFix        = sys.argv[2]
ffile         = int(sys.argv[3])
step          = int(sys.argv[4])
minBiasPreFix = sys.argv[5]
avgPU         = float(sys.argv[6])
outputDir='./'
if len(sys.argv)>7 : outputDir=sys.argv[7]

# Input source
from UserCode.HGCanalysis.storeTools_cff import fillFromStore
process.source = cms.Source("PoolSource",
                            secondaryFileNames = cms.untracked.vstring(),
                            fileNames=cms.untracked.vstring()
                            )
if preFix.find('/store')>=0:
    process.source.fileNames=fillFromStore(preFix,ffile,step)
else:
    process.source.fileNames=fillFromStore('/store/cmst3/group/hgcal/CMSSW/%s'%preFix,ffile,step)
print process.source.fileNames

process.options = cms.untracked.PSet(    )

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('step2 nevts:10'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    fileName = cms.untracked.string('file:%s/Events_%d_PU%d.root'%(outputDir,ffile,avgPU)),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW')
    )
)

#mixing
process.mix.input.nbPileupEvents.averageNumber = cms.double(avgPU)
process.mix.bunchspace = cms.int32(25)
process.mix.minBunch = cms.int32(-12)
process.mix.maxBunch = cms.int32(3)
print minBiasPreFix
mixFileNames=fillFromStore('/store/cmst3/group/hgcal/CMSSW/%s'%minBiasPreFix,0,-1)
import random
random.shuffle(mixFileNames)
process.mix.input.fileNames = cms.untracked.vstring(mixFileNames)
process.mix.digitizers = cms.PSet(process.theDigitizersValid)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'DES23_62_V1::All', '')

# Path and EndPath definitions
process.digitisation_step = cms.Path(process.pdigi_valid)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.output_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(process.digitisation_step,process.L1simulation_step,process.digi2raw_step,process.endjob_step,process.output_step)

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023HGCalMuon 

#call to customisation function cust_2023HGCalMuon imported from SLHCUpgradeSimulations.Configuration.combinedCustoms
process = cust_2023HGCalMuon(process)

#custom digitization
process.mix.digitizers.hgceeDigitizer.digiCfg.feCfg.tdcResolutionInPs = cms.double( 0.00001 )
process.mix.digitizers.hgchefrontDigitizer.digiCfg.feCfg.tdcResolutionInPs = cms.double( 0.00001 )


print 'Will digitize with the following parameters'
print process.source.fileNames
print 'Sample %s starting at %s and processing %d files'%(preFix,process.source.fileNames[0],step)
print 'MinBias from %s will be used to generate <PU>=%f starting with %s'%(minBiasPreFix,avgPU,mixFileNames[0])
print 'Output will be store in %s'%process.output.fileName




# End of customisation functions
