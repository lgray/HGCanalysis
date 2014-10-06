import FWCore.ParameterSet.Config as cms

process = cms.Process('REDIGIANDMIX')

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
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

#configure from command line
import os,sys
if(len(sys.argv)<7):
    print '\ncmsRun pulseShapeTau signalTag first_file step minBiasTag outputDir\n'
    sys.exit() 
pulseShapeTau = float(sys.argv[2])
preFix        = sys.argv[3]
ffile         = int(sys.argv[4])
step          = int(sys.argv[5])
minBiasPreFix = sys.argv[6]
outputDir='./'
if len(sys.argv)>6 : outputDir=sys.argv[7]
print outputDir

# Input source
from UserCode.HGCanalysis.storeTools_cff import fillFromStore
process.source = cms.Source("PoolSource",
                            secondaryFileNames = cms.untracked.vstring(),
                            fileNames=cms.untracked.vstring()
                            )
process.source.fileNames=fillFromStore('/store/cmst3/group/hgcal/CMSSW/%s'%preFix,ffile,step)

process.options = cms.untracked.PSet(    )

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('step2 nevts:10'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    fileName = cms.untracked.string('file:%s/HGCEvents_%s_%d_%d.root'%(outputDir,preFix,pulseShapeTau,ffile)),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW')
    )
)

#mixing
process.mix.input.nbPileupEvents.averageNumber = cms.double(140.000000)
process.mix.bunchspace = cms.int32(25)
process.mix.minBunch = cms.int32(-12)
process.mix.maxBunch = cms.int32(3)
mixFileNames=fillFromStore('/store/cmst3/group/hgcal/CMSSW/%s'%minBiasPreFix,0,-1)
import random
random.shuffle(mixFileNames)
process.mix.input.fileNames = cms.untracked.vstring(mixFileNames)
process.mix.digitizers = cms.PSet(process.theDigitizersValid)

print process.source.fileNames
print process.mix.input.fileNames[0]

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'DES23_62_V1::All', '')

# Path and EndPath definitions
process.digitisation_step = cms.Path(process.pdigi_valid)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# Schedule definition
process.schedule = cms.Schedule(process.digitisation_step,process.L1simulation_step,process.digi2raw_step,process.endjob_step,process.FEVTDEBUGHLToutput_step)

# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023HGCalMuon 

#call to customisation function cust_2023HGCalMuon imported from SLHCUpgradeSimulations.Configuration.combinedCustoms
process = cust_2023HGCalMuon(process)

process.mix.digitizers.hgceeDigitizer.digiCfg.shaperTau      = cms.double(pulseShapeTau)
process.mix.digitizers.hgchebackDigitizer.digiCfg.shaperTau  = cms.double(pulseShapeTau)
process.mix.digitizers.hgchefrontDigitizer.digiCfg.shaperTau = cms.double(pulseShapeTau)


# End of customisation functions
