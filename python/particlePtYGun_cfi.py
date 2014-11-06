import FWCore.ParameterSet.Config as cms
from Configuration.Generator.PythiaUEZ2starLEPSettings_cfi import pythiaUESettingsBlock

generator = cms.EDProducer("Pythia6PtYDistGun",
                           maxEventsToPrint = cms.untracked.int32(5),
                           pythiaHepMCVerbosity = cms.untracked.bool(True),
                           pythiaPylistVerbosity = cms.untracked.int32(1),
                           PGunParameters = cms.PSet( ParticleID = cms.vint32(0),
                                                      kinematicsFile = cms.FileInPath('UserCode/HGCanalysis/data/flatptygun.root'),                   
                                                      PtBinning = cms.int32(200),
                                                      YBinning = cms.int32(100),
                                                      MinPt = cms.double(0.0),
                                                      MaxPt = cms.double(200.0),
                                                      MinY = cms.double(1.5),
                                                      MaxY = cms.double(3.0),
                                                      MinPhi = cms.double(-3.14159265359),
                                                      MaxPhi = cms.double(3.14159265359),
                                                      ),                    
                           PythiaParameters = cms.PSet( pythiaDefaultBlock,
                                                        parameterSets = cms.vstring('pythiaDefault')
                                                        )
                           )
