import FWCore.ParameterSet.Config as cms
from Configuration.Generator.PythiaUEZ2starLEPSettings_cfi import pythiaUESettingsBlock

generator = cms.EDProducer("Pythia6EGun",
                           PGunParameters = cms.PSet( MaxEta     = cms.double(3.0),
                                                      MinEta     = cms.double(1.5),
                                                      ParticleID = cms.vint32(0), #do not change this line without changing scripts/generateParticleGun.sh!
                                                      MinE       = cms.double(0),       #idem
                                                      MaxE       = cms.double(0),       #ibidem
                                                      MinPhi     = cms.double(-3.1415),
                                                      MaxPhi     = cms.double(3.1415),
                                                      AddAntiParticle = cms.bool(True)
                                                      ),
                           PythiaParameters = cms.PSet( pythiaUESettingsBlock,
                                                        parameterSets = cms.vstring('pythiaUESettings')
                                                        )
                           )
