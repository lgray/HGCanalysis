import FWCore.ParameterSet.Config as cms

generator = cms.EDProducer("Pythia6JetGun",
                           PGunParameters = cms.PSet(ParticleID = cms.vint32(0),
                                                     # this defines "absolute" energy spead of particles in the jet
                                                     MinE = cms.double(0),
                                                     MaxE = cms.double(0),
                                                     # the following params define the boost
                                                     MinP = cms.double(0),
                                                     MaxP = cms.double(0),
                                                     MinPhi = cms.double(-3.1415926535),
                                                     MaxPhi = cms.double(+3.1415926535),
                                                     MinEta = cms.double(1.5),
                                                     MaxEta = cms.double(3.0)
                                                     ),
                           # no detailed pythia6 settings necessary
                           PythiaParameters = cms.PSet( parameterSets = cms.vstring() )
                           )
