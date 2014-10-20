import FWCore.ParameterSet.Config as cms

analysis = cms.EDAnalyzer("HGCSimHitsAnalyzer",
                          genSource        = cms.untracked.string("genParticles"),
                          geometrySource   = cms.untracked.vstring('HGCalEESensitive','HGCalHESiliconSensitive',  'HGCalHEScintillatorSensitive'),
                          hitCollections   = cms.untracked.vstring('HGCHitsEE',       'HGCHitsHEfront',           'HGCHitsHEback'               ),
                          )
