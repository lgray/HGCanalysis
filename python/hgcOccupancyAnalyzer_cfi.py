import FWCore.ParameterSet.Config as cms

analysis = cms.EDAnalyzer("HGCOccupancyAnalyzer",
                          geometrySource   = cms.untracked.vstring('HGCalEESensitive','HGCalHESiliconSensitive',  'HGCalHEScintillatorSensitive'),
                          digiCollections  = cms.untracked.vstring('HGCDigisEE',      'HGCDigisHEfront',          'HGCDigisHEback'              )
                          #geometrySource   = cms.untracked.vstring('HGCalHEScintillatorSensitive'),
                          #digiCollections  = cms.untracked.vstring('HGCDigisHEback'              )
                          )
