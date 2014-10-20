import FWCore.ParameterSet.Config as cms

analysis = cms.EDAnalyzer("HGCSimHitsAnalyzer",
                          genSource      = cms.untracked.string ("genParticles"),
                          geometrySource = cms.untracked.vstring('HGCalEESensitive','HGCalHESiliconSensitive',  'HGCalHEScintillatorSensitive'),
                          hitCollections = cms.untracked.vstring('HGCHitsEE',       'HGCHitsHEfront',           'HGCHitsHEback'               ),
                          mipEn          = cms.untracked.vdouble(55.1e-6,           85e-6,                      1498.4e-6),
                          thrList        = cms.untracked.vdouble(0.5, 1.0, 2.0, 5.0)
                          )
