import FWCore.ParameterSet.Config as cms

analysis = cms.EDAnalyzer("HGCSimHitsAnalyzer",
                          genSource             = cms.untracked.string ("genParticles"),
                          geometrySource        = cms.untracked.vstring('HGCalEESensitive',        'HGCalHESiliconSensitive',  'HGCalHEScintillatorSensitive'),
                          hitCollections        = cms.untracked.vstring('HGCHitsEE',               'HGCHitsHEfront',           'HGCHitsHEback'               ),
                          recHitCollections     = cms.untracked.vstring('HGCEERecHits',            'HGCHEFRecHits',            'HGCHEBRecHits'),
                          pfClustersCollections = cms.untracked.vstring('particleFlowClusterHGCEE','particleFlowClusterHGCHEF','particleFlowClusterHGCHEB'),
                          mipEn                 = cms.untracked.vdouble(55.1,                       85,                        1498.4),
                          thrList               = cms.untracked.vdouble(0.5)
                          )
