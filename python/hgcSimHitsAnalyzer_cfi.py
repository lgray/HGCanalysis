import FWCore.ParameterSet.Config as cms

analysis = cms.EDAnalyzer("HGCSimHitsAnalyzer",
                          genSource                = cms.untracked.string ("genParticles"),
                          geometrySource           = cms.untracked.vstring('HGCalEESensitive',        'HGCalHESiliconSensitive',  'HGCalHEScintillatorSensitive'),
                          hitCollections           = cms.untracked.vstring('HGCHitsEE',               'HGCHitsHEfront',           'HGCHitsHEback'               ),
                          recHitCollections        = cms.untracked.vstring('HGCEERecHits',            'HGCHEFRecHits',            'HGCHEBRecHits'),
                          mipEn                    = cms.untracked.vdouble(55.1,                       85,                        1498.4),
                          g4TracksSource           = cms.untracked.string('g4SimHits'),
                          g4VerticesSource         = cms.untracked.string('g4SimHits'),
                          pfCandAssociationCone    = cms.untracked.double(0.3),
                          pfClusterAssociationCone = cms.untracked.double(0.3),
                          #pfClustersCollection     = cms.untracked.string('particleFlowClusterHGCEE'),
                          pfClustersCollection     = cms.untracked.string('pandorapfanew'),
                          emPFClustersCollection   = cms.untracked.string('particleFlowSuperClusterHGCEE')
                          )

