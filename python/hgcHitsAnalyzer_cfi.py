import FWCore.ParameterSet.Config as cms

analysis = cms.EDAnalyzer("HGCHitsAnalyzer",
                          genSource        = cms.untracked.string("genParticles"),
                          genJetsSource    = cms.untracked.string("kt4GenJets"),
                          geometrySource   = cms.untracked.vstring('HGCalEESensitive','HGCalHESiliconSensitive',  'HGCalHEScintillatorSensitive'),
                          hitCollections   = cms.untracked.vstring('HGCEERecHits',    'HGCHEFRecHits',            'HGCHEBRecHits'),
                          mipEn            = cms.untracked.vdouble(55.1,               85,                        1498.4),
                          taggingMode      = cms.untracked.bool(False), #True),
                          roipuParamFile   = cms.untracked.FileInPath('UserCode/HGCanalysis/data/ROIPUparams.root')
                          )
