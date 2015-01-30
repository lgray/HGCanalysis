import FWCore.ParameterSet.Config as cms

analysis = cms.EDAnalyzer("HGCROIAnalyzer",
                          genSource        = cms.untracked.string("genParticles"),
                          genJetsSource    = cms.untracked.string("kt4GenJets"),
                          geometrySource   = cms.untracked.vstring('HGCalEESensitive','HGCalHESiliconSensitive',  'HGCalHEScintillatorSensitive'),
                          hitCollections   = cms.untracked.vstring('HGCEERecHits',    'HGCHEFRecHits',            'HGCHEBRecHits'),
                          mipEn            = cms.untracked.vdouble(55.1,               85,                        1498.4),
                          trackJetCollection  = cms.untracked.string('ak5TrackJets'),
                          taggingMode      = cms.untracked.bool(False), #True),
                          saveHitTree      = cms.untracked.bool(False),
                          roipuParamFile   = cms.untracked.FileInPath('UserCode/HGCanalysis/data/ROIPUparams.root')
                          )
