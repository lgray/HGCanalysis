import FWCore.ParameterSet.Config as cms

analysis = cms.EDAnalyzer("HGCJetAnalyzer",
                          genSource        = cms.untracked.string("genParticles"),
                          genJetsSource    = cms.untracked.string("ak4GenJets"),
                          pfJetsSource     = cms.untracked.string("ak4PFJets"),
                          eeRecHitsSource  = cms.untracked.string('HGCEERecHits'),
                          hefRecHitsSource = cms.untracked.string('HGCHEFRecHits')
                          )

