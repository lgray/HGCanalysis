import FWCore.ParameterSet.Config as cms

analysis_kt4 = cms.EDAnalyzer("HGCJetAnalyzer",
                              genJetsSource = cms.untracked.string("kt4GenJets"),
                              pfJetsSource  = cms.untracked.string("kt4PFJetsPandora")
                              )

analysis_ak4 = cms.EDAnalyzer("HGCJetAnalyzer",
                              genJetsSource = cms.untracked.string("ak4GenJets"),
                              pfJetsSource  = cms.untracked.string("ak4PFJetsPandora")
                              )

analysis_ak5 = cms.EDAnalyzer("HGCJetAnalyzer",
                              genJetsSource = cms.untracked.string("ak5GenJets"),
                              pfJetsSource  = cms.untracked.string("ak5PFJetsPandora")
                              )
