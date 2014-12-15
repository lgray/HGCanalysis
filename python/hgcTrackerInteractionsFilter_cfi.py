import FWCore.ParameterSet.Config as cms

trackerIntFilter = cms.EDFilter("HGCTrackerInteractionsFilter",
                                g4TracksSource    = cms.untracked.string('g4SimHits'),
                                g4VerticesSource  = cms.untracked.string('g4SimHits')
                                )
