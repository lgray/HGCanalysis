import FWCore.ParameterSet.Config as cms

hgcGenParticleFilter = cms.EDFilter("CandViewSelector",
                                    src = cms.InputTag("genParticles"),
                                    cut = cms.string("((abs(pdgId())==22 && status()==1) || (abs(pdgId())<6 && status()==3)) && pt>20 && abs(eta)>1.5 && abs(eta)<3.0"),
                                    filter =  cms.bool(True)
                                    )

hgcCandFilter = cms.EDFilter("CandViewCountFilter",
                             src = cms.InputTag("hgcGenParticleFilter"),
                             minNumber = cms.uint32(1)
                             )

candInHGCSequence=cms.Sequence(hgcGenParticleFilter*hgcCandFilter)
