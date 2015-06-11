import FWCore.ParameterSet.Config as cms

from Configuration.Generator.PythiaUEZ2starSettings_cfi import *

generator = cms.EDFilter("Pythia6GeneratorFilter",
                         pythiaPylistVerbosity = cms.untracked.int32(1),
                         # put here the efficiency of your filter (1. if no filter)
                         filterEfficiency = cms.untracked.double(1.0),
                         pythiaHepMCVerbosity = cms.untracked.bool(False),
                         # put here the cross section of your process (in pb)
                         crossSection = cms.untracked.double(2.41),
                         comEnergy = cms.double(14000.0),
                         maxEventsToPrint = cms.untracked.int32(3),
                         PythiaParameters = cms.PSet( pythiaUESettingsBlock,
                                                      processParameters = cms.vstring('PMAS(25,1)=125.0 !mass of Higgs', 
                                                                                      'MSEL=0           !user selection for process', 
                                                                                      'MSUB(102)=0      !ggH', 
                                                                                      'MSUB(123)=1      !ZZ fusion to H', 
                                                                                      'MSUB(124)=1      !WW fusion to H', 
                                                                                      'MSUB(24)=0       !ZH production', 
                                                                                      'MSUB(26)=0       !WH production', 
                                                                                      'MSUB(121)=0      !gg to ttH', 
                                                                                      'MSUB(122)=0      !qq to ttH', 
                                                                                      'MDME(210,1)=0    !Higgs decay into dd', 
                                                                                      'MDME(211,1)=0    !Higgs decay into uu', 
                                                                                      'MDME(212,1)=0    !Higgs decay into ss', 
                                                                                      'MDME(213,1)=0    !Higgs decay into cc', 
                                                                                      'MDME(214,1)=0    !Higgs decay into bb', 
                                                                                      'MDME(215,1)=0    !Higgs decay into tt', 
                                                                                      'MDME(216,1)=0    !Higgs decay into', 
                                                                                      'MDME(217,1)=0    !Higgs decay into Higgs decay', 
                                                                                      'MDME(218,1)=0    !Higgs decay into e nu e', 
                                                                                      'MDME(219,1)=0    !Higgs decay into mu nu mu', 
                                                                                      'MDME(220,1)=0    !Higgs decay into tau nu tau', 
                                                                                      'MDME(221,1)=0    !Higgs decay into Higgs decay', 
                                                                                      'MDME(222,1)=0    !Higgs decay into g g', 
                                                                                      'MDME(223,1)=0    !Higgs decay into gam gam', 
                                                                                      'MDME(224,1)=0    !Higgs decay into gam Z', 
                                                                                      'MDME(225,1)=1    !Higgs decay into Z Z', 
                                                                                      'MDME(226,1)=0    !Higgs decay into W W',
                                                                                      'MDME(174,1)=0    !Z decay to d dbar',
                                                                                      'MDME(175,1)=0    !Z decay to u ubar',
                                                                                      'MDME(176,1)=0    !Z decay to s sbar',
                                                                                      'MDME(177,1)=0    !Z decay to c cbar',
                                                                                      'MDME(178,1)=0    !Z decay to b bar',
                                                                                      'MDME(179,1)=0    !Z decay to t tbar',
                                                                                      'MDME(182,1)=0    !Z decay to e- e+',
                                                                                      'MDME(183,1)=1    !Z decay to nu_e nu_ebar',
                                                                                      'MDME(184,1)=0    !Z decay to mu- mu+',
                                                                                      'MDME(185,1)=1    !Z decay to nu_mu nu_mubar',
                                                                                      'MDME(186,1)=0    !Z decay to tau- tau+',
                                                                                      'MDME(187,1)=1    !Z decay to nu_tau nu_taubar'
                                                                                      ),
                                                      # This is a vector of ParameterSet names to be read, in this order
                                                      parameterSets = cms.vstring('pythiaUESettings', 
                                                                                  'processParameters')
                                                      )
                         )
