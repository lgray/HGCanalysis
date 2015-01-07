#! /usr/bin/env python
# python fastJetResol.py inputFiles=root://eoscms//eos/cms/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC21/RECO-PU0/Events_211_50_77.root,root://eoscms//eos/cms/store/cmst3/group/hgcal/CMSSW/Events_211_50_37.root

import ROOT
import sys
from DataFormats.FWLite import Events, Handle



from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.parseArguments()

# use Varparsing object
events = Events (options)

# create handle outside of loop
gjhandle = Handle ("std::vector<reco::GenJet>")
gjlabel = ("kt4GenJets")
jhandle  = Handle ("std::vector<reco::PFJet>")
jlabel = ("kt4PFJets")

# Create histograms, etc.
ROOT.gROOT.SetBatch()        # don't pop up canvases
ROOT.gROOT.SetStyle('Plain') # white background
histos={
    'gjcount' : ROOT.TH1F ("gjcount",  ";Pseudo-rapidity;Jets",50,1.5,3.0),
    'jcount'  : ROOT.TH1F ("jcount",  ";Pseudo-rapidity;Jets",50,1.5,3.0),
    'eresol'  : ROOT.TH1F ("eresol",  ";[ E(rec)-E(gen) ] / E(gen);Jets",50,-2,2),
    'ptresol' : ROOT.TH1F ("ptresol", ";[ p_{T}(rec)-p_{T}(gen) ] / p_{T}(gen);Jets",50,-2,2),
    'ptresp'  : ROOT.TH1F ("ptresp",  ";p_{T}(rec) / p_{T}(gen);Jets",50,0,2),
    'eresp'   : ROOT.TH1F ("eresp",   ";E(rec) / E(gen);Jets",50,0,2),
    'efrac'   : ROOT.TH1F ("efrac",   "electrons;Energy fraction;Jets",10,0,1.05),
    'gfrac'   : ROOT.TH1F ("gfrac",   "photons;Energy fraction;Jets",10,0,1.05),
    'nhfrac'  : ROOT.TH1F ("nhfrac",  "neutral had.;Energy fraction;Jets",10,0,1.05),
    'chfrac'  : ROOT.TH1F ("chfrac",  "charged had.;Energy fraction;Jets",10,0,1.05),
    'mufrac'  : ROOT.TH1F ("mufrac",  "muons;Energy fraction;Jets",10,0,1.05),
    'erespvsgfrac'  : ROOT.TH2F ("erespvsgfrac",  ";E(rec) / E(gen);Photon energy fraction;Jets",25,0,2,5,0,1.05),
    }
for h in histos:
    histos[h].SetDirectory(0)
    histos[h].Sumw2()

# loop over events
for event in events:
    event.getByLabel (jlabel, jhandle)
    event.getByLabel (gjlabel, gjhandle)

    #loop over gen jets
    for genjet in gjhandle.product():
        genP4=ROOT.TLorentzVector (genjet.px(), genjet.py(), genjet.pz(), genjet.energy())

        histos['gjcount'].Fill(ROOT.TMath.Abs(genjet.eta()))
        
        #matched reconstructed jet
        minDR=0.4
        matchedjet=None
        for jet in jhandle.product():
            recP4=ROOT.TLorentzVector (jet.px(), jet.py(), jet.pz(), jet.energy())
            dR=genP4.DeltaR(recP4)
            if dR<minDR:
                matchedjet=jet
                minDR=dR
        if not matchedjet: continue

        #fill histos
        histos['jcount'].Fill(ROOT.TMath.Abs(genjet.eta()))
        histos['eresol'].Fill( (genjet.energy()-matchedjet.energy()) / genjet.energy() )
        histos['erespvsgfrac'].Fill( matchedjet.energy() / genjet.energy(),matchedjet.photonEnergyFraction() )
        histos['ptresol'].Fill( (genjet.pt()-matchedjet.pt()) / genjet.pt() )
        histos['eresp'].Fill( matchedjet.energy() / genjet.energy() )
        histos['ptresp'].Fill( matchedjet.pt() / genjet.pt() )
        histos['efrac'].Fill( matchedjet.electronEnergyFraction())
        histos['chfrac'].Fill(matchedjet.chargedHadronEnergyFraction())
        histos['nhfrac'].Fill(matchedjet.neutralHadronEnergyFraction())
        histos['gfrac'].Fill(matchedjet.photonEnergyFraction())
        histos['mufrac'].Fill(matchedjet.muonEnergyFraction())


#dump to file
fOut=ROOT.TFile.Open('FastJetResol.root','RECREATE')
for h in histos:
    histos[h].Write()
fOut.Close()
