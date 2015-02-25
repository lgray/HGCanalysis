#! /usr/bin/env python
# python runFastJetStudy.py inputFiles=root://eoscms//eos/cms/store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC21/RECO-PU0/Events_211_50_77.root,root://eoscms//eos/cms/store/cmst3/group/hgcal/CMSSW/Events_211_50_37.root

import ROOT
import sys
from DataFormats.FWLite import Events, Handle
import array

ROOT.gSystem.Load( "libUserCodeHGCanalysis")
ROOT.AutoLibraryLoader.enable()


from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register ('traceTkInt',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,"Activate tracing of particles interacting in tracker")
options.parseArguments()

# use Varparsing object
events = Events (options)

# create handle outside of loop
gjhandle = Handle ("std::vector<reco::GenJet>")
jhandle  = Handle ("std::vector<reco::PFJet>")
gjlabel,jlabel = ("kt4GenJets"), ("kt4PFJetsPandora")

#optional, for single particle gun samples
sthandle  = Handle('std::vector<SimTrack>')
stlabel   = ('g4SimHits')
svhandle  = Handle('std::vector<SimVertex>')
svlabel   = ('g4SimHits')
gbchandle = Handle('std::vector<int>')
gbclabel  = ('genParticles')
gchandle  = Handle('std::vector<reco::GenParticle>')
gclabel   = ('genParticles')

print 'File Run Lumi Event'

# Create histograms, etc.
ROOT.gROOT.SetBatch()        # don't pop up canvases
ROOT.gROOT.SetStyle('Plain') # white background
baseHistos={
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

histos={}
for h in baseHistos:
    histos[h]=baseHistos[h]
    histos[h].SetDirectory(0)
    histos[h].Sumw2()    
    if options.traceTkInt:
        histos[h+'_notkint']=baseHistos[h].Clone(h+'_notkint')
        histos[h+'_notkint'].SetDirectory(0)
        histos[h+'_notkint'].SetLineColor(ROOT.kRed)

# loop over events
for event in events:
    event.getByLabel(jlabel, jhandle)
    event.getByLabel(gjlabel, gjhandle)
    event.getByLabel(gclabel, gchandle)
    
    #trace which particles interact in tracker
    genCandidateInteractsInTracker=[False]*gchandle.product().size()
    if options.traceTkInt:
        event.getByLabel(stlabel,sthandle)
        event.getByLabel(svlabel,svhandle)
        event.getByLabel(gbclabel,gbchandle)
        for igp in xrange(0,gbchandle.product().size()):
            intPos=ROOT.getInteractionPosition(sthandle.product(), svhandle.product(),gbchandle.product().at(igp))
            if ROOT.TMath.Abs(intPos.pos.z())<317 : genCandidateInteractsInTracker[igp]=True

    #loop over gen jets
    for genjet in gjhandle.product():
        genP4=ROOT.TLorentzVector (genjet.px(), genjet.py(), genjet.pz(), genjet.energy())
        
        #matched gen particle
        minDR=0.4
        matchedgenp=None
        hasTkInt=False
        for igenp in xrange(0,gchandle.product().size()):
            genp=gchandle.product().at(igenp)
            p4=ROOT.TLorentzVector (genp.px(), genp.py(), genp.pz(), genp.energy())
            dR=genP4.DeltaR(p4)
            if dR<minDR:
                matchedgenp=genp
                hasTkInt=genCandidateInteractsInTracker[igenp]
                minDR=dR

        cats=['']
        if not hasTkInt and matchedgenp and options.traceTkInt: cats.append('_notkint')

        for c in cats: histos['gjcount'+c].Fill(ROOT.TMath.Abs(genjet.eta()))

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
        for c in cats:
            histos['jcount'+c].Fill(ROOT.TMath.Abs(genjet.eta()))
            histos['eresol'+c].Fill((matchedjet.energy()-genjet.energy()) / genjet.energy())
            histos['erespvsgfrac'+c].Fill( matchedjet.energy() / genjet.energy(),matchedjet.photonEnergyFraction() )
            histos['ptresol'+c].Fill( (matchedjet.pt()-genjet.pt()) / genjet.pt() )
            histos['eresp'+c].Fill( matchedjet.energy() / genjet.energy() )
            histos['ptresp'+c].Fill( matchedjet.pt() / genjet.pt() )
            histos['efrac'+c].Fill( matchedjet.electronEnergyFraction())
            histos['chfrac'+c].Fill(matchedjet.chargedHadronEnergyFraction())
            histos['nhfrac'+c].Fill(matchedjet.neutralHadronEnergyFraction())
            histos['gfrac'+c].Fill(matchedjet.photonEnergyFraction())
            histos['mufrac'+c].Fill(matchedjet.muonEnergyFraction())


        eresp=matchedjet.energy() / genjet.energy()
        if matchedgenp and  options.traceTkInt and not hasTkInt :
            if (matchedgenp.pdgId()==22 and eresp>1.4) or (ROOT.TMath.Abs(matchedgenp.pdgId())==211 and eresp>1.6) :
                evaux=event.eventAuxiliary()
                print options.inputFiles[event.fileIndex()],evaux.run(),evaux.luminosityBlock(),evaux.event()


#dump to file
fOut=ROOT.TFile.Open('FastJetResol.root','RECREATE')
for h in histos:
    histos[h].Write()
fOut.Close()
