#!/usr/bin/env python

import ROOT
from UserCode.HGCanalysis.PlotUtils import *
from array import array
import os,sys
import optparse
import commands

"""
Loops over the trees and profiles the showers
"""
def runGenLevelAnalysis(opt):
    
    #book histograms
    histos={
        'cutflow'    : ROOT.TH1F('cutflow','Selection;Events',5,0,5),
        'mjj'        : ROOT.TH1F('mjj',';Dijet invariant mass [GeV];Events',20,0,2000),
        'detajj'     : ROOT.TH1F('detajj',';#Delta#eta;Events',10,0,7),
        'etaLead'    : ROOT.TH1F('etaLead',';Pseudo-rapidity;Events',20,0,5),
        'etaTrailer' : ROOT.TH1F('etaTrailer',';Pseudo-rapidity;Events',20,0,5),
        'enLead'     : ROOT.TH1F('enLead',';Energy [GeV];Events',20,0,3000),
        'enTrailer'  : ROOT.TH1F('enTrailer',';Energy [GeV];Events',20,0,3000),
        'emFrac'     : ROOT.TH2F('emFrac',';Energy [GeV];e.m. fraction;Jets',20,0,1000,10,0,1),
        'hadFrac'    : ROOT.TH2F('hadFrac',';Energy [GeV];hadronic fraction;Jets',20,0,1000,10,0,1),
        'invFrac'    : ROOT.TH2F('invFrac',';Energy [GeV];invisible fraction;Jets',20,0,1000,10,0,1),
        'enHits005'  : ROOT.TH3F('enHits005',';HGC layer;Energy [MIP];Pseudo-rapidity;Jets',54,0,54,100,0,500,10,1.5,3),
        'enHits01'   : ROOT.TH3F('enHits01',';HGC layer;Energy [MIP];Pseudo-rapidity;Jets',54,0,54,100,0,500,10,1.5,3),
        'nHits005'   : ROOT.TH3F('nHits005',';HGC layer;# hits (>25 MIP);Pseudo-rapidity;Jets',54,0,54,10,0,10,10,1.5,3),
        'nHits01'    : ROOT.TH3F('nHits01',';HGC layer;# hits (>25 MIP);Pseudo-rapidity;Jets',54,0,54,10,0,10,10,1.5,3)
    }
    histos['cutflow'].GetXaxis().SetBinLabel(1,'Total')
    histos['cutflow'].GetXaxis().SetBinLabel(2,'#geq 2j')
    histos['cutflow'].GetXaxis().SetBinLabel(3,'M_{jj}, #Delta#eta')
    histos['cutflow'].GetXaxis().SetBinLabel(4,'1j in HGC')
    histos['cutflow'].GetXaxis().SetBinLabel(5,'2j in HGC')
    for var in histos:
        histos[var].Sumw2()
        histos[var].SetDirectory(0)

    #fill histos
    nEvts=0
    fIn=ROOT.TFile.Open(opt.input)
    HGC=fIn.Get('analysis/HGC')
    for i in xrange(0,HGC.GetEntriesFast()):
        HGC.GetEntry(i)

        if i%10==0 : drawProgressBar(float(i)/float(HGC.GetEntriesFast()))

        histos['cutflow'].Fill(0)

        #find the leading-energy jets
        jetP4=[]
        nLepMatches=0
        for nj in xrange(0,HGC.njgen):
            if HGC.genj_pt[nj]<30 : continue
            if ROOT.TMath.Abs(HGC.genj_eta[nj])>4.7 : continue

            #check if not matched to a tau
            closestMatchId=0
            closestDR=1.0
            for ngen in xrange(0,HGC.ngen):
                deta=HGC.genj_eta[nj]-HGC.gen_eta[ngen]
                dphi=ROOT.TVector2.Phi_mpi_pi(HGC.genj_phi[nj]-HGC.gen_phi[ngen])
                dr=ROOT.TMath.Sqrt(deta*deta+dphi*dphi)
                if dr>0.5 : continue
                if dr>closestDR : continue
                closestDR=dr
                closestMatchId=HGC.gen_id[ngen]
            if ROOT.TMath.Abs(closestMatchId)>=11 and ROOT.TMath.Abs(closestMatchId)<=16 : 
                nLepMatches+=1
                continue

            nselJets=len(jetP4)
            jetP4.append([nj,ROOT.TLorentzVector()])
            jetP4[nselJets][1].SetPtEtaPhiE( HGC.genj_pt[nj], HGC.genj_eta[nj], HGC.genj_phi[nj], HGC.genj_en[nj] )
        if len(jetP4)<2 : continue
        histos['cutflow'].Fill(1)
        
        #define the tag jets with highest Mjj
        tagJet1Idx, tagJet2Idx = 0, 1
        for nj in xrange(0,len(jetP4)):
            for mj in xrange(nj+1,len(jetP4)):
                newDijet=jetP4[nj][1]+jetP4[mj][1]
                dijet=jetP4[tagJet1Idx][1]+jetP4[tagJet2Idx][1]                
                if newDijet.M()<dijet.M() : 
                    continue
                tagJet1Idx=nj
                tagJet2Idx=mj

        tagJet1,tagJet2=jetP4[tagJet1Idx][1], jetP4[tagJet2Idx][1]
        dijet=tagJet1+tagJet2
        deta=ROOT.TMath.Abs(tagJet1.Eta()-tagJet2.Eta())
        histos['mjj'].Fill(dijet.M())
        histos['detajj'].Fill(deta)
        if dijet.M()<250 or deta<2.5 : continue
        histos['cutflow'].Fill(2)
        
        #requirement for jets in HGC
        njetsInHGC=0
        for jetIdx in [jetP4[tagJet1Idx][0],jetP4[tagJet2Idx][0]]:

            hitsPerLayer01=[0]*54
            hitsPerLayer005=[0]*54
            absEta=ROOT.TMath.Abs(HGC.genj_eta[jetIdx])
            if absEta>1.5 and absEta<3 : 
                njetsInHGC+=1
                histos['emFrac'].Fill(HGC.genj_en[jetIdx],HGC.genj_emfrac[jetIdx])
                histos['hadFrac'].Fill(HGC.genj_en[jetIdx],HGC.genj_hadfrac[jetIdx])
                histos['invFrac'].Fill(HGC.genj_en[jetIdx],HGC.genj_invfrac[jetIdx])

                for ihit in xrange(0,HGC.nhits):
                    deta=-HGC.genj_eta[jetIdx]-HGC.hit_eta[ihit]
                    dphi=ROOT.TVector2.Phi_mpi_pi(HGC.genj_phi[jetIdx]-HGC.hit_phi[ihit])
                    dr=ROOT.TMath.Sqrt(deta*deta+dphi*dphi)
                    
                    if dr>0.1 : continue
                    histos['enHits01'].Fill(HGC.hit_layer[ihit]-1,HGC.hit_edep[ihit],absEta)
                    if HGC.hit_edep[ihit]>50 : hitsPerLayer01[HGC.hit_layer[ihit]-1]+=1
                    
                    if dr>0.05 : continue
                    histos['enHits005'].Fill(HGC.hit_layer[ihit]-1,HGC.hit_edep[ihit],absEta)
                    if HGC.hit_edep[ihit]>50 : hitsPerLayer005[HGC.hit_layer[ihit]-1]+=1
                
                for ilay in xrange(0,len(hitsPerLayer01)):
                    histos['nHits01'].Fill(ilay,hitsPerLayer01[ilay],absEta)
                    histos['nHits005'].Fill(ilay,hitsPerLayer005[ilay],absEta)

        if njetsInHGC==1 : histos['cutflow'].Fill(3)
        if njetsInHGC==2 : histos['cutflow'].Fill(4)
        histos['etaLead'].Fill( ROOT.TMath.Max( ROOT.TMath.Abs(tagJet1.Eta()), ROOT.TMath.Abs(tagJet2.Eta()) ) )
        histos['etaTrailer'].Fill( ROOT.TMath.Min( ROOT.TMath.Abs(tagJet1.Eta()), ROOT.TMath.Abs(tagJet2.Eta()) ) )
        histos['enLead'].Fill( ROOT.TMath.Max(tagJet1.E(),tagJet2.E()) )
        histos['enTrailer'].Fill( ROOT.TMath.Min(tagJet1.E(),tagJet2.E()) )
        
    #all done here
    fIn.Close()

    #dump to output
    if opt.output!='./' : 
        os.system('mkdir -p %s'%opt.output)
    fOut=ROOT.TFile.Open('%s/VBFanalysis.root'%opt.output,'RECREATE')
    for var in histos:
        histos[var].SetDirectory(fOut)
        histos[var].Write()
    fOut.Close()


"""
steer 
"""
def main():
    
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i',      '--in' ,      dest='input',        help='Input file',              default=None)
    parser.add_option('-o',      '--out' ,     dest='output',       help='Output directory',        default='./')
    (opt, args) = parser.parse_args()
                                       
    #check inputs
    if opt.input is None:
        parser.print_help()
        sys.exit(1)

    #basic ROOT customization
    customROOTstyle()
    ROOT.gROOT.SetBatch(False)
    #ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    runGenLevelAnalysis(opt)
    
if __name__ == "__main__":
    sys.exit(main())
