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
def prepareQGdiscriminatorTree(opt) :

    #prepare output
    fOut = ROOT.TFile("QGTree.root","RECREATE")
    fOut.cd()

    #create ntuple
    varNames  = 'gen_pt:gen_eta:width_EE_0.05:width_EE_0.1:width_EE_0.2:width_FH_0.05:width_FH_0.1:width_FH_0.2'
    varNames += ':nhits_EE_0.05:nhits_EE_0.1:nhits_EE_0.2:nhits_FH_0.05:nhits_FH_0.1:nhits_FH_0.2'
    varNames += ':enratio_0.05to0.2:enratio_0.1to0.2'
    varNames += ':weight'
    varList =  varNames.split(':')
    TreeS=ROOT.TNtuple('TreeS','TreeS',varNames)
    TreeB=ROOT.TNtuple('TreeB','TreeB',varNames)


    #book histograms
    histos={
        'q_ptvseta' : ROOT.TH2F('q_ptvseta',';p_{T} [GeV];#eta',50,0,250,5,1.5,3.0),
        'g_ptvseta' : ROOT.TH2F('g_ptvseta',';p_{T} [GeV];#eta',50,0,250,5,1.5,3.0)
    }
    for var in histos:
        histos[var].Sumw2()
        histos[var].SetDirectory(fOut)

    #prepare chain
    HGCROI=ROOT.TChain('analysis/HGCROI')
    for f in opt.input.split(',') : HGCROI.AddFile(f)

    #first pass fill histos used to derive weights
    nentries=HGCROI.GetEntries()
    print 'Preparing weights'
    for i in xrange(0,nentries):
        HGCROI.GetEntry(i)
        if i%10==0 : drawProgressBar(float(i)/float(nentries))
        key='q_ptvseta'
        if HGCROI.gen_id==0 or HGCROI.gen_id==21 : key='g_ptvseta'
        absEta=ROOT.TMath.Abs(HGCROI.gen_eta)
        if absEta>3 or absEta<1.5 : continue
        genPt=HGCROI.gen_pt
        if genPt>250 : genPt=250
        histos[key].Fill(genPt,absEta)

    #second pass, prepare training tree
    print '\nPreparing tree'
    for i in xrange(0,nentries):
        HGCROI.GetEntry(i)
        if i%10==0 : drawProgressBar(float(i)/float(nentries))
        
        key,isquark='q_ptvseta',True
        if HGCROI.gen_id==0 or HGCROI.gen_id==21 : key,isquark='g_ptvseta',False

        #derive weight
        absEta=ROOT.TMath.Abs(HGCROI.gen_eta)
        if absEta>3 or absEta<1.5 : continue
        genPt=HGCROI.gen_pt
        if genPt>250 : genPt=250
        xbin=histos[key].GetXaxis().FindBin(genPt)
        ybin=histos[key].GetYaxis().FindBin(absEta)
        weight=1./histos[key].GetBinContent(xbin,ybin)

        en005=HGCROI.wgt_en[3*0+0]+HGCROI.wgt_en[3*0+1]
        en01=HGCROI.wgt_en[3*1+0]+HGCROI.wgt_en[3*1+1]
        en02=HGCROI.wgt_en[3*2+0]+HGCROI.wgt_en[3*2+1]
        if en02==0 : continue

        #prepare values for the ntuple
        values=[]
        for v in varList :
            if v=='gen_pt'        : values.append( HGCROI.gen_pt )
            if v=='gen_eta'       : values.append( HGCROI.gen_eta )
            if v=='weight'        : values.append( weight )
            if v=='width_EE_0.05' : values.append( HGCROI.wgt_width[3*0+0] ) 
            if v=='width_EE_0.1'  : values.append( HGCROI.wgt_width[3*1+0] ) 
            if v=='width_EE_0.2'  : values.append( HGCROI.wgt_width[3*2+0] ) 
            if v=='width_FH_0.05' : values.append( HGCROI.wgt_width[3*0+1] ) 
            if v=='width_FH_0.1'  : values.append( HGCROI.wgt_width[3*1+1] ) 
            if v=='width_FH_0.2'  : values.append( HGCROI.wgt_width[3*2+1] ) 
            if v=='nhits_EE_0.05' : values.append( HGCROI.wgt_nhits[3*0+0] ) 
            if v=='nhits_EE_0.1'  : values.append( HGCROI.wgt_nhits[3*1+0] ) 
            if v=='nhits_EE_0.2'  : values.append( HGCROI.wgt_nhits[3*2+0] ) 
            if v=='nhits_FH_0.05' : values.append( HGCROI.wgt_nhits[3*0+1] ) 
            if v=='nhits_FH_0.1'  : values.append( HGCROI.wgt_nhits[3*1+1] ) 
            if v=='nhits_FH_0.2'  : values.append( HGCROI.wgt_nhits[3*2+1] ) 
            if v=='enratio_0.05to0.2':  values.append( en005/en02 )
            if v=='enratio_0.1to0.2' :  values.append( en01/en02 )

        #fill the appropriate tree
        if isquark : TreeS.Fill( array("f",values) )
        else       : TreeB.Fill( array("f",values) )

    #all done
    print '\nS=%d B=%d'%(TreeS.GetEntriesFast(),TreeB.GetEntriesFast())
    print 'Output available in %s'%fOut.GetName()
    fOut.Write()
    fOut.Close()


"""
steer 
"""
def main():
    
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i',      '--in' ,      dest='input',        help='Input files (CSV)',              default=None)
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

    prepareQGdiscriminatorTree(opt)
    
if __name__ == "__main__":
    sys.exit(main())
