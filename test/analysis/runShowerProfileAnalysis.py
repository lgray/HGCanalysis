#!/usr/bin/env python

import ROOT
from UserCode.HGCanalysis.PlotUtils import *
from array import array
import os,sys
import optparse
import commands

"""
Draws the final plots
"""
def showProfiles(histos,norm,genEn,output):
    
    stepTitle={'sim':'SimHit','rec':'RecHit','clus':'Cluster','pf':'PF candidate'}
    
    #show the histograms
    canvas     = ROOT.TCanvas('c','c',500,500)
    profCanvas = ROOT.TCanvas('profc','profc',500,500)
    for var in histos:

        #inclusive distributions
        canvas.Clear()
        leg=ROOT.TLegend(0.2,0.8,0.5,0.94)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.035)
        varProj=[]
        for step in histos[var]:
            nproj=len(varProj)
            varProj.append( histos[var][step].ProjectionY('%s_%s_proj'%(var,step),1,histos[var][step].GetXaxis().GetNbins()) )
            varProj[nproj].SetMarkerStyle(20+nproj)
            varProj[nproj].SetMarkerColor(37+nproj)
            varProj[nproj].SetLineColor(37+nproj)
            varProj[nproj].SetTitle(stepTitle[step])
            varProj[nproj].SetDirectory(0)
            varProj[nproj].GetYaxis().SetTitle('PDF')
            fixExtremities(varProj[nproj])
            varProj[nproj].Scale(norm)
            leg.AddEntry(varProj[nproj],varProj[nproj].GetTitle(),'p')
            if nproj==0 :
                varProj[nproj].Draw('e1')
                varProj[nproj].GetYaxis().SetRangeUser(0,varProj[nproj].GetMaximum()*1.5)
            else:
                varProj[nproj].Draw('e1same')
        leg.Draw()
        simpt=MyPaveText('#bf{CMS} #it{simulation}   E=%3.1f GeV'%genEn)
        simpt.SetTextSize(0.035)
        simpt.SetTextFont(42)    

        #profiles
        profCanvas.Clear()
        profleg=ROOT.TLegend(0.2,0.8,0.5,0.94) 
        profleg.SetFillStyle(0)
        profleg.SetBorderSize(0)
        profleg.SetTextFont(42)
        profleg.SetTextSize(0.035)
        varProfs=[]
        probSum = array('d', [0.5])
        quantiles = array('d', [0.0])
        for step in histos[var]:

            nprof=len(varProfs)
            varProfs.append( ROOT.TGraphErrors() )
            varProfs[nprof].SetMarkerStyle(20+nprof)
            varProfs[nprof].SetMarkerColor(37+nprof)
            varProfs[nprof].SetLineColor(37+nprof)
            varProfs[nprof].SetTitle(stepTitle[step])
            for xbin in xrange(1, histos[var][step].GetXaxis().GetNbins() ):

                #coordinates for this projection
                xcen=histos[var][step].GetXaxis().GetBinCenter(xbin)
                xerr=0.5*histos[var][step].GetXaxis().GetBinWidth(xbin)
                
                #compute the median
                yproj=histos[var][step].ProjectionY('yproj',xbin,xbin)
                fixExtremities(yproj)
                yproj.GetQuantiles(1,quantiles,probSum)

                #add to graph
                np=varProfs[nprof].GetN()
                varProfs[nprof].SetPoint(np,xcen,quantiles[0])
                varProfs[nprof].SetPointError(np,xerr,1.253*yproj.GetMeanError())
        
            profleg.AddEntry(varProfs[nprof],varProfs[nprof].GetTitle(),'p')
            if nprof==0 :
                varProfs[nprof].Draw('ap')
                varProfs[nprof].GetYaxis().SetRangeUser(0,varProfs[nprof].GetYaxis().GetXmax()*2)
                varProfs[nprof].GetYaxis().SetTitle('Median %s'%histos[var][step].GetYaxis().GetTitle())
                varProfs[nprof].GetXaxis().SetTitle(histos[var][step].GetXaxis().GetTitle())
            else:
                varProfs[nprof].Draw('p')

        profleg.Draw()
        profpt=MyPaveText('#bf{CMS} #it{simulation}   E=%3.1f GeV'%genEn)
        profpt.SetTextSize(0.035)
        profpt.SetTextFont(42)    

        #raw_input()
        #save to png
        for c in [canvas,profCanvas]:
            c.Modified()
            c.cd()
            c.SaveAs('%s/%s_%s_en%d.png'%(output,var,c.GetName(),genEn))


"""
Loops over the trees and profiles the showers
"""
def runShowerProfileAnalysis(opt) :
    
    #prepare histograms
    steps=['sim','rec']
    baseHistos={
        'length': ROOT.TH2F('length',';Pseudo-rapidity;Length [#lambda];Events',5,1.5,3.0,10,0,15),
        'width':ROOT.TH2F('width',';HGC layer;Area [cm^{2}];Events',54,0,54,25,0,200),
        'rho':ROOT.TH2F('rho',';HGC layer;<#rho> [cm];Events',54,0,54,20,0,50),
        'volume':ROOT.TH2F('volume',';Pseudo-rapidity;Volume [cm^{2}#lambda];Events',   5,1.5,3.0,50,0,10000),
        'edep':     ROOT.TH2F('edep',   ';HGC layer;Energy [MIP];Events ',              54,0,54,100,0,1e3),
        'nhits':    ROOT.TH2F('nhits',  ';HGC layer;Number of hits;Events',             54,0,54,100,0,100),
        'sihih':    ROOT.TH2F('sihih',  ';HGC layer;#sigma(#eta,#eta);hits / event',    54,0,54,50,0,0.2),
        'sipip':    ROOT.TH2F('sipip',  ';HGC layer;#sigma(#phi,#phi);hits / event',    54,0,54,50,0,0.2),
        'sipih':    ROOT.TH2F('sipih',  ';HGC layer;#sigma(#eta,#phi);hits / event',    54,0,54,50,0,0.2)
        }
    histos={}
    for var in baseHistos:
        histos[var]={}
        baseHistos[var].Sumw2()
        for step in steps:
            histos[var][step]=baseHistos[var].Clone('%s_%s'%(var,step))
            histos[var][step].SetDirectory(0)

    #fill histos
    nEvts=0
    fIn=ROOT.TFile.Open(opt.input)
    HGC=fIn.Get('analysis/HGC')
    for i in xrange(0,HGC.GetEntriesFast()):
        HGC.GetEntry(i)

        #kinematics filter
        if HGC.genEn>opt.en+1 or HGC.genEn<opt.en-1 : continue 

        #interaction in HGC
        if HGC.hasInteractionBeforeHGC : continue

        #fully contained shower in HGC
        sumBackHEB=0
        for ilayer in [51,52,53]:
            sumBackHEB+=(getattr(HGC,'edep_sim'))[ilayer-1]
        if sumBackHEB>3 : continue

        #count selected event
        nEvts=nEvts+1

        for step in histos['length'] :
            showerLength=getattr(HGC,'totalLength_%s'%(step))
            histos['length'][step].Fill(ROOT.TMath.Abs(HGC.genEta),showerLength)
            showerVol=getattr(HGC,'totalVolume_%s'%(step))
            histos['volume'][step].Fill(ROOT.TMath.Abs(HGC.genEta),showerVol)

        for ilay in xrange(0,HGC.nlay):
            for step in histos['width']:
            
                #require some energy deposit in the layer
                if getattr(HGC,'edep_%s'%(step))[ilay]==0: continue

                #width
                layerWidth=getattr(HGC,'edepArea_%s'%(step))[ilay]
                histos['width'][step].Fill(ilay,layerWidth)
                
                #average distance
                rho=getattr(HGC,'edepdR_%s'%(step))[ilay]
                histos['rho'][step].Fill(ilay,rho)

                #other simple variables
                for var in ['edep','nhits','sihih','sipip','sipih']:
                    histos[var][step].Fill(ilay,getattr(HGC,'%s_%s'%(var,step))[ilay])

    showProfiles(histos=histos,norm=1./nEvts,genEn=opt.en,output=opt.output)


"""
steer 
"""
def main():
    
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i',      '--in' ,      dest='input',        help='Input file',              default=None)
    parser.add_option('-o',      '--out' ,     dest='output',       help='Output directory',        default='./')
    parser.add_option('-e' ,     '--en' ,      dest='en',           help='energy to profile [GeV]', default=None, type=float)
    (opt, args) = parser.parse_args()
                                       
    #check inputs
    if opt.input is None and opt.en is None:
        parser.print_help()
        sys.exit(1)

    #basic ROOT customization
    customROOTstyle()
    #ROOT.gROOT.SetBatch(False)
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    runShowerProfileAnalysis(opt)
    
if __name__ == "__main__":
    sys.exit(main())
