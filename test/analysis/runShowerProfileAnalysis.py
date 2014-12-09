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
    
    stepTitle={'sim':'SimHit',
               'rec':'RecHit',
               'clus':'Cluster',
               'pf':'PF candidate'}
    
    #show the histograms
    canvas     = ROOT.TCanvas('c','c',500,500)
    profCanvas = ROOT.TCanvas('profc','profc',500,500)
    canvas2D   = ROOT.TCanvas('c2d','c2d',1000,1000)
    for var in histos:

        #show 2D
        canvas2D.Clear()
        canvas2D.Divide(2,2)
        i2dCtr=0
        for step in histos[var]:
            i2dCtr=i2dCtr+1
            canvas2D.cd(i2dCtr).SetLogz()
            histos[var][step].Draw('colz')
            MyPaveText(stepTitle[step],0.5,0.8,0.9,0.9).SetTextSize(0.05)
            if i2dCtr==0 : MyPaveText('#bf{CMS} #it{simulation}   E=%3.1f GeV'%genEn).SetTextSize(0.05)

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
            if varProj[nproj].Integral()==0 : continue
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
                if yproj.Integral()==0 : continue
                fixExtremities(yproj)
                yproj.GetQuantiles(1,quantiles,probSum)

                #add to graph
                np=varProfs[nprof].GetN()
                varProfs[nprof].SetPoint(np,xcen,quantiles[0])
                varProfs[nprof].SetPointError(np,xerr,1.253*yproj.GetMeanError())
        
            profleg.AddEntry(varProfs[nprof],varProfs[nprof].GetTitle(),'p')
            if nprof==0 :
                varProfs[nprof].Draw('ap')
                varProfs[nprof].GetYaxis().SetRangeUser(0,varProfs[nprof].GetYaxis().GetXmax()*1.5)
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
        for c in [canvas,profCanvas,canvas2D]:
            c.Modified()
            c.cd()
            c.SaveAs('%s/%s_%s_en%d.png'%(output,var,c.GetName(),genEn))

"""
"""
def getLayerLambda(layerIdx) :
    if layerIdx==0:                   return 0.01
    if layerIdx>=1 and layerIdx<=10:  return 0.036
    if layerIdx>=11 and layerIdx<=20: return 0.043
    if layerIdx>=21 and layerIdx<=29: return 0.056
    if layerIdx==30:                  return 0.338
    if layerIdx>=31 and layerIdx<=41: return 0.273
    if layerIdx>=42 and layerIdx<=53: return 0.475
    return 0
  


"""
Loops over the trees and profiles the showers
"""
def runShowerProfileAnalysis(opt,en,steps,shower_tuple):
    
    maxLayers=54
    #if opt.input.find('22')>=0 or opt.input.find('photon')>=0 or opt.input.find('11')>=0 : maxLayers=35

    #prepare histograms
    baseHistos={
        'length': ROOT.TH2F('length',';Pseudo-rapidity;Length [#lambda];Events',5,1.5,3.0,10,0,15),
        'width':ROOT.TH2F('width',';HGC layer;Area [cm^{2}];Events',maxLayers,0,maxLayers,25,0,200),
        'rho':ROOT.TH2F('rho',';HGC layer;<#rho> [cm];Events',maxLayers,0,maxLayers,20,0,50),
        'volume':ROOT.TH2F('volume',';Pseudo-rapidity;Volume [cm^{2}#lambda];Events',   5,1.5,3.0,50,0,10000),
        'edep':     ROOT.TH2F('edep',   ';HGC layer;Energy [MIP];Events ',              maxLayers,0,maxLayers,1000,0,5e3),
        'nhits':    ROOT.TH2F('nhits',  ';HGC layer;Number of hits;Events',             maxLayers,0,maxLayers,1000,0,1000),
        'sihih':    ROOT.TH2F('sihih',  ';HGC layer;#sigma(#eta,#eta);hits / event',    maxLayers,0,maxLayers,50,0,0.2),
        'sipip':    ROOT.TH2F('sipip',  ';HGC layer;#sigma(#phi,#phi);hits / event',    maxLayers,0,maxLayers,50,0,0.2),
        'sipih':    ROOT.TH2F('sipih',  ';HGC layer;#sigma(#eta,#phi);hits / event',    maxLayers,0,maxLayers,50,0,0.2)
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
    url=opt.input
    if url.find('/store/')>=0 : url='root://eoscms//eos/cms/%s'%url
    fIn=ROOT.TFile.Open(url)
    HGC=fIn.Get('analysis/HGC')
    for i in xrange(0,HGC.GetEntriesFast()):
        HGC.GetEntry(i)

        #kinematics filter
        if HGC.genEn>en+1 or HGC.genEn<en-1 : continue 

        #interaction in HGC
        if HGC.hasInteractionBeforeHGC : continue

        #fully contained shower in HGC
        sumBackHEB=0
        for ilayer in [51,52,53]:
            sumBackHEB+=(getattr(HGC,'edep_sim'))[ilayer-1]
        if sumBackHEB>3 : continue

        #count selected event
        nEvts=nEvts+1

        totalEn={}
        totalEn95={}
        length95={}
        for step in histos['length'] :
            showerLength=getattr(HGC,'totalLength_%s'%(step))
            histos['length'][step].Fill(ROOT.TMath.Abs(HGC.genEta),showerLength)
            showerVol=getattr(HGC,'totalVolume_%s'%(step))
            histos['volume'][step].Fill(ROOT.TMath.Abs(HGC.genEta),showerVol)
            totalEn[step]=0
            totalEn95[step]=0
            length95[step]=0

        for ilay in xrange(0,HGC.nlay):
            for step in steps:

                #require some energy deposit in the layer and 1 hit
                if getattr(HGC,'edep_%s'%(step))[ilay]<0.5: continue
                if getattr(HGC,'nhits_%s'%(step))[ilay]==0: continue

                totalEn[step]+=getattr(HGC,'edep_%s'%(step))[ilay]
                
                #width
                layerWidth=getattr(HGC,'edepArea_%s'%(step))[ilay]
                histos['width'][step].Fill(ilay,layerWidth)
                
                #average distance
                rho=getattr(HGC,'edepdR_%s'%(step))[ilay]
                histos['rho'][step].Fill(ilay,rho)

                #other simple variables
                for var in ['edep','nhits','sihih','sipip','sipih']:
                    histos[var][step].Fill(ilay,getattr(HGC,'%s_%s'%(var,step))[ilay])

        for ilay in xrange(0,HGC.nlay):
            for step in steps:
                #require some energy deposit in the layer and 1 hit
                if getattr(HGC,'edep_%s'%(step))[ilay]<0.5: continue
                if getattr(HGC,'nhits_%s'%(step))[ilay]==0: continue
                totalEn95[step]+=getattr(HGC,'edep_%s'%(step))[ilay]
                if totalEn95[step]>0.95*totalEn[step]: break
                length95[step]+=getLayerLambda(ilay)
        
        #add to tuple
        values=[HGC.genEn,HGC.genEta]
        for step in steps : values.append( length95[step] )
        shower_tuple.Fill(array("f",values))

    showProfiles(histos=histos,norm=1./nEvts,genEn=en,output=opt.output)
    fIn.Close()

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
    if opt.input is None and opt.en is None:
        parser.print_help()
        sys.exit(1)

    #basic ROOT customization
    customROOTstyle()
    #ROOT.gROOT.SetBatch(False)
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    #enRanges=[5, 10, 20, 40, 50, 100, 175, 250]
    enRanges=[20]
    steps=['sim','rec','clus']

    #prepare output
    fOut=ROOT.TFile.Open('%s/ShowerProfiles.root'%opt.output,'RECREATE')
    fOut.cd()
    shower_tuple=ROOT.TNtuple("shower","shower","en:eta:length95_sim:length95_rec:length95_clus")

    for ien in xrange(0,len(enRanges)):
        print 'Starting %f'%enRanges[ien]
        runShowerProfileAnalysis(opt,enRanges[ien],steps,shower_tuple)
    
    #dump to a file
    fOut.cd()
    shower_tuple.Write()
    fOut.Close()

if __name__ == "__main__":
    sys.exit(main())
