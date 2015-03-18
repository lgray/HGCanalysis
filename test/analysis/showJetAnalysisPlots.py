#!/usr/bin/env python

import ROOT
import sys
import math
import optparse
import os
from UserCode.HGCanalysis.PlotUtils import *

ALLPLOTS=[]

"""
"""
def showEfficiencyMap(algo,fIn,outputDir,hasPU) :

    gen=fIn.Get('analysis_%s/kin_gen'%algo)
    recmatch=fIn.Get('analysis_%s/kin_gen_matched'%algo)
    rec=fIn.Get('analysis_%s/kin_rec'%algo)
    rec_unm=fIn.Get('analysis_%s/kin_rec_unm'%algo)

    #efficiency map
    c=ROOT.TCanvas('c','c',500,500)
    c.SetRightMargin(0.15)
    eff=recmatch.Clone('rec_eff')
    eff.Divide(gen)
    eff.Draw('colz')
    eff.GetZaxis().SetTitle("Reconstruction efficiency     ");
    eff.GetZaxis().SetRangeUser(0,1)
    eff.GetZaxis().SetTitleOffset(-0.5)
    eff.GetYaxis().SetTitleOffset(1.2)
    MyPaveText('#bf{CMS} #it{simulation}    %s jets'%algo)
    c.Modified()
    c.Update()
    c.SaveAs('%s/receff_%s.png'%(outputDir,algo))

    #show projected kinematics
    for i in xrange(0,2):

        binRanges=[(0,-1)]
        if i==1:
            binRanges=[(0,gen.GetXaxis().FindBin(50)),
                       (gen.GetXaxis().FindBin(50)+1,gen.GetXaxis().FindBin(300)),
                       (gen.GetXaxis().FindBin(300)+1,-1)]
            
        for binmin,binmax in binRanges:

            xmin,xmax=gen.GetXaxis().GetBinLowEdge(binmin),gen.GetXaxis().GetBinUpEdge(binmax)
            if binmin==0 : xmin=gen.GetXaxis().GetXmin()
            if binmax==-1: xmax=gen.GetXaxis().GetXmax()


            c.Clear()
            c.SetRightMargin(0.04)
            c.SetLogy()        
            genprof = gen.ProjectionX('genx',binmin,binmax) if i==0 else gen.ProjectionY('geny',binmin,binmax)
            fixExtremities(genprof)
            genprof.SetLineColor(1)
            genprof.SetFillStyle(1001)
            genprof.SetFillColor(ROOT.kCyan-10)
            genprof.SetTitle('gen')
            genprof.Draw('hist')
            genprof.GetYaxis().SetTitle('Jets')
            genprof.GetYaxis().SetTitleOffset(1.2)
            genprof.GetYaxis().SetLabelSize(0.04)
            genprof.GetYaxis().SetTitleSize(0.04)
            genprof.GetYaxis().SetRangeUser(1,genprof.GetMaximum()*2)
        
            #if i==1 and hasPU : genprof.GetYaxis().SetRangeUser(1,genprof.GetMaximum()*2)

            recprof = rec.ProjectionX('recx',binmin,binmax) if i==0 else rec.ProjectionY('recy',binmin,binmax)
            fixExtremities(recprof)
            recprof.SetLineColor(1)
            recprof.SetLineWidth(2)
            recprof.SetTitle('reco (matched)')
            recprof.Draw('histsame')
        
            recunmprof = rec_unm.ProjectionX('recunmx',binmin,binmax) if i==0 else rec_unm.ProjectionY('recunmy',binmin,binmax)
            fixExtremities(recunmprof)
            recunmprof.SetLineColor(ROOT.kRed)
            recunmprof.SetLineWidth(2)
            recunmprof.SetTitle('reco (unmatched)')
            recunmprof.Draw('histsame')
            
            leg=c.BuildLegend(0.6,0.8,0.95,0.95)
            leg.SetFillStyle(0)
            leg.SetBorderSize(0)
            leg.SetTextFont(42)
            leg.SetTextSize(0.035)
            
            if binmin!=0 or binmax!=-1:
                if binmax==-1 : MyPaveText('[ >%3.2f ]'%(xmin),0.2,0.85,0.4,0.9).SetTextSize(0.035)
                else: MyPaveText('[ %3.2f-%3.2f ]'%(xmin,xmax),0.2,0.85,0.4,0.9).SetTextSize(0.035)
            else:
                MyPaveText('[ inclusive ]',0.2,0.85,0.4,0.9).SetTextSize(0.035)
            MyPaveText('#bf{CMS} #it{simulation}    %s jets'%algo)

            c.Modified()
            c.Update()
            c.SaveAs('%s/genkin_%d_%d_%d_%s.png'%(outputDir,i,binmin,binmax,algo))
    
"""
"""
def showConstituentsProfiles(algo,fIn,outputDir) :

    varMap={
        'enfrac':'<Energy fraction>',
        'en':'<Total energy> [GeV]',
        'inden':'<Energy> [GeV]',
        'mult':'<Multiplicity>'
    }
    
    c=ROOT.TCanvas('c','c',500,500)
    for jetType in ['','_unm']:
        for profvar in ['pt','eta']:

            for var in varMap:

                binCats=['']
                binCatTitles={'':'[ inclusive ]'}
                if profvar=='eta' : 
                    binCats+=['_pt30','_pt70','_pt200']
                    binCatTitles['_pt30']='[ 30<p_{T}<70 ]'
                    binCatTitles['_pt70']='[ 70<p_{T}<20 ]'
                    binCatTitles['_pt200']='[ p_{T}>200 ]'
                    
                for binCat in binCats:
                    
                    ch2d=fIn.Get('analysis_%s/chf_%s%s_%s%s'%(algo,profvar,binCat,var,jetType))
                    xmin=ch2d.GetXaxis().GetXmin()
                    xmin=ch2d.GetXaxis().GetXmax()

                    ch=ch2d.ProfileX('ch_prof')
                    ch.SetMarkerStyle(20)
                    ch.SetTitle('charged had.')

                    gamma=fIn.Get('analysis_%s/gamma_%s%s_%s%s'%(algo,profvar,binCat,var,jetType)).ProfileX('gamma_prof')
                    gamma.SetMarkerStyle(24)
                    gamma.SetMarkerColor(ROOT.kRed)
                    gamma.SetLineColor(ROOT.kRed)
                    gamma.SetTitle('photons')
                    
                    nh=fIn.Get('analysis_%s/nhf_%s%s_%s%s'%(algo,profvar,binCat,var,jetType)).ProfileX('nh_prof')
                    nh.SetMarkerStyle(23)
                    nh.SetMarkerColor(30)
                    nh.SetLineColor(30)
                    nh.SetTitle('neutral had.')

                    plotTag='const%s_%s_%s%s_%s'%(jetType,var,profvar,binCat,algo)

                    iplot=len(ALLPLOTS)
                    ALLPLOTS.append((plotTag,[ch.Clone(),gamma.Clone(),nh.Clone()]))
                    for plot in ALLPLOTS[iplot][1]: plot.SetDirectory(0)

                    c.Clear()
                    c.SetRightMargin(0.04)
                    ch.Draw('e1')
                    ch.GetYaxis().SetLabelSize(0.04)
                    ch.GetYaxis().SetTitleSize(0.04)
                    ch.GetYaxis().SetTitleOffset(1.4)
                    ch.GetYaxis().SetTitle(varMap[var])
                    ch.GetYaxis().SetRangeUser(ch2d.GetYaxis().GetXmin(),ch2d.GetYaxis().GetXmax())
                    gamma.Draw('e1same')
                    nh.Draw('e1same')
                    
                    leg=c.BuildLegend(0.6,0.8,0.95,0.95)
                    leg.SetFillStyle(0)
                    leg.SetBorderSize(0)
                    leg.SetTextFont(42)
                    leg.SetTextSize(0.035)

                    jetTitle=algo
                    if jetType!='' : jetTitle='unmatched %s'%algo
                    MyPaveText('#bf{CMS} #it{simulation}    %s jets'%jetTitle)
                    MyPaveText('%s'%binCatTitles[binCat],0.2,0.85,0.4,0.9).SetTextSize(0.035)
                    c.Modified()
                    c.Update()
                    c.SaveAs('%s/const%s_%s_%s%s_%s.png'%(outputDir,jetType,var,profvar,binCat,algo))


"""
"""
def showResponseSummary(algo,fIn,outputDir,hasPU) :

    ROOT.gStyle.SetOptFit(1111)
    
    for profvar in ['pt','eta']:

        step=3
        varTitle="p_{T}"
        if profvar=="eta" : 
            varTitle="#eta"
            step=1

        respGr=[]
        resGr=[]
        binCats=['']
        binCatTitles={'':'inclusive'}
        if profvar=='eta' : 
            binCats+=['_pt30','_pt70','_pt200']
            binCatTitles['_pt30']='30<p_{T}<70'
            binCatTitles['_pt70']='70<p_{T}<20'
            binCatTitles['_pt200']='p_{T}>200'
                
        binCatCtr=-1
        for binCat in binCats:
            binCatCtr+=1
            ptresp=fIn.Get('analysis_%s/ptresp_%s%s'%(algo,profvar,binCat))
        
            profs=[]
            respGr.append( ROOT.TGraphErrors() )
            respGr[binCatCtr].SetName('resp_%s%s_%s'%(profvar,binCat,algo))
            respGr[binCatCtr].SetMarkerStyle(19+binCatCtr)
            respGr[binCatCtr].SetMarkerSize(2.5)
            respGr[binCatCtr].SetMarkerColor(1+binCatCtr)
            respGr[binCatCtr].SetLineColor(1+binCatCtr)
            respGr[binCatCtr].SetFillColor(ROOT.kCyan-10)
            respGr[binCatCtr].SetFillStyle(1001)
            respGr[binCatCtr].SetTitle(binCatTitles[binCat])
            resGr.append( respGr[binCatCtr].Clone('res_%s%s_%s'%(profvar,binCat,algo)) )

            for xbin in xrange(1,ptresp.GetXaxis().GetNbins(),step):

                np=len(profs)
                xbinini=xbin
                xbinfin=xbin+(step-1)
                xmin,xmax=ptresp.GetXaxis().GetBinCenter(xbinini),ptresp.GetXaxis().GetBinCenter(xbinfin)
                x=0.5*(xmin+xmax)
                xerr=(xmax-xmin)*0.5
                profs.append( ptresp.ProjectionY('prof_%d_%d_%d_%d'%(xbin,xbinini,xbinfin,binCatCtr),xbinini,xbinfin ) )

                #require at least 5 jets for a gaussian fit
                if profs[np].Integral()<5 :
                    profs.pop()
                    continue

                fixExtremities(profs[np])
                profTitle='%3.1f < %s < %3.1f'%(xmin,varTitle,xmax)
                if binCat!='' : profTitle += ', %s'%binCatTitles[binCat]
                profs[np].SetTitle(profTitle)

                maxBin=profs[np].GetMaximumBin()
                centerX=profs[np].GetXaxis().GetBinCenter(maxBin)
                if centerX>3 or centerX<0.5 : continue
            
                profs[np].Fit('gaus','LMRQ+','',centerX-0.4, centerX+0.4)
                gaus=profs[np].GetFunction('gaus')
                mean, mean_err = gaus.GetParameter(1), gaus.GetParError(1)
                sigma, sigma_err = gaus.GetParameter(2), gaus.GetParError(2)

                respGr[binCatCtr].SetPoint(np,x,mean)
                respGr[binCatCtr].SetPointError(np,xerr,mean_err)

                resGr[binCatCtr].SetPoint(np,x,sigma/mean)
                resGr[binCatCtr].SetPointError(np,xerr,ROOT.TMath.Sqrt(ROOT.TMath.Power(sigma*mean_err,2)+ROOT.TMath.Power(sigma_err*mean,2))/(mean*mean))

            if len(profs)==0 : continue

            ndiv=int(math.floor(ROOT.TMath.Sqrt(len(profs))))
            while ndiv*ndiv < len(profs) : ndiv+=1
            ww=int(ndiv*300)
            c=ROOT.TCanvas('c','c',ww,ww)
            c.Divide(ndiv,ndiv)
            for np in xrange(0,len(profs)):
                c.cd(np+1)
                profs[np].Draw('e1')
                profs[np].GetYaxis().SetTitle('Jets')
                profs[np].GetXaxis().SetTitleSize(0.05)
                profs[np].GetYaxis().SetTitleSize(0.05)
                profs[np].GetXaxis().SetLabelSize(0.05)
                profs[np].GetYaxis().SetLabelSize(0.05)
                MyPaveText('[ %s %s ]'%(profs[np].GetTitle(),binCatTitles[binCat]),0.45,0.5,0.9,0.6).SetTextSize(0.05)
            if np==0:
                MyPaveText('#bf{CMS} #it{simulation}    %s jets'%algo).SetTextSize(0.055)

            c.Modified()
            c.Update()
            c.SaveAs('%s/jetresponse_%s%s_%s.png'%(outputDir,profvar,binCat,algo))
        
        c.SetWindowSize(500,500)
        grctr=0
        for grColl,grtitle in [ (respGr,'<p_{T}(reco)/p_{T}(gen)>'), (resGr,'#sigma_{p_{T}(reco)/p_{T}(gen)}>/<p_{T}(reco)/p_{T}(gen)>') ]:
            grctr+=1
            c.Clear()
            c.SetLeftMargin(0.15)

            leg=ROOT.TLegend(0.6,0.75,0.95,0.95)
            leg.SetFillStyle(0)
            leg.SetBorderSize(0)
            leg.SetTextFont(42)
            leg.SetTextSize(0.035)

            igr=0
            for gr in grColl:
                if gr.GetN()==0 : continue
                leg.AddEntry(gr,gr.GetTitle(),'p')
                if igr==0 : 
                    if len(grColl)==1: gr.Draw('ap')
                    else : gr.Draw('a3')
                    gr.GetXaxis().SetTitle(varTitle)
                    gr.GetYaxis().SetTitle(grtitle)
                    gr.GetXaxis().SetTitleSize(0.05)
                    gr.GetYaxis().SetTitleSize(0.05)
                    gr.GetYaxis().SetTitleOffset(1.4)
                    gr.GetXaxis().SetLabelSize(0.04)
                    gr.GetYaxis().SetLabelSize(0.04)
                    if profvar=='eta': gr.GetXaxis().SetRangeUser(1.5,3.0)
                    if grctr==1 : 
                        gr.GetYaxis().SetRangeUser(0.5,2)
                    else        : 
                        gr.GetYaxis().SetRangeUser(0,0.3)
                else:
                    gr.Draw('p')
                igr+=1
            plotTag='jetresponseevol_%s_%d_%s.png'%(profvar,grctr,algo)
            ALLPLOTS.append((plotTag,grColl))

            leg.Draw()
            MyPaveText('#bf{CMS} #it{simulation}    %s jets'%algo)
            c.Modified()
            c.Update()
            c.SaveAs('%s/jetresponseevol_%s_%d_%s.png'%(outputDir,profvar,grctr,algo))


"""
steer 
"""
def main():
    
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i',      '--in' ,      dest='input',    help='Input file',         default=None)
    parser.add_option('-o',      '--out' ,     dest='output',   help='Output directory',   default='plots')
    (opt, args) = parser.parse_args()

    #check inputs    
    if opt.input is None :
        parser.print_help()
        sys.exit(1)

    header=os.path.splitext( os.path.basename(opt.input) )[0]
    outputDir=opt.output+'/'+header
    os.system('mkdir -p %s'%outputDir)
    cmsswBase=os.environ['CMSSW_BASE']
    os.system('sed s/HEADERGOESHERE/%s/ %s/src/UserCode/HGCanalysis/test/analysis/jetanalysis_index.html > %s/index.html'%(header,cmsswBase,outputDir))

    print '[showJetAnalysisPlots] results from %s will be available in %s/index.html'%(opt.input,outputDir)

    customROOTstyle()
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    hasPU=False
    if outputDir.find('PU140') : hasPU=True

    fIn=ROOT.TFile.Open(opt.input)
    #greedy scan...
    for algo in ['kt4','ak2','ak4','ak5','ak4_cmspf']:
        showEfficiencyMap(algo,fIn,outputDir,hasPU)
        showConstituentsProfiles(algo,fIn,outputDir)
        showResponseSummary(algo,fIn,outputDir,hasPU)
    fIn.Close()


    fOut=ROOT.TFile('%s/plots.root'%outputDir,'RECREATE')
    for plotTag,plotColl in ALLPLOTS:
        fOut.cd()
        fOut.mkdir(plotTag)
        fOut.cd(plotTag)
        for plot in plotColl: plot.Write()
    fOut.Close()

    

    
if __name__ == "__main__":
    sys.exit(main())


