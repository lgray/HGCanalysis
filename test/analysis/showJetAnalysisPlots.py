#!/usr/bin/env python

import ROOT
import sys
import math
import optparse
import os
from UserCode.HGCanalysis.PlotUtils import *

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
        c.Clear()
        c.SetRightMargin(0.04)
        c.SetLogy()
        
        genprof = gen.ProjectionX('genx') if i==0 else gen.ProjectionY('geny')
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
        
        if i==1 and hasPU : genprof.GetYaxis().SetRangeUser(1,genprof.GetMaximum()*4)

        recprof = rec.ProjectionX('recx') if i==0 else rec.ProjectionY('recy')
        fixExtremities(recprof)
        recprof.SetLineColor(1)
        recprof.SetLineWidth(2)
        recprof.SetTitle('reco (matched)')
        recprof.Draw('histsame')
        
        recunmprof = rec_unm.ProjectionX('recunmx') if i==0 else rec_unm.ProjectionY('recunmy')
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

        MyPaveText('#bf{CMS} #it{simulation}    %s jets'%algo)
        c.Modified()
        c.Update()
        c.SaveAs('%s/genkin_%d_%s.png'%(outputDir,i,algo))
    
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
                ch2d=fIn.Get('analysis_%s/chf_%s_%s%s'%(algo,profvar,var,jetType))

                ch=ch2d.ProfileX('ch_prof')
                ch.SetMarkerStyle(20)
                ch.SetTitle('charged had.')

                gamma=fIn.Get('analysis_%s/gamma_%s_%s%s'%(algo,profvar,var,jetType)).ProfileX('gamma_prof')
                gamma.SetMarkerStyle(24)
                gamma.SetMarkerColor(ROOT.kRed)
                gamma.SetLineColor(ROOT.kRed)
                gamma.SetTitle('photons')
                
                nh=fIn.Get('analysis_%s/nhf_%s_%s%s'%(algo,profvar,var,jetType)).ProfileX('nh_prof')
                nh.SetMarkerStyle(23)
                nh.SetMarkerColor(30)
                nh.SetLineColor(30)
                nh.SetTitle('neutral had.')
                
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
                c.Modified()
                c.Update()
                c.SaveAs('%s/const%s_%s_%s_%s.png'%(outputDir,jetType,var,profvar,algo))


"""
"""
def showResponseSummary(algo,fIn,outputDir,hasPU) :

    ROOT.gStyle.SetOptFit(1111)
    
    respGr={}
    resGr={}
    for profvar in ['pt','eta']:

        varTitle="p_{T}"
        if profvar=="eta" : varTitle="#eta"

        ptresp=fIn.Get('analysis_%s/ptresp_%s'%(algo,profvar))

        profs=[]
        respGr[profvar]=ROOT.TGraphErrors()
        respGr[profvar].SetName('resp_%s_%s'%(profvar,algo))
        respGr[profvar].SetMarkerStyle(20)
        respGr[profvar].SetTitle(algo)
        resGr[profvar]=respGr[profvar].Clone('res_%s_%s'%(profvar,algo))
        for xbin in xrange(1,ptresp.GetXaxis().GetNbins()):
            x=ptresp.GetXaxis().GetBinCenter(xbin)
            xerr=ptresp.GetXaxis().GetBinWidth(xbin)
            np=len(profs)
            profs.append( ptresp.ProjectionY('prof_%d'%xbin,xbin,xbin) )

            #require at least 20 jets for a gaussian fit
            if profs[np].Integral()<20 :
                profs.pop()
                continue

            xmin=ptresp.GetXaxis().GetBinLowEdge(xbin)
            xmax=ptresp.GetXaxis().GetBinUpEdge(xbin)

            fixExtremities(profs[np])
            profs[np].SetTitle('%3.1f < %s < %3.1f'%(xmin,varTitle,xmax))
            
            maxBin=profs[np].GetMaximumBin()
            centerX=profs[np].GetXaxis().GetBinCenter(maxBin)
            if centerX>3 or centerX<0.5 : continue
            
            profs[np].Fit('gaus','LMRQ+','',centerX-0.4, centerX+0.4)
            gaus=profs[np].GetFunction('gaus')
            mean, mean_err = gaus.GetParameter(1), gaus.GetParError(1)
            sigma, sigma_err = gaus.GetParameter(2), gaus.GetParError(2)

            respGr[profvar].SetPoint(np,x,mean)
            respGr[profvar].SetPointError(np,xerr*0.5,mean_err)

            resGr[profvar].SetPoint(np,x,sigma/mean)
            resGr[profvar].SetPointError(np,xerr*0.5,ROOT.TMath.Sqrt(ROOT.TMath.Power(sigma*mean_err,2)+ROOT.TMath.Power(sigma_err*mean,2))/(mean*mean))


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
            MyPaveText('[ %s ]'%profs[np].GetTitle(),0.6,0.5,0.9,0.6).SetTextSize(0.05)
            if np==0:
                MyPaveText('#bf{CMS} #it{simulation}    %s jets'%algo).SetTextSize(0.055)

        c.Modified()
        c.Update()
        c.SaveAs('%s/jetresponse_%s_%s.png'%(outputDir,profvar,algo))


        c.SetWindowSize(500,500)
        grctr=0
        for gr,grtitle in [ (respGr[profvar],'<p_{T}(reco)/p_{T}(gen)>'), (resGr[profvar],'#sigma_{p_{T}(reco)/p_{T}(gen)}>/<p_{T}(reco)/p_{T}(gen)>') ]:
            grctr+=1
            c.Clear()
            c.SetLeftMargin(0.15)
            gr.Draw('ap')
            gr.SetMarkerStyle(20)
            gr.GetXaxis().SetTitle(varTitle)
            gr.GetYaxis().SetTitle(grtitle)
            gr.GetXaxis().SetTitleSize(0.05)
            gr.GetYaxis().SetTitleSize(0.05)
            gr.GetYaxis().SetTitleOffset(1.4)
            gr.GetXaxis().SetLabelSize(0.04)
            gr.GetYaxis().SetLabelSize(0.04)
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
    for algo in ['kt4','ak2','ak3','ak4','ak5','ak4_cmspf']:
        showEfficiencyMap(algo,fIn,outputDir,hasPU)
        showConstituentsProfiles(algo,fIn,outputDir)
        showResponseSummary(algo,fIn,outputDir,hasPU)

    
if __name__ == "__main__":
    sys.exit(main())


