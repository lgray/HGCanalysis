#! /usr/bin/env python

import ROOT
import sys
from UserCode.HGCanalysis.PlotUtils import *

customROOTstyle()
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)


fin=ROOT.TFile.Open(sys.argv[1])
c=ROOT.TCanvas('c','c',500,500)
for h in ['ptresol','eresol','ptresp','eresp']:
    histo=fin.Get(h)
    histonotkint=fin.Get(h+'_notkint')
    c.Clear()
    histo.Draw('hist')
    histo.SetTitle('inclusive')
    if histonotkint:
        histonotkint.SetTitle('interacting in HGC')
        histonotkint.SetLineColor(ROOT.kRed)
        histonotkint.Draw('histsame')
        leg=c.BuildLegend(0.2,0.8,0.4,0.95)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.03)
    MyPaveText('#bf{CMS} #it{simulation}')
    c.SaveAs('fastjet_%s.png'%h)

gr={}
for cat in ['','_notkint']:
    c.Clear()
    fracplots=[]
    for name in ['nhfrac','mufrac','efrac','chfrac','gfrac']:
        np=len(fracplots)
        h=fin.Get(name+cat) 
        if not h: 
            continue
        fracplots.append( h )
        fracplots[np].SetMarkerStyle(20+np)
        color=np+1
        if color>=5: color+=1
        fracplots[np].SetLineColor(color)
        fracplots[np].SetMarkerColor(color)
        fracplots[np].SetLineWidth(2)
        fracplots[np].GetYaxis().SetRangeUser(1,fracplots[np].GetMaximum()*1.2)
        if np==0 : fracplots[np].Draw('e1')
        else : fracplots[np].Draw('e1same')
    c.SetLogy()
    MyPaveText('#bf{CMS} #it{simulation}')
    leg=c.BuildLegend(0.3,0.8,0.7,0.95)    
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.03)
    leg.SetNColumns(2)
    c.SaveAs('fastjet_constfrac%s.png'%cat)

    c.Clear()
    c.SetLogy(False)
    gjcount=fin.Get('gjcount'+cat)
    jcount=fin.Get('jcount'+cat)
    if not gjcount or not jcount: continue
    gr[cat]=ROOT.TGraphAsymmErrors()
    gr[cat].BayesDivide(jcount,gjcount)
    gr[cat].SetMarkerStyle(20)
    gr[cat].Draw('ap')
    gr[cat].GetYaxis().SetTitle("Efficiency")
    gr[cat].GetXaxis().SetTitle("Pseudo-rapidity")
    gr[cat].GetYaxis().SetRangeUser(0.7,1.05)
    MyPaveText('#bf{CMS} #it{simulation}')
    c.SaveAs('fastjet_eff%s.png'%cat)

    c.Clear()
    c.SetRightMargin(0.12)
    h2d=fin.Get('erespvsgfrac'+cat)
    if not h2d : continue
    h2d.Draw('colz')
    h2d.GetZaxis().SetTitleOffset(-0.3)
    h2dprof=h2d.ProfileX('erespvsgfracprof')
    h2dprof.SetMarkerStyle(20)
    h2dprof.Draw('e1same')
    MyPaveText('#bf{CMS} #it{simulation}')
    c.SaveAs('fastjet_erespvsgfrac%s.png'%cat)

fin.Close()
