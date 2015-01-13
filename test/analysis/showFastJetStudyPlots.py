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

c.Clear()
gr={}
gjcount=fin.Get('gjcount')
jcount=fin.Get('jcount')
gjcount_notkint=fin.Get('gjcount_notkint')
jcount_notkint=fin.Get('jcount_notkint')
if gjcount and jcount:
    gr['inc']=ROOT.TGraphAsymmErrors()
    gr['inc'].BayesDivide(jcount,gjcount)
    gr['inc'].SetMarkerStyle(20)
    gr['inc'].SetFillStyle(0)
    gr['inc'].SetTitle('inclusive')
    gr['inc'].Draw('ap')
    gr['inc'].GetYaxis().SetTitle("Efficiency")
    gr['inc'].GetXaxis().SetTitle("Pseudo-rapidity")
    gr['inc'].GetYaxis().SetRangeUser(0.7,1.05)
    gr['notkint']=ROOT.TGraphAsymmErrors()
    gr['notkint'].BayesDivide(jcount_notkint,gjcount_notkint)
    gr['notkint'].SetMarkerStyle(20)
    gr['notkint'].SetMarkerColor(ROOT.kRed)
    gr['notkint'].SetLineColor(ROOT.kRed)
    gr['notkint'].SetFillStyle(0)
    gr['notkint'].SetTitle('interacting in HGC')
    gr['notkint'].Draw('p')
    leg=c.BuildLegend(0.3,0.3,0.5,0.5)    
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.03)
    MyPaveText('#bf{CMS} #it{simulation}')
    c.SaveAs('fastjet_eff.png')

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
