#! /usr/bin/env python

import ROOT
import sys
from UserCode.HGCanalysis.PlotUtils import *

customROOTstyle()
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)


fin=ROOT.TFile.Open('FastJetResol.root')
c=ROOT.TCanvas('c','c',500,500)
for h in ['ptresol','eresol','ptresp','eresp']:
    histo=fin.Get(h)
    c.Clear()
    histo.Draw('hist')
    MyPaveText('#bf{CMS} #it{simulation}')
    c.SaveAs('fastjet_%s.png'%h)

c.Clear()
fracplots=[]
for name in ['nhfrac','mufrac','efrac','chfrac','gfrac']:
    np=len(fracplots)
    fracplots.append( fin.Get(name) )
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
leg=c.BuildLegend()
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextFont(42)
c.SaveAs('fastjet_constfrac.png')

c.Clear()
c.SetLogy(False)
gjcount=fin.Get('gjcount')
jcount=fin.Get('jcount')
gr=ROOT.TGraphAsymmErrors()
gr.BayesDivide(jcount,gjcount)
gr.SetMarkerStyle(20)
gr.Draw('ap')
gr.GetYaxis().SetTitle("Efficiency")
gr.GetXaxis().SetTitle("Pseudo-rapidity")
gr.GetYaxis().SetRangeUser(0.7,1.05)
MyPaveText('#bf{CMS} #it{simulation}')
c.SaveAs('fastjet_eff.png')

c.Clear()
c.SetRightMargin(0.12)
h2d=fin.Get('erespvsgfrac')
h2d.Draw('colz')
h2d.GetZaxis().SetTitleOffset(-0.3)
h2dprof=h2d.ProfileX('erespvsgfracprof')
h2dprof.SetMarkerStyle(20)
h2dprof.Draw('e1same')
MyPaveText('#bf{CMS} #it{simulation}')
c.SaveAs('fastjet_erespvsgfrac.png')

fin.Close()
