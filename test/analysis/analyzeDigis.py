#!/usr/bin/env python

import ROOT
import io,os,sys
import optparse
import commands
from array import array

from UserCode.HGCanalysis.PlotUtils import *

def analyzeDigis(opt) :

	#prepare output
	os.system('mkdir -p %s'%opt.output)

	fIn=ROOT.TFile.Open(opt.input)

	hnames={}
	hnames['mip']=['Energy [MIP]',         0,50]
	#hnames['tot']=['Integration time [bx]',1,5]
	#hnames['busy']=['Busy state',1,2]

	c=ROOT.TCanvas('c','c',1200,800)
	allProjs=[]
	allLegs=[]
	#quantileName='99% quantile'
	#probSum = array('d', [0.99])
	#quantileName='90% quantile'
	#probSum = array('d', [0.90])
	quantileName='Median'
	probSum = array('d', [0.5])
	quantiles = array('d', [0.0])
	for hname in hnames:
		xtitle=hnames[hname][0]
		xmin=hnames[hname][1]
		xmax=hnames[hname][2]
		profGr={}
		for layer in xrange(1,31):
			hch=fIn.Get("analysis/sd_0_layer%d_nch"%layer)
			h=fIn.Get("analysis/sd_0_layer%d_%s"%(layer,hname))

			#prepare canvas
			c.Clear()
			netaBins=h.GetYaxis().GetNbins()
			divx=int(netaBins/3)
			if divx*3 < netaBins : divx+=1
			c.Divide(divx,2)
			allPads=[]
			for ybin in xrange(1,h.GetYaxis().GetNbins()):
				nch=2*hch.GetBinContent(ybin)
				if not (ybin in profGr): profGr[ybin]=[]

				allPads.append( c.cd(ybin) )				
				allPads[ybin-1].SetLogy()
				etaMin=h.GetYaxis().GetBinLowEdge(ybin)
				etaMax=h.GetYaxis().GetBinUpEdge(ybin)
				if ybin==1 :
					ileg = len(allLegs)
					allLegs.append( ROOT.TLegend(0.6,0.6,0.9,0.9) )
					allLegs[ileg].SetTextSize(0.05)
					allLegs[ileg].SetBorderSize(0)
					allLegs[ileg].SetFillStyle(0)
					allLegs[ileg].SetTextFont(42)
					allLegs[ileg].SetTextSize(0.05)
				for zbin in xrange(1,h.GetZaxis().GetNbins()+1):
					drawOpt=''
					if zbin>1 : drawOpt+='same'
					iproj=len(allProjs)
					allProjs.append( h.ProjectionX('proj%d%d'%(ybin,zbin),ybin,ybin,zbin,zbin) )
					fixExtremities(allProjs[iproj])
					
					#profile plots
					if len(profGr[ybin])<zbin:
						profGr[ybin].append( ROOT.TGraphErrors() )
						profGr[ybin][zbin-1].SetMarkerStyle(20+zbin)
						profGr[ybin][zbin-1].SetMarkerColor(42-zbin/2+1)
						profGr[ybin][zbin-1].SetLineColor(42-zbin/2+1)
						profGr[ybin][zbin-1].SetLineWidth(zbin%2+1)	
						profGr[ybin][zbin-1].SetTitle('bx=-%d'%(h.GetZaxis().GetNbins()-zbin))
					allProjs[iproj].GetQuantiles(1,quantiles,probSum)
					np=profGr[ybin][zbin-1].GetN()
					offset=allProjs[iproj].GetXaxis().GetBinCenter(1)
					profGr[ybin][zbin-1].SetPoint(np,layer,quantiles[0]-offset)
					profGr[ybin][zbin-1].SetPointError(np,0,1.253*allProjs[iproj].GetMeanError())

					allProjs[iproj].Scale(1./nch)
					allProjs[iproj].SetDirectory(0)
					allProjs[iproj].Draw(drawOpt)
					allProjs[iproj].SetMarkerStyle(20+zbin)
					allProjs[iproj].SetMarkerColor(42-zbin/2+1)
					allProjs[iproj].SetLineColor(42-zbin/2+1)
					allProjs[iproj].SetLineWidth(zbin%2+1)	
					allProjs[iproj].SetTitle('bx=-%d'%(h.GetZaxis().GetNbins()-zbin))
					allProjs[iproj].GetYaxis().SetTitle('Fraction of cells / Event')
					allProjs[iproj].GetXaxis().SetTitle(xtitle)
					if ybin==1 : allLegs[ len(allLegs)-1 ].AddEntry( allProjs[iproj], allProjs[iproj].GetTitle(), 'lp')
					if zbin==1 :
						allProjs[iproj].GetXaxis().SetRangeUser(xmin,xmax)
						allProjs[iproj].GetYaxis().SetTitleSize(0.06)
						allProjs[iproj].GetXaxis().SetTitleSize(0.06)
						allProjs[iproj].GetYaxis().SetTitleOffset(0.9)
						allProjs[iproj].GetXaxis().SetTitleOffset(0.9)
						allProjs[iproj].GetYaxis().SetLabelSize(0.05)
						allProjs[iproj].GetXaxis().SetLabelSize(0.05)
						allProjs[iproj].GetYaxis().SetRangeUser(1e-6,2)
				MyPaveText('[%3.2f<|#eta|<%3.2f]'%(etaMin,etaMax),0.6,0.9,0.9,0.95).SetTextSize(0.05)
				if ybin==1 : 
					MyPaveText('#bf{CMS} #it{simulation}').SetTextSize(0.05)
					MyPaveText('#it{Layer %d}'%layer,0.75,0.95,0.95,0.99).SetTextSize(0.05)
					allLegs[ybin-1].Draw()	

			c.Modified()
			c.Update()
			c.SaveAs('%s/%s_layer%d.png'%(opt.output,hname,layer))
			c.SaveAs('%s/%s_layer%d.C'%(opt.output,hname,layer))

		#now show profiles
		for ybin in profGr:
			pad=c.cd(ybin)
			pad.Clear()
			pad.SetLogy(False)
			etaMin=h.GetYaxis().GetBinLowEdge(ybin)
			etaMax=h.GetYaxis().GetBinUpEdge(ybin)
			if ybin==1:
				ileg = len(allLegs)
				allLegs.append( ROOT.TLegend(0.6,0.6,0.9,0.9) )
				allLegs[ileg].SetTextSize(0.05)
				allLegs[ileg].SetBorderSize(0)
				allLegs[ileg].SetFillStyle(0)
				allLegs[ileg].SetTextFont(42)
				allLegs[ileg].SetTextSize(0.05)
			for igr in xrange(0,len(profGr[ybin])):
				if ybin==1:
					allLegs[ileg].AddEntry(profGr[ybin][igr],profGr[ybin][igr].GetTitle(),'p')
				if igr==0:
					profGr[ybin][igr].Draw('ap')
					profGr[ybin][igr].GetXaxis().SetTitle('Layer #')
					profGr[ybin][igr].GetYaxis().SetTitle('%s %s'%(quantileName,xtitle))
					profGr[ybin][igr].GetYaxis().SetRangeUser(profGr[ybin][igr].GetYaxis().GetXmin(),profGr[ybin][igr].GetYaxis().GetXmax()*1.5)
				else:
					profGr[ybin][igr].Draw('p')
			MyPaveText('[%3.2f<|#eta|<%3.2f]'%(etaMin,etaMax),0.6,0.9,0.9,0.95).SetTextSize(0.05)
			if ybin==1:
				allLegs[ileg].Draw()
				MyPaveText('#bf{CMS} #it{simulation}').SetTextSize(0.05)

		c.Modified()
		c.Update()
		c.SaveAs('%s/%s_profile.png'%(opt.output,hname))
		c.SaveAs('%s/%s_profile.C'%(opt.output,hname))



"""
steer 
"""
def main():
    
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i',      '--in' ,      dest='input',        help='Input file',   default=None)
    parser.add_option('-o',      '--out' ,     dest='output',       help='Output dir',   default=None)
    (opt, args) = parser.parse_args()

     #check inputs                                                                                                                                                                                                  
    if opt.input is None and opt.wsUrl is None:
        parser.print_help()
        sys.exit(1)

    #basic ROOT customization
    customROOTstyle()
    #ROOT.gROOT.SetBatch(False)
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    
    analyzeDigis(opt)
    
if __name__ == "__main__":
    sys.exit(main())


