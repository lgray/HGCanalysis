import ROOT

from UserCode.HGCanalysis.PlotUtils import *

"""
"""
def showResults(gr,grecalhb,ymin,ymax,ytitle,tag):
        #show results
	customROOTstyle()
	ROOT.gROOT.SetBatch(False)

	canvas=ROOT.TCanvas('c','c',500,500)
	canvas.SetRightMargin(0.05)
	canvas.SetTopMargin(0.05)
	canvas.SetLeftMargin(0.15)

	grframe=ROOT.TGraph()
	grframe.SetPoint(0,0,0)
	grframe.SetPoint(1,1000,2)
	grframe.SetMarkerStyle(1)

	allLegends=[]
	for key in gr:
		canvas.Clear()
		canvas.SetGridx()
		canvas.SetGridy()
		canvas.SetLogx()
		grframe.Draw('ap')
		grframe.GetXaxis().SetTitle('Energy [GeV]')
		grframe.GetYaxis().SetTitle(ytitle)
		grframe.GetYaxis().SetTitleOffset(1.3)
		grframe.GetYaxis().SetRangeUser(ymin,ymax)
		grframe.GetXaxis().SetRangeUser(5,300)

		np=len(allLegends)
		allLegends.append(ROOT.TLegend(0.3,0.7,0.9,0.94))
		allLegends[np].SetFillStyle(1001)
		allLegends[np].SetFillColor(0)
		allLegends[np].SetTextFont(42)
		allLegends[np].SetBorderSize(0)
		allLegends[np].SetTextSize(0.028)

		if grecalhb is None:
			print 'Nothing to overlay'
		else:
			grecalhb.Draw('3p')
			allLegends[np].AddEntry(grecalhb,grecalhb.GetTitle(),'f')

		for ip in xrange(0,len(gr[key])): 
			gr[key][ip].Draw('p')
			allLegends[np].AddEntry(gr[key][ip],gr[key][ip].GetTitle(),'p')

		allLegends[np].Draw()
	
		MyPaveText('#bf{CMS} #it{preliminary} [%s]'%titles[key]).SetTextSize(0.035)

		canvas.Modified()
		canvas.Update()
		canvas.SaveAs('%s_%s.png'%(tag,key))


"""
Pion resolution for Test beam
"""
def getTBresults():
	grecalhb=ROOT.TGraphErrors()
	grecalhb.SetName('ecalhb')
	grecalhb.SetMarkerColor(ROOT.kGray)
	grecalhb.SetMarkerStyle(0)
	grecalhb.SetFillColor(ROOT.kGray)
	grecalhb.SetFillStyle(1001)
	grecalhb.SetLineColor(ROOT.kGray)
	grecalhb.SetTitle('CMS EB+HB testbeam - EPJC 60 (2009)')
	
	for val in [[0.057077912884806, 0.08013245033112593],[0.07130452564848754, 0.09072847682119206],[0.0817369555300373, 0.0986754966887417],[0.10070157776344564, 0.11456953642384105],[0.1414748863475484, 0.14900662251655628],[0.1831888754306209, 0.18741721854304633],[0.22116216513819198, 0.20066225165562912],[0.354936999575278, 0.2814569536423841],[0.37761715247518524, 0.33311258278145695],[0.4079611772663636, 0.3582781456953642],[0.4477906592628714, 0.3900662251655629]] :
		np=grecalhb.GetN()
		grecalhb.SetPoint(np,ROOT.TMath.Power(1./val[0],2),val[1])
		grecalhb.SetPointError(np,0,0.05*val[1])
	return grecalhb


#photons
#resolModelType=0
#overlayTestBeam=False
#toPlot=[
#	['Single22_CMSSW_6_2_0_SLHC23_patch1_RECO-PU0_SimHits/edep_rec/calib__calib_uncalib.root',    'simple', 'Perfect cluster', True],
#	['Single22_CMSSW_6_2_0_SLHC23_patch1_RECO-PU0_SimHits/edep_clus/calib__calib_uncalib.root',   'simple', 'CMS PF',          False],
#	['Single22_CMSSW_6_2_0_SLHC23_patch1_pandoraRECO_SimHits/edep_clus/calib__calib_uncalib.root','simple', 'Pandora',         False],
#	['Single22_CMSSW_6_2_0_SLHC23_patch1_mqRECO_SimHits/edep_clus/calib__calib_uncalib.root',     'simple', 'Arbor',           False]
#]
#ymin,ymax=0,0.25

resolModelType=1
overlayTestBeam=True
#toPlot=[
#	['Single211_CMSSW_6_2_0_SLHC23_patch1_RECO-PU0_SimHits/edep_rec/calib__calib_uncalib.root',    'simple', 'Perfect cluster', True],
#	['Single211_CMSSW_6_2_0_SLHC23_patch1_RECO-PU0_SimHits/edep_clus/calib__calib_uncalib.root',   'simple', 'CMS PF',          False],
#	['Single211_CMSSW_6_2_0_SLHC23_patch1_pandoraRECO_SimHits/edep_clus/calib__calib_uncalib.root','simple', 'Pandora',         False],
#	['Single211_CMSSW_6_2_0_SLHC23_patch1_mqRECO_SimHits/edep_clus/calib__calib_uncalib.root',     'simple', 'Arbor',           False]
#	]

#toPlot=[
#	['Single211_CMSSW_6_2_0_SLHC23_patch1_RECO-PU0_SimHits/edep_rec/calib__calib_uncalib.root',    'simple', 'Perfect cluster', False],
#	['Single211_CMSSW_6_2_0_SLHC23_patch1_RECO-PU0_SimHits/edep_rec/calib__calib_uncalib.root',   'gc', 'Perfect cluster (GC)',          False],
#	['Single211_CMSSW_6_2_0_SLHC23_patch1_pandoraRECO_SimHits/edep_clus/calib__calib_uncalib.root','simple', 'Pandora',         False],
#	['Single211_CMSSW_6_2_0_SLHC23_patch1_pandoraRECO_SimHits/edep_clus/calib__calib_uncalib.root','gc', 'Pandora (GC)',         False]
#	]

#toPlot=[
#	['Single130-FixE_CMSSW_6_2_0_SLHC23_patch2_RECO-PU0_SimHits/edep_rec/calib__calib_uncalib.root',    'simple', 'Perfect cluster', True],
#	['Single130-FixE_CMSSW_6_2_0_SLHC23_patch2_RECO-PU0_SimHits/edep_clus/calib__calib_uncalib.root',   'simple', 'CMS PF',          False],
#	['Single130-FixE_CMSSW_6_2_0_SLHC23_patch2_pandoraRECO_SimHits/edep_clus/calib__calib_uncalib.root','simple', 'Pandora',         False],
#	['Single130-FixE_CMSSW_6_2_0_SLHC23_patch2_mqRECO_SimHits/edep_clus/calib__calib_uncalib.root',     'simple', 'Arbor',           False]
#	]

#toPlot=[
#	['Single130-FixE_CMSSW_6_2_0_SLHC23_patch2_RECO-PU0_SimHits/edep_rec/calib__calib_uncalib.root',    'simple', 'Perfect cluster', False],
#	['Single130-FixE_CMSSW_6_2_0_SLHC23_patch2_RECO-PU0_SimHits/edep_rec/calib__calib_uncalib.root',   'gc', 'Perfect cluster (GC)',          False],
#	['Single130-FixE_CMSSW_6_2_0_SLHC23_patch2_pandoraRECO_SimHits/edep_clus/calib__calib_uncalib.root','simple', 'Pandora',         False],
#	['Single130-FixE_CMSSW_6_2_0_SLHC23_patch2_pandoraRECO_SimHits/edep_clus/calib__calib_uncalib.root','gc', 'Pandora (GC)',         False]
#	]

#toPlot=[
#	['~/www/pion/rec/calib__calib_uncalib.root',           'simple', 'RecHit with scaling',  False],
#	['~/www/pion/rec-nocomp/calib__calib_uncalib.root',    'simple', 'RecHit',               False],
#	['~/www/pion/clus/calib__calib_uncalib.root',          'simple', 'Pandora with scaling', False],
#	['~/www/pion/clus-nocomp/calib__calib_uncalib.root',   'simple', 'Pandora',              False]
#	]

#toPlot=[
#	['~/www/pion/rec/calib__calib_uncalib.root',           'gc', 'RecHit with scaling + GC',  False],
#	['~/www/pion/rec-nocomp/calib__calib_uncalib.root',    'gc', 'RecHit+GC',               False],
#	['~/www/pion/clus/calib__calib_uncalib.root',          'gc', 'Pandora with scaling+GC', False],
#	['~/www/pion/clus-nocomp/calib__calib_uncalib.root',   'gc', 'Pandora+GC',              False]
#	]

toPlot=[
	['Single211_CMSSW_6_2_0_SLHC23_patch2_RECO-PU0_SimHits/edep_rec-TrivialComb/calib__calib_uncalib.root',  'simple', 'no #pi/e corr.',        False],
	['Single211_CMSSW_6_2_0_SLHC23_patch2_RECO-PU0_SimHits/edep_rec-IndComb/calib__calib_uncalib.root',      'simple', '#pi/e corr.',      False],
	['Single211_CMSSW_6_2_0_SLHC23_patch2_RECO-PU0_SimHits/edep_rec-IndComb/calib__calib_uncalib.root',      'gc',     '#pi/e corr.+GC', True]
	]

ROOT.gROOT.SetBatch(True)


ymin,ymax=0.05,0.5

gr={}
gr_ratio={}
titles={}
for ieta in xrange(0,6):
	gr[ieta]=[]
	gr_ratio[ieta]=[]
	pcount=0
	for f,algo,title, doFit in toPlot:
		fIn=ROOT.TFile.Open(f)
		obj=fIn.Get('res_%d_%s'%(ieta,algo))
		titles[ieta]=obj.GetTitle()
		
		obj.SetFillStyle(0)
		obj.SetMarkerStyle(20+4*(pcount%2))
		obj.SetMarkerColor(1+pcount/2)
		obj.SetLineColor(1+pcount/2)
		fullTitle=title
		if doFit:
			resolModel=None
			if resolModelType==0 :
				resolModel=ROOT.TF1('resolmodel_%d_%s'%(ieta,algo),"sqrt([0]*[0]/x+[1]*[1])",0,1000)
				resolModel.SetParameter(0,0.2);
				resolModel.SetParLimits(0,0,2);
				resolModel.SetParameter(1,0);
				resolModel.SetParLimits(1,0,1.0);
			else:
				resolModel=ROOT.TF1('resolmodel_%d_%s'%(ieta,algo),"sqrt([0]*[0]/x+[1]*[1]+[2]*[2]/(x*x))",0,1000)
				resolModel.SetParameter(2,0)
				resolModel.SetParLimits(2,0.01,0.1)
				resolModel.SetParameter(0,0.2);
				resolModel.SetParLimits(0,0,2);
				resolModel.SetParameter(1,0);
				resolModel.SetParLimits(1,0,1.0);
			resolModel.SetLineColor(1+pcount/2)
			resolModel.SetLineStyle(1+(pcount%2))
			obj.Fit(resolModel,'MRQ+')
			if resolModelType==0 :
				fullTitle = '%s :  #frac{%3.3f}{#sqrt{E}}#oplus%3.3f'%(title,resolModel.GetParameter(0),resolModel.GetParameter(1))
			else :
				fullTitle = '%s :  #frac{%3.3f}{#sqrt{E}}#oplus%3.3f#oplus#frac{%3.3f}{E}'%(title,resolModel.GetParameter(0),resolModel.GetParameter(1),resolModel.GetParameter(2))
		
		obj.SetTitle(fullTitle)
		gr[ieta].append(obj)
		ratioToFirstGr=getRatio(gr[ieta][0],obj)
		ratioToFirstGr.SetTitle(title)
		gr_ratio[ieta].append(ratioToFirstGr)
		fIn.Close()
		pcount=pcount+1

#show all
grecalhb=getTBresults()
if overlayTestBeam :  showResults(gr,grecalhb,ymin,ymax,'#sigma_{E}/E','resol')
else : showResults(gr,None,ymin,ymax,'#sigma_{E}/E','resol')
showResults(gr_ratio,None,0.5,1.5,'(#sigma_{E}/E)_{algo} / (#sigma_{E}/E)_{perfect cluster}','relresol')

