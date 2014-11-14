import ROOT
from UserCode.HGCanalysis.PlotUtils import *
from array import array

#prepare histograms
steps=['sim','rec','clus']
stepTitle=['SimHit','RecHit','PFCluster']
#hitEnArray=[ i*0.5 for i in xrange(0,20)] + [i*5+10 for i in xrange(0,18)] + [i*50+100 for i in xrange(0,19)]
#nEnPts=len(hitEnArray)-1

baseHistos={
    'edep':     ROOT.TH2F('edep',   ';HGC layer;Energy [MIP];hits / event ',               54,0,54,100,0,50),
    'emeanPhi': ROOT.TH2F('dphi',   ';HGC layer;#phi(E wgt)-#phi(gen) [rad];hits / event', 54,0,54,50,-3,3),
    'emeanEta': ROOT.TH2F('deta',   ';HGC layer;#eta(E wgt)-#eta(gen);hits / event',       54,0,54,50,-3,3),
    'nhits':    ROOT.TH2F('nhits',  ';HGC layer;Number of hits;hits / event',              54,0,54,100,0,100),
    'edepdR':   ROOT.TH2F('edepdR', ';HGC layer;<#Delta R>',                               54,0,54,100,0,100)
    #'sihih':    ROOT.TH2F('sihih',  ';HGC layer;#sigma(#eta,#eta);hits / event',           54,0,54,50,0,2),
    #'sipip':    ROOT.TH2F('sipip',  ';HGC layer;#sigma(#phi,#phi);hits / event',           54,0,54,50,0,20),
    #'sipih':    ROOT.TH2F('sipih',  ';HGC layer;#sigma(#eta,#phi);hits / event',           54,0,54,50,-2,2)
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
fIn=ROOT.TFile.Open('Single211_CMSSW_6_2_0_SLHC20_lgray_SimHits_0.root')
HGC=fIn.Get('analysis/HGC')
for i in xrange(0,HGC.GetEntriesFast()):
    HGC.GetEntry(i)
    if HGC.genEn>11 or HGC.genEn<9 : continue
    #if HGC.genEn>31 or HGC.genEn<29 : continue
    #if HGC.genEn>51 or HGC.genEn<49 : continue
    nEvts=nEvts+1
    for ilay in xrange(0,HGC.nlay):
        for var in histos:
            for step in histos[var]:
                if getattr(HGC,'edep_%s'%(step))[ilay]==0: continue
                val=getattr(HGC,'%s_%s'%(var,step))[ilay]
                if var=='emeanPhi' : val=ROOT.TVector2.Phi_mpi_pi(val-HGC.genPhi)
                if var=='emeanEta' : val=val-HGC.genEta
                weight=1
                if var.find('edep')>=0 :
                    ybinToFill=histos[var][step].GetYaxis().FindBin(val)
                    weight=1./histos[var][step].GetYaxis().GetBinWidth(ybinToFill)
                histos[var][step].Fill(ilay, val, weight)

#show these profiles
customROOTstyle()
ROOT.gROOT.SetBatch(False)


layersToProfile=[1,15,31,44] #10,20,31,35,44,50]
allLegs=[]
allProjs=[]

#show the histograms
canvas = ROOT.TCanvas('c','c',800,800)
for var in histos:
    canvas.Clear()
    canvas.Divide(2,2)

    ipad=0
    longProfiles=[]
    for step in histos[var]:
        ipad=ipad+1

        histos[var][step].Scale(1./nEvts)

        #main distribution
        pad=canvas.cd(ipad)
        pad.SetLogz()
        pad.SetRightMargin(0.12)
        histos[var][step].Draw('colz')
        histos[var][step].GetYaxis().SetTitleOffset(1.7)
        histos[var][step].GetZaxis().SetTitleOffset(-0.3)
        MyPaveText(stepTitle[ipad-1],0.8)
        if ipad==1:
            simpt=MyPaveText('#bf{CMS} #it{simulation}')
            simpt.SetTextSize(0.035)
            simpt.SetTextFont(42)
        
        #profile
        #pad=canvas.cd(4)
        #pad.SetLogy()
        #if var.find('edep')>=0 : pad.SetLogx()
        longProfiles.append( ROOT.TGraphErrors() )
        longProfiles[ipad-1].SetName(step)
        longProfiles[ipad-1].SetTitle(stepTitle[ipad-1])
        longProfiles[ipad-1].SetMarkerStyle(20+ipad)
        longProfiles[ipad-1].SetFillStyle(0)
        #allLegs.append( ROOT.TLegend(0.7,0.65,0.95,0.9) )
        #nLegs=len(allLegs)-1
        #allLegs[nLegs].SetFillStyle(0)
        #allLegs[nLegs].SetBorderSize(0)
        #allLegs[nLegs].SetTextFont(42)
        #allLegs[nLegs].SetTextSize(0.035)
        for xbin in xrange(1,histos[var][step].GetXaxis().GetNbins()+1):
            projH=histos[var][step].ProjectionY('pfy_%d_%d'%(xbin,ipad),xbin,xbin)
            fixExtremities(projH)
            np=longProfiles[ipad-1].GetN()
            longProfiles[ipad-1].SetPoint(np,xbin,projH.GetMean())
            longProfiles[ipad-1].SetPointError(np,0,projH.GetMeanError())

            #draw in profile bins
            #if xbin in layersToProfile:
            #    iproj=len(allProjs)
            #    allProjs.append(projH.Clone())
            #    allProjs[iproj].SetLineColor(iproj%len(layersToProfile)+1)
            #    allProjs[iproj].SetLineWidth(iproj%2+2)
            #    if xbin==1 : 
            #        allProjs[iproj].Draw("hist")
            #        allProjs[iproj].GetYaxis().SetLabelSize(0.04)
            #        allProjs[iproj].GetYaxis().SetTitleSize(0.04)
            #        allProjs[iproj].GetYaxis().SetTitle("#hits / event")
            #    else : 
            #        allProjs[iproj].Draw('histsame')
            #    allProjs[iproj].SetTitle('Layer #%d'%xbin)
            #    allLegs[nLegs].AddEntry(allProjs[iproj],allProjs[iproj].GetTitle(),'l')
        #allLegs[nLegs].Draw()

    #draw the profiles on the last canvas
    pad=canvas.cd(4)
    longProfiles[0].Draw('ap')
    longProfiles[0].GetXaxis().SetTitleSize(0.04)
    longProfiles[0].GetXaxis().SetLabelSize(0.04)
    longProfiles[0].GetYaxis().SetTitleSize(0.04)
    longProfiles[0].GetYaxis().SetLabelSize(0.04)
    longProfiles[0].GetYaxis().SetTitle('<%s>'%histos[var]['sim'].GetYaxis().GetTitle())
    longProfiles[0].GetYaxis().SetTitleOffset(1.2)
    longProfiles[0].GetXaxis().SetTitle('HGC layer')
    for il in xrange(1,len(longProfiles)) : longProfiles[il].Draw('p')
    leg=pad.BuildLegend()
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.04)

    canvas.cd()
    canvas.Modified()
    canvas.Update()
    raw_input()
