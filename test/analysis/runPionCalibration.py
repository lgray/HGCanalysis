#!/usr/bin/env python

import ROOT
import numpy as numpy
import array
import io, os, sys
import optparse
import commands
from UserCode.HGCanalysis.PlotUtils import *
from UserCode.HGCanalysis.HGCTree2Workspace import *

"""
draws the correlation for different sub-dets
"""
def computeSubdetectorWeights(enRanges,etaRanges,ws,xaxis,yaxis,outDir) :

    divx=len(etaRanges)/2
    while True:
        if divx*2 >= len(etaRanges) : break
        divx+=1

    #decode the combination which goes to the yaxis
    yAxisVars=[]
    yAxisWeights=[]
    for var in yaxis.split('+') :
        varComps=var.split('/')
        yAxisVars.append( varComps[0] )
        if len(varComps)!=2 : yAxisWeights.append( 1.0 )
        else: yAxisWeights.append( float(varComps[1]) )
    yaxis=''.join(yAxisVars)

    canvas=ROOT.TCanvas('c','c',350*divx,700)
    all2DScatters=[]
    all2DProfiles=[]
    etaSliceCombGr={}
    for ien in xrange(0,len(enRanges)):
        genEn_min=enRanges[ien][0]
        genEn_max=enRanges[ien][1]
        genEn_mean=0.5*(genEn_max+genEn_min)

        canvas.Clear()
        canvas.Divide(divx,2)   
        for ieta in xrange(0,len(etaRanges)):
            genEta_min=etaRanges[ieta][0]
            genEta_max=etaRanges[ieta][1]

            if not ieta in etaSliceCombGr:
                etaSliceCombGr[ieta]=ROOT.TGraphErrors()
                etaSliceCombGr[ieta].SetName('comb_%s'%ieta)
                etaSliceCombGr[ieta].SetTitle('%3.1f<|#eta|<%3.1f'%(genEta_min,genEta_max))
                etaSliceCombGr[ieta].SetLineColor(49-ieta)
                etaSliceCombGr[ieta].SetMarkerColor(49-ieta)
                etaSliceCombGr[ieta].SetMarkerStyle(ieta+20)

            #fill a new profile
            redData=ws.data('data_uncalib').reduce('en>=%f && en<=%f && eta>=%f && eta<=%f'%(genEn_min,genEn_max,genEta_min,genEta_max))
            if redData.numEntries()<10 : continue

            iprof=len(all2DScatters)
            all2DScatters.append( ROOT.TH2F('hprof%d%d'%(ieta,ien),';%s energy/E_{beam};%s energy/E_{beam};Events'%(xaxis,yaxis),75,0,1.5,75,0,1.5) )
            all2DScatters[iprof].SetDirectory(0)
            all2DScatters[iprof].Sumw2()
            for ientry in xrange(0,redData.numEntries()):
                entryVars=redData.get(ientry)
                xval=entryVars.find('lambda_emEn_%s'%xaxis).getVal()/genEn_mean
                yval=0
                for iyvar in xrange(0,len(yAxisVars)):
                    yval+=entryVars.find('lambda_emEn_%s'%yAxisVars[iyvar]).getVal()/yAxisWeights[iyvar]
                yval/=genEn_mean
                all2DScatters[iprof].Fill(xval,yval)
    
            #show profile and fit
            p=canvas.cd(ieta+1)
            p.SetLogz()
            p.SetRightMargin(0.1)
            
            #find fit range from projection quantiles to avoid biases
            projY=all2DScatters[iprof].ProjectionX('tmpproj',2)
            probSum = array.array('d', [0.05,0.4])
            if len(yAxisVars)==2 : probSum = array.array('d', [1.0,0.7])
            quantiles = array.array('d', [0.0]*len(probSum))
            projY.GetQuantiles(len(probSum), quantiles, probSum)
            projY.Delete()
            all2DProfiles.append( all2DScatters[iprof].ProfileX('%s_prof'%all2DScatters[iprof].GetName()) )
            all2DProfiles[iprof].SetMarkerStyle(20)
            all2DScatters[iprof].Draw('colz')
            all2DScatters[iprof].GetZaxis().SetTitleOffset(-0.5)            
            all2DProfiles[iprof].Draw('e1same')
            all2DProfiles[iprof].Fit('pol1','RQM+','same',quantiles[0],quantiles[1])
            combSlope,combSlope_err=all2DProfiles[iprof].GetFunction('pol1').GetParameter(0),all2DProfiles[iprof].GetFunction('pol1').GetParError(0)
            MyPaveText('#it{%3.1f<#eta<%3.1f}\\%s = (%3.2f#pm%3.2f) %s'%(genEta_min,genEta_max,yaxis,combSlope,combSlope_err,xaxis),0.4,0.95,0.8,0.8).SetTextSize(0.04)
            if ieta==0 : MyPaveText('#bf{CMS} #it{simulation}    Energy=%d GeV'%genEn_mean)

            if combSlope<=0: continue
            np=etaSliceCombGr[ieta].GetN()
            etaSliceCombGr[ieta].SetPoint(np,genEn_mean,combSlope)
            etaSliceCombGr[ieta].SetPointError(np,0,combSlope_err)
    
        #save
        canvas.Modified()
        canvas.Update()
        canvas.SaveAs('%s/%s_%s_%d_profile.png'%(outDir,xaxis,yaxis,ien))


    #show combination coefficient evolution
    ccombCanvas=ROOT.TCanvas('ccomb','ccomb',500,500)
    combGr=ROOT.TMultiGraph()
    combGr.SetName('%s_%s_comb'%(xaxis,yaxis))
    for ieta in xrange(0,len(etaRanges)): combGr.Add(etaSliceCombGr[ieta],'p')      
    combGr.Draw('a')
    combGr.GetYaxis().SetTitle('%s/%s response'%(yaxis,xaxis))
    combGr.GetXaxis().SetTitle('Beam energy [GeV]')
    combGr.GetYaxis().SetRangeUser(0,2)
    leg=ROOT.TLegend(0.6,0.6,0.9,0.94)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.03)
    for gr in combGr.GetListOfGraphs(): leg.AddEntry(gr,gr.GetTitle(),"p")
    leg.Draw()
    combGr.Fit('pol0')
    finalCombSlope, finalCombSlope_err=combGr.GetFunction('pol0').GetParameter(0),combGr.GetFunction('pol0').GetParError(0)
    MyPaveText('#bf{CMS} #it{simulation}').SetTextSize(0.04)
    ccombCanvas.Modified()
    ccombCanvas.Update()
    ccombCanvas.SaveAs('%s/%s_%s_combFactor.png'%(outDir,xaxis,yaxis))
   
    #show the combined energy sum
    allSums=[]
    allWeightedSums=[]
    for ien in xrange(0,len(enRanges)):
        genEn_min=enRanges[ien][0]
        genEn_max=enRanges[ien][1]
        genEn_mean=0.5*(genEn_max+genEn_min)

        canvas.Clear()
        canvas.Divide(divx,2)   
        
        for ieta in xrange(0,len(etaRanges)):
            genEta_min=etaRanges[ieta][0]
            genEta_max=etaRanges[ieta][1]

            if not ieta in etaSliceCombGr:
                etaSliceCombGr[ieta]=ROOT.TGraphErrors()
                etaSliceCombGr[ieta].SetName('comb_%s'%ieta)
                etaSliceCombGr[ieta].SetTitle('%3.1f<|#eta|<%3.1f'%(genEta_min,genEta_max))
                etaSliceCombGr[ieta].SetLineColor(49-ieta)
                etaSliceCombGr[ieta].SetMarkerColor(49-ieta)
                etaSliceCombGr[ieta].SetMarkerStyle(ieta+20)

            #fill a new profile
            redData=ws.data('data_uncalib').reduce('en>=%f && en<=%f && eta>=%f && eta<=%f'%(genEn_min,genEn_max,genEta_min,genEta_max))
            if redData.numEntries()<10 : continue            
            iprof=len(allSums)
            allSums.append( ROOT.TH1F('hsum%d%d'%(ieta,ien),';%s+%s energy/E_{beam};Events'%(xaxis,yaxis),25,0,2) )   
            allSums[iprof].SetDirectory(0)
            allSums[iprof].Sumw2()
            allWeightedSums.append( allSums[iprof].Clone('hwgtsum%d%d'%(ieta,ien)) )
            allWeightedSums[iprof].SetLineColor(ROOT.kRed)
            allWeightedSums[iprof].SetLineWidth(2)
            for ientry in xrange(0,redData.numEntries()):
                entryVars=redData.get(ientry)
                xval=entryVars.find('lambda_emEn_%s'%xaxis).getVal()/genEn_mean
                yval=0
                for iyvar in xrange(0,len(yAxisVars)):
                    yval+=entryVars.find('lambda_emEn_%s'%yAxisVars[iyvar]).getVal()/yAxisWeights[iyvar]
                yval/=genEn_mean
                allSums[iprof].Fill(xval+yval)
                allWeightedSums[iprof].Fill(xval+yval/finalCombSlope)

            #show raw sum and after weighted combination
            pad=canvas.cd(ieta+1)
            allSums[iprof].Draw('hist')
            allSums[iprof].GetYaxis().SetRangeUser(0,allSums[iprof].GetMaximum()*2)
            allSums[iprof].SetTitle("trivial combination")
            allWeightedSums[iprof].Draw('histsame')
            allWeightedSums[iprof].SetTitle("weighted combination")
            MyPaveText('#it{Energy=%d GeV, %3.1f<#eta<%3.1f}'%(genEn_mean,genEta_min,genEta_max),0.15,0.95,0.5,0.9).SetTextSize(0.04)
            if ieta==0 :
                MyPaveText('#bf{CMS} #it{simulation}')
                leg=ROOT.TLegend(0.6,0.6,0.9,0.94)
                leg.SetFillStyle(0)
                leg.SetBorderSize(0)
                leg.SetTextFont(42)
                leg.SetTextSize(0.03)
                leg.AddEntry(allSums[iprof],allSums[iprof].GetTitle(),'l')
                leg.AddEntry(allWeightedSums[iprof],allWeightedSums[iprof].GetTitle(),'l')
                leg.Draw()

        canvas.Modified()
        canvas.Update()
        canvas.SaveAs('%s/%s_%s_%d_comb.png'%(outDir,xaxis,yaxis,ien))

    #free memory
    for p in all2DScatters: p.Delete()
    for p in all2DProfiles: p.Delete()
    for p in allSums: p.Delete()
    for p in allWeightedSums : p.Delete()

    #save combination coefficient
    fOut=ROOT.TFile.Open('%s/%s%s_comb.root'%(outDir,xaxis,yaxis),'RECREATE')
    combGr.Write()
    combGr.GetFunction('pol0').Clone().Write('%s%s_combfunc'%(xaxis,yaxis))
    print '%s vs %s combination written @ %s'%(xaxis,yaxis,fOut.GetName())
    fOut.Close()

    return finalCombSlope, finalCombSlope_err





"""
draws the correlation between showerVolume and energy estimator
"""
def computeCompensationWeights(enRanges,etaRanges,ws,outDir):

    divx=len(etaRanges)/2
    while True:
        if divx*2 >= len(etaRanges) : break
        divx+=1

    #optimize in energy density ranges
    #uRanges=[[0,0.1],[0.1,0.2],[0.2,0.3],[0.3,0.5]]
    uRanges=[[0,0.2],[0.2,0.4],[0.4,0.6],[0.6,0.8]]

    weightOptimGr=[]
    aGr=ROOT.TGraphErrors()
    aGr.SetMarkerStyle(20)
    aGr.SetName('swweights_slope')
    bGr=aGr.Clone()
    bGr.SetName('swweights_offset')

    canvas=ROOT.TCanvas('c','c',350*divx,700)
    all2DScatters=[]
    all2DProfiles=[]
    all1DProfiles=[]
    for ien in xrange(0,len(enRanges)):
        genEn_min=enRanges[ien][0]
        genEn_max=enRanges[ien][1]
        genEn_mean=0.5*(genEn_max+genEn_min)

        #sums for weight optimization
        nSum=[0]*len(uRanges)
        densSum=[0]*len(uRanges)
        eSum=[0]*len(uRanges)
        e2Sum=[0]*len(uRanges)

        iprof1d=len(all1DProfiles)
        all1DProfiles.append( ROOT.TH1F('hprof1d%d'%ien,'Energy=%d GeV;U (shower density) [GeV/Volume];PDF'%genEn_mean,50,0,1) )
        all1DProfiles[iprof1d].SetDirectory(0)
        all1DProfiles[iprof1d].Sumw2()

        canvas.Clear()
        canvas.Divide(divx,2)   
        
        for ieta in xrange(0,len(etaRanges)):
            genEta_min=etaRanges[ieta][0]
            genEta_max=etaRanges[ieta][1]

            #fill a new profile
            redData=ws.data('data_uncalib_final').reduce('en>=%f && en<=%f && eta>=%f && eta<=%f'%(genEn_min,genEn_max,genEta_min,genEta_max))
            if redData.numEntries()<10 : continue

            iprof=len(all2DScatters)
            all2DScatters.append( ROOT.TH2F('hprof%d%d'%(ieta,ien),';E_{rec} [GeV];Shower density [GeV/Volume];Events',30,0.5*genEn_mean,3*genEn_mean,50,0,1) )
            all2DScatters[iprof].SetDirectory(0)
            all2DScatters[iprof].Sumw2()
            for ientry in xrange(0,redData.numEntries()):
                entryVars=redData.get(ientry)
                egen=entryVars.find('en').getVal()                
                erec=entryVars.find('lambda_emEn').getVal()                
                vol=entryVars.find('volume').getVal()
                dens=0
                if vol>0 : dens=erec/vol 
                all2DScatters[iprof].Fill(erec,dens)
                all1DProfiles[iprof1d].Fill(dens,1./redData.numEntries())

                #compute sums for maximization
                for iuRange in xrange(0,len(uRanges)):
                    if dens<uRanges[iuRange][0] or dens>uRanges[iuRange][1] : continue
                    nSum[iuRange]    += 1
                    densSum[iuRange] += dens
                    eSum[iuRange]    += erec/egen
                    e2Sum[iuRange]   += ROOT.TMath.Power(erec/egen,2)

            #show profile
            p=canvas.cd(ieta+1)
            p.SetLogz()
            p.SetRightMargin(0.1)
            
            all2DProfiles.append( all2DScatters[iprof].ProfileX('%s_prof'%all2DScatters[iprof].GetName()) )
            all2DProfiles[iprof].SetMarkerStyle(20)
            all2DScatters[iprof].Draw('colz')
            all2DScatters[iprof].GetZaxis().SetTitleOffset(-0.5)            
            all2DProfiles[iprof].Draw('e1same')
            MyPaveText('#it{%3.1f<#eta<%3.1f}'%(genEta_min,genEta_max),0.4,0.95,0.8,0.8).SetTextSize(0.04)
            if ieta==0 : MyPaveText('#bf{CMS} #it{simulation}    Energy=%d GeV'%genEn_mean)
            
        #save
        canvas.Modified()
        canvas.Update()
        canvas.SaveAs('%s/showerdens_%d_profile.png'%(outDir,ien))

        #now optimize weights for this energy
        weightOptimGr.append( ROOT.TGraph() )
        weightOptimGr[ien].SetName('sweights%d'%ien)
        weightOptimGr[ien].SetTitle('Energy=%d GeV'%genEn_mean)
        weightOptimGr[ien].SetMarkerStyle(20)
        for iuRange in xrange(0,len(uRanges)):
            if nSum[iuRange]<2 : continue
            np=weightOptimGr[ien].GetN()
            weightOptimGr[ien].SetPoint(np,densSum[iuRange]/nSum[iuRange],eSum[iuRange]/e2Sum[iuRange])
        weightOptimGr[ien].Fit('pol1','MQR+')

        #save evolution of the parameters for this energy
        a,a_err=weightOptimGr[ien].GetFunction('pol1').GetParameter(1),weightOptimGr[ien].GetFunction('pol1').GetParError(1)
        aGr.SetPoint(ien,genEn_mean,a)
        aGr.SetPointError(ien,0,a_err)
        b,b_err=weightOptimGr[ien].GetFunction('pol1').GetParameter(0),weightOptimGr[ien].GetFunction('pol1').GetParError(0)
        bGr.SetPoint(ien,genEn_mean,b)
        bGr.SetPointError(ien,0,b_err)
        
    #compare shower densities for different events
    canvas.Clear()
    canvas.SetWindowSize(500,500)
    canvas.SetCanvasSize(500,500)
    canvas.SetLogy()
    leg=ROOT.TLegend(0.6,0.6,0.9,0.94)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.03)
    for iprof1d in xrange(0,len(all1DProfiles)):
        all1DProfiles[iprof1d].SetLineColor(1+iprof1d)
        if iprof1d==0 : 
            all1DProfiles[iprof1d].Draw('hist')
            all1DProfiles[iprof1d].GetYaxis().SetRangeUser(1e-3,1.0)
        else:
            all1DProfiles[iprof1d].Draw('histsame')
        leg.AddEntry(all1DProfiles[iprof1d],all1DProfiles[iprof1d].GetTitle(),'l')
    leg.Draw()
    MyPaveText('#bf{CMS} #it{simulation}')
    canvas.Modified()
    canvas.Update()
    canvas.SaveAs('%s/showerdens.png'%outDir)

    #compare evolution
    canvas.Clear()
    canvas.SetWindowSize(1000,500)
    canvas.SetCanvasSize(1000,500)
    canvas.Divide(2,1)
    canvas.cd(2)    
    aGr.Draw('ap')
    aGr.GetXaxis().SetTitle('Energy [GeV]')
    aGr.GetYaxis().SetTitle('Weight function slope')
    aGr.GetYaxis().SetTitleOffset(1.4)
    aFunc=ROOT.TF1('afunc','[0]+[1]*exp(x*[2])',0,1000)
    aFunc.SetParLimits(0,-1,2)
    aFunc.SetParLimits(1,-2,0)
    aFunc.SetParLimits(2,-2,0)
    aGr.Fit(aFunc,'MR+')
    MyPaveText('slope(E)=q_{1}+q_{2}e^{(q_{3}E)}\\q_{1}=%3.4f#pm%3.4f\\q_{2}=%3.4f#pm%3.4f\\q_{3}=%3.4f#pm%3.4f'
            %(aFunc.GetParameter(0),aFunc.GetParError(0),
              aFunc.GetParameter(1),aFunc.GetParError(1),
              aFunc.GetParameter(2),aFunc.GetParError(2)),
               0.4,0.3,0.9,0.6)

    canvas.cd(1)
    bGr.Draw('ap')
    bGr.GetXaxis().SetTitle('Energy [GeV]')
    bGr.GetYaxis().SetTitle('Weight function offset')
    bGr.GetYaxis().SetTitleOffset(1.4)
    bFunc=ROOT.TF1('bfunc','[0]*(1-exp(x*[1]))+[2]',0,1000)
    bFunc.SetParLimits(0,-1,1)
    bFunc.SetParLimits(1,-1,0)
    bFunc.SetParLimits(2,-1,1)
    bGr.Fit(bFunc,'MR+')
    MyPaveText('offset(E)=p_{1}[1-e^{(p_{2}E)}]+p_{3}\\p_{1}=%3.4f#pm%3.4f\\p_{2}=%3.4f#pm%3.4f\\p_{3}=%3.4f#pm%3.4f'
               %(bFunc.GetParameter(0),bFunc.GetParError(0),
                 bFunc.GetParameter(1),bFunc.GetParError(1),
                 bFunc.GetParameter(2),bFunc.GetParError(2)),
               0.4,0.3,0.9,0.6)
    MyPaveText('#bf{CMS} #it{simulation}')

    canvas.Modified()
    canvas.Update()
    canvas.SaveAs('%s/swcompweights.png'%outDir)
    canvas.SaveAs('%s/swcompweights.C'%outDir)

    #save weights to file
    fOut=ROOT.TFile('%s/swcompweights.root'%outDir,'RECREATE')
    fOut.cd()
    aGr.Write()
    bGr.Write()
    for w in weightOptimGr : w.Write()
    fOut.Close()
    






"""
Adapts the workspace for pion calibration
"""
def adaptWorkspaceForPionCalibration(opt,outDir):
    url               = opt.input
    wsUrl             = opt.wsUrl
    treeVarName       = opt.treeVarName
    vetoTrackInt      = opt.vetoTrackInt
    vetoHEBLeaks      = opt.vetoHEBLeaks
    emCalibUrls       = {}
    if opt.emCalibUrl :
        for iurl in opt.emCalibUrl.split(','):
            key_url=iurl.split(':')
            emCalibUrls[ key_url[0] ] = key_url[1]

    #init weights and integration ranges
    subDets            = ['EE','HEF','HEB']
    integRanges        = [[1,1],[2,11],[12,21],[22,30],[31,31],[32,42],[43,54]]
    subDetRanges       = ['EE', 'EE',  'EE',   'EE',   'HEF',  'HEF',  'HEB']

    weights            = {}
    #weights["lambda"]    = [0.01, 0.036, 0.043,  0.056,  0.338,  0.273,  0.476]
    weights["lambda_em"] = [0.01, 0.036, 0.043,  0.056,  0.338,  0.273,  0.476]
    weightTitles={}
    #weightTitles["lambda"]  = "#lambda-based weights"
    weightTitles["lambda_em"]  = "#lambda-based + e.m. scale weights"

    #prepare workspace (if needed) and output
    if wsUrl is None :
        wsUrl=prepareWorkspace(url=url,integRanges=integRanges,vetoTrackInt=vetoTrackInt,vetoHEBLeaks=vetoHEBLeaks,treeVarName=treeVarName,addRaw=False)
    
    #get the workspace from the file
    wsOutF=ROOT.TFile.Open(wsUrl)
    ws=wsOutF.Get('w')
    wsOutF.Close()

    #readout material overburden file
    matParamBeforeHGCMap={}
    matF=ROOT.TFile.Open("/afs/cern.ch/user/p/psilva/work/HGCal/Calib/CMSSW_6_2_0_SLHC20/src/UserCode/HGCanalysis/data/HGCMaterialOverburden.root")
    matParamBeforeHGCGr=matF.Get("lambdaOverburden")
    
    #readout e.m. calibration for sub-detectors
    emCalibMap={}
    emCalibMapRes={}
    for key in emCalibUrls:
        if key in subDets:
            calibF=ROOT.TFile.Open(emCalibUrls[key])
            if calibF:
                print 'Reading out e.m. calibration for %s from %s'%(key,emCalibUrls[key])
                emCalibMap[key]=calibF.Get('lambda_calib').Clone('%s_lambda_calib'%key)
                emCalibMapRes[key]=calibF.Get('lambda_calib_res').Clone('%s_lambda_calib_res'%key)
                calibF.Close()

    #global calibration (tbu when we read a calibration file)
    calibPostFix='uncalib'
    
    #prepare energy estimators
    for subDet in subDets:
 
        funcArgs='{'
        funcFormula={}
        for wType in weights: funcFormula[wType]='('
        varCtr=0
        for ireg in xrange(0,len(integRanges)) :

            #require only sub-detector being looped over
            if subDetRanges[ireg] != subDet : continue
            
            #update the list of arguments and the formula for the new variable
            if varCtr>0 : funcArgs += ','
            funcArgs += 'edep%d'%ireg
            for wType in weights:
                if varCtr>0 : funcFormula[wType] += '+'
                funcFormula[wType] += '%f*@%d'%(weights[wType][ireg],varCtr)
            varCtr+=1

        funcArgs += '}'
        for wType in weights:

            #scale up to the e.m. scale, if available parameterization
            if subDet in emCalibMap and wType.find('_em')>0 :
                subdet_emscale_calib_offset = emCalibMap[subDet].GetParameter(1)
                subdet_emscale_calib_slope  = emCalibMap[subDet].GetParameter(0)
                funcFormula[wType] += '-%f)/%f'%(subdet_emscale_calib_offset,subdet_emscale_calib_slope)
            else:
                funcFormula[wType] += ')'

            theFunc=ws.factory("RooFormulaVar::%sEn_%s_func('%s',%s)"%(wType,subDet,funcFormula[wType],funcArgs))
            ws.data('data').addColumn( theFunc ) 


    #finally create dataset for calibration (only sub-detector sums at this point)
    uncalibDataVars=ROOT.RooArgSet(ws.var('en'), ws.var('eta'), ws.var('phi'), ws.var('length'), ws.var('volume'))
    for wType in weights: 
        for subDet in subDets:
            ws.factory('%sEn_%s[0,0,999999999]'%(wType,subDet))
            ws.var('%sEn_%s'%(wType,subDet)).SetTitle( '%s, %s'%(subDet,weightTitles[wType] ) )
            uncalibDataVars.add( ws.var('%sEn_%s'%(wType,subDet) ) )
    getattr(ws,'import')( ROOT.RooDataSet('data_%s'%calibPostFix,'data_%s'%calibPostFix,uncalibDataVars) )

    #fill the dataset
    for ientry in xrange(0,ws.data('data').numEntries()):

        entryVars=ws.data('data').get(ientry)
        newEntry=ROOT.RooArgSet()

        for baseVar in ['en','eta','phi','length','volume']: 
            ws.var(baseVar).setVal( entryVars.find(baseVar).getVal() )
            newEntry.add( ws.var(baseVar) )

        #apply material correction before HGC
        if matParamBeforeHGCGr:
            genEta = ROOT.TMath.Abs(ws.var('eta').getVal())
            edep0_corr  = ws.var('edep0').getVal()/ROOT.TMath.TanH(genEta) 
            edep0_corr  *= matParamBeforeHGCGr.Eval(genEta)*int(vetoTrackInt) 
            ws.var('edep0').setVal(edep0_corr)
                
        for wType in weights :
                    
            #per sub-detector
            totalEnVal=0
            for subDet in subDets:
                ienVal=entryVars.find('%sEn_%s_func'%(wType,subDet)).getVal()

                #compute residual, if available -> this should change to shower mean eta...
                if wType in emCalibMapRes:
                    emcalib_residual=calibMapRes[wType].Eval(ROOT.TMath.Abs(ws.var('eta').getVal()))
                    ienVal*=(1-emcalib_residual)
                    
                #add corrected value
                ws.var('%sEn_%s'%(wType,subDet)).setVal(ienVal)
                newEntry.add(ws.var('%sEn_%s'%(wType,subDet)))
                    
                #increment to the total
                totalEnVal+=ienVal
            
        #all filled, add new row
        ws.data('data_%s'%calibPostFix).add( newEntry )

    #save to file
    newWsUrl=wsUrl.replace('.root','_%s_pion.root'%calibPostFix)
    ws.writeToFile(newWsUrl,True)
    print 'Workspace adapted for pion calibration, stored @ %s'%newWsUrl
    return newWsUrl

"""
Steer the calibration study
"""
def runCalibrationStudy(opt):

    #prepare the output, get the workspace
    wsUrl=opt.wsUrl
    outDir="./"
    if wsUrl is None:
        outDir=os.path.basename(opt.input).replace('.root','')
        os.system('mkdir -p '+outDir)
    else:
        outDir=os.path.dirname(wsUrl)
    if wsUrl is None or wsUrl.find('pion')<0 : wsUrl=adaptWorkspaceForPionCalibration(opt,outDir=outDir)
    wsOutF=ROOT.TFile.Open(wsUrl)
    ws=wsOutF.Get('w')
    wsOutF.Close()

    #init phase space regions of interest
    etaRanges = [[1.6,1.75],[1.75,2.0],[2.0,2.25],[2.25,2.5],[2.5,2.75]]
    enRanges  = [[4,6],[9,11],[19,21],[39,41],[49,51],[74,75],[99,101],[124,126],[174,176],[249,251],[399,401]]
    #enRanges  = [[9,11],[19,21],[39,41],[49,51],[74,75],[99,101],[149,151],[249,251]]
    #enRanges = [[9,11],[19,21],[39,41],[49,51],[74,75],[99,101],[249,251]]
    #small stats for test
    #etaRanges = [[1.75,2.5]]
    #enRanges  = [[9,11],[29,31],[49,51]]

    weightTitles={}
    #weightTitles["lambda"]  = "#lambda-based weights"
    weightTitles["lambda_em"]  = "rescaled #lambda-based weights"


    #determine how to combine HEF and HEB
    hefhebCombSlope, hefhebCombSlope_err = 1.0, 0.0
    try:
        hefhebFin=ROOT.TFile.Open(opt.hefhebCombUrl)
        hefhebCombSlope, hefhebCombSlope_err = hefhebFin.Get('HEFHEB_combfunc').GetParameter(0), hefhebFin.Get('HEFHEB_combfunc').GetParError(0)
        hefhebFin.Close()
    except:
        print 'Will compute HEF vs HEB combination slope'
        hefhebCombSlope, hefhebCombSlope_err = computeSubdetectorWeights(enRanges=enRanges,etaRanges=etaRanges,ws=ws,xaxis='HEF',yaxis='HEB',outDir=outDir)
    ws.data('data_uncalib').addColumn( ws.factory("RooFormulaVar::lambda_emEn_HEFHEB('@0+@1/%f',{lambda_emEn_HEF,lambda_emEn_HEB})"%(hefhebCombSlope)) )

    #determine how to combine EE and HE(F+B)
    ehCombSlope, ehCombSlope_err = 1.0, 0.0
    if opt.noEE:
        for wgt in weightTitles:
            ws.data('data_uncalib').addColumn( ws.factory("RooFormulaVar::%sEnFunc('0*@0+@1',{%sEn_EE,%sEn_HEFHEB})"%(wgt,wgt,wgt) ) )
    else:
        try:
            ehFin=ROOT.TFile.Open(opt.ehCombUrl)
            ehCombSlope, ehCombSlope_err = ehFin.Get('EEHEFHEB_combfunc').GetParameter(0), ehFin.Get('EEHEFHEB_combfunc').GetParError(0)
            ehFin.Close()
        except:
            print 'Will compute EE vs HE(F+B) combination slope using %f coefficient for HEB'%(1./hefhebCombSlope)
            ehCombSlope, ehCombSlope_err = computeSubdetectorWeights(enRanges=enRanges,etaRanges=etaRanges,ws=ws,xaxis='EE',yaxis='HEF+HEB/%3.4f'%hefhebCombSlope,outDir=outDir)
        for wgt in weightTitles:
            ws.data('data_uncalib').addColumn( ws.factory("RooFormulaVar::%sEnFunc('@0+@1/%f',{%sEn_EE,%sEn_HEFHEB})"%(wgt,hefhebCombSlope,wgt,wgt)) )

    #create the final dataset for calibration
    print 'Will use the following combination of sub-detectors %d x EE + %3.4f x (HEF + %3.4f x HEB)'%(1-int(opt.noEE),1./hefhebCombSlope,1./ehCombSlope)
    uncalibDataVars=ROOT.RooArgSet(ws.var('en'), ws.var('eta'), ws.var('phi'),ws.var('length'),ws.var('volume'))
    for wType in weightTitles: 
        ws.factory('%sEn[0,0,999999999]'%wType)
        ws.var('%sEn'%wType).SetTitle( '%s'%weightTitles[wType] )
        uncalibDataVars.add( ws.var('%sEn'%wType) )
    getattr(ws,'import')( ROOT.RooDataSet('data_uncalib_final','data_uncalib_final',uncalibDataVars) )

    #read calibrations from file
    calibPostFix='uncalib'
    calibMap={}
    calibMapRes={}
    try:
        calibPostFix=os.path.basename(opt.calibUrl)
        calibPostFix='_'+calibPostFix.replace('.root','')
        calibF=ROOT.TFile.Open(opt.calibUrl)
        print 'Reading out calibrations from %s'%opt.calibUrl
        for wType in weightTitles:
            calibMap[wType]=calibF.Get('%s_calib'%wType).Clone()
            calibMapRes[wType]=calibF.Get('%s_calib_res'%wType).Clone()
        calibF.Close()
    except:
        print 'No calibration will be applied'

    #read sw compensation weights and calibrations from file
    swCompParams={}
    swCompCalibMap={}
    swCompCalibMapRes={}
    try:
        swF=ROOT.TFile.Open(opt.compWeights)
        print 'Reading out compensation weights from %s'%(opt.compWeights)
        swCompParams['slope']=swF.Get('swweights_slope').GetFunction('afunc')
        swCompParams['offset']=swF.Get('swweights_offset').GetFunction('bfunc')
        swF.Close()
        calibPostFix+='_swcomp'
        try:            
            calibF=ROOT.TFile.Open(opt.compCalib)
            print 'Reading out sw compensation calibrations from %s'%opt.compCalib
            for wType in weightTitles:
                swCompCalibMap[wType]=calibF.Get('%s_calib'%wType).Clone('swcompcalib')
                swCompCalibMapRes[wType]=calibF.Get('%s_calib_res'%wType).Clone('swcompcalibres')
            calibF.Close()
            calibPostFix+='_calib'
        except:
            print 'No calibration will be applied to sw compensation'
    except:
        print 'No compensation weights have been found'

    #fill the dataset
    for ientry in xrange(0,ws.data('data_uncalib').numEntries()):

        entryVars=ws.data('data_uncalib').get(ientry)
        newEntry=ROOT.RooArgSet()

        for baseVar in ['en','eta','phi','length','volume']: 
            ws.var(baseVar).setVal( entryVars.find(baseVar).getVal() )
            newEntry.add( ws.var(baseVar) )
            
        shVol=ws.var('volume').getVal()
    
        for wType in weightTitles :

            ienVal=entryVars.find('%sEnFunc'%wType).getVal()

            #apply calibration
            calib_offset, calib_slope = 0.0, 1.0
            if wType in calibMap:
                calib_offset=calibMap[wType].GetParameter(1)
                calib_slope=calibMap[wType].GetParameter(0)
            
            #compute residual, if available -> this should change to shower mean eta...
            calib_residual = 0.0
            if wType in calibMapRes:
                calib_residual=calibMapRes[wType].Eval(ROOT.TMath.Abs(ws.var('eta').getVal()))
                
            #calibrated energy estimator
            ienVal=((ienVal-calib_offset)/calib_slope)*(1-calib_residual)

            #software compensation weight
            if wType.find('_em')>0 and len(swCompParams):
                swCompWeight=1
                if shVol>0:
                    shDensity=ROOT.TMath.Min(ienVal/shVol,0.3)
                    #shDensity=ienVal/shVol
                    swCompWeight=swCompParams['slope'].Eval(ienVal)*shDensity+swCompParams['offset'].Eval(ienVal)
                ienVal=ienVal*swCompWeight
                swcomp_calib_offset, swcomp_calib_slope, swcomp_calib_residual = 0.0, 1.0, 0.0
                if wType in swCompCalibMap:
                    swcomp_calib_offset   = swCompCalibMap[wType].GetParameter(1)
                    swcomp_calib_slope    = swCompCalibMap[wType].GetParameter(0)
                    swcomp_calib_residual = swCompCalibMapRes[wType].Eval(ROOT.TMath.Abs(ws.var('eta').getVal()))
                ienVal=((ienVal-swcomp_calib_offset)/swcomp_calib_slope)*(1-swcomp_calib_residual)
                
            #add corrected value
            ws.var('%sEn'%wType).setVal(ienVal)
            newEntry.add(ws.var('%sEn'%wType))

        #all filled, add new row
        ws.data('data_uncalib_final').add(newEntry)


    #calibrate the energy estimators (split up in different energies and pseudo-rapidity ranges)
    nSigmasToFit=3.0
    calibGr={}
    resGr={}
    for wType in weightTitles:
        calibGr[wType]=ROOT.TMultiGraph()
        calibGr[wType].SetName('calib_%s'%wType)
        calibGr[wType].SetTitle( weightTitles[wType] )
        resGr[wType]=calibGr[wType].Clone('res_%s'%wType)

    for iEtaRange in xrange(0,len(etaRanges)):
        
        genEta_min=etaRanges[iEtaRange][0]
        genEta_max=etaRanges[iEtaRange][1]
        genEta_mean=0.5*(genEta_max+genEta_min)
        
        #keep track of eta slices separately
        etaSliceCalibGr={}
        etaSliceResGr={}
        iwtype=0
        for wType in weightTitles:
            iwtype=iwtype+1
            etaSliceCalibGr[wType]=ROOT.TGraphErrors()
            etaSliceCalibGr[wType].SetName('calib_%s_%s'%(iEtaRange,wType))
            etaSliceCalibGr[wType].SetTitle('%3.1f<|#eta|<%3.1f'%(genEta_min,genEta_max))
            etaSliceCalibGr[wType].SetLineColor(iwtype)
            etaSliceCalibGr[wType].SetMarkerColor(iwtype)
            etaSliceCalibGr[wType].SetMarkerStyle(iEtaRange+20)
            etaSliceResGr[wType]=etaSliceCalibGr[wType].Clone('res_%s_%s'%(iEtaRange,wType))

        for iEnRange in xrange(0,len(enRanges)):
            genEn_min = enRanges[iEnRange][0]
            genEn_max = enRanges[iEnRange][1]
            redData=ws.data('data_uncalib_final').reduce('en>=%f && en<=%f && eta>=%f && eta<=%f'%(genEn_min,genEn_max,genEta_min,genEta_max))
            if redData.numEntries()<10 : continue
            genEn_mean, genEn_sigma = redData.mean(ws.var('en')),  redData.sigma(ws.var('en'))

            for wType in weightTitles :

                #prepare the fit to this slice
                vName = '%sEn'%wType
                v_mean, v_sigma    = redData.mean(ws.function(vName)), redData.sigma(ws.function(vName))
                #v_min, v_max       = ROOT.TMath.Max(v_mean*0.5,v_mean-5*v_sigma), v_mean+5*v_sigma
                #v_fitMin, v_fitMax = ROOT.TMath.Max(1,v_mean-nSigmasToFit*v_sigma), v_mean+nSigmasToFit*v_sigma
                v_min, v_max       = ROOT.TMath.Max(1,v_mean-5*v_sigma), v_mean+5*v_sigma
                v_fitMin, v_fitMax = ROOT.TMath.Max(1,v_mean-nSigmasToFit*v_sigma), v_mean+nSigmasToFit*v_sigma

                #define PDF
                fitName          = 'range%d%d_%s'%(iEtaRange,iEnRange,vName)            
                ws.var(vName).setRange(fitName,v_min,v_max)
                ws.var(vName).setRange('fit_%s'%fitName,v_fitMin, v_fitMax)
                ws.factory('RooCBShape::resol_%s(%s,mean_%s[%f,%f,%f],sigma_%s[%f,%f,%f],alpha_%s[-1.0,-20.0,-0.001],n_%s[2,1,10])'%
                           (fitName,vName,
                            fitName,v_mean,v_min,v_max,
                            fitName,v_sigma,v_sigma*0.001, v_sigma*2,
                            fitName,
                            fitName)
                           )
                
                #fit
                fres = ws.pdf('resol_%s'%fitName).fitTo( redData, ROOT.RooFit.Range('fit_%s'%fitName), ROOT.RooFit.Save(True) )
                meanFit, meanFit_error   = ws.var('mean_%s'%fitName).getVal(), ws.var('mean_%s'%fitName).getError()
                sigmaFit, sigmaFit_error = ws.var('sigma_%s'%fitName).getVal(), ws.var('sigma_%s'%fitName).getError()

                #save results
                np=etaSliceCalibGr[wType].GetN()
                etaSliceCalibGr[wType].SetPoint(np,genEn_mean,meanFit)
                etaSliceCalibGr[wType].SetPointError(np,0,meanFit_error)
                etaSliceResGr[wType].SetPoint(np,genEn_mean,sigmaFit/meanFit)
                etaSliceResGr[wType].SetPointError(np,0,ROOT.TMath.Sqrt( ROOT.TMath.Power(meanFit*sigmaFit_error,2)+
                                                                          ROOT.TMath.Power(meanFit_error*sigmaFit,2) ) /
                                                    ROOT.TMath.Power(meanFit,2) )

                
                #for debug purposes
                showCalibrationFitResults( theVar=ws.var(vName),
                                           theData=redData,
                                           thePDF=ws.pdf('resol_%s'%fitName),
                                           theLabel='#it{Energy=%d GeV, %3.1f<#eta<%3.1f}\\#mu=%3.2f#pm%3.2f\\#sigma=%3.2f#pm%3.2f'%(genEn_mean,genEta_min,genEta_max,meanFit,meanFit_error,sigmaFit,sigmaFit_error),
                                           fitName=fitName,
                                           outDir=outDir)

        for wType in weightTitles: 
            calibGr[wType].Add(etaSliceCalibGr[wType],'p')
            resGr[wType].Add(etaSliceResGr[wType],'p')

    #derive calibration
    calibModel=ROOT.TF1('calibmodel',"[0]*x+[1]",0,800)
    #calibModel=ROOT.TF1('calibmodel',"x>100 ? [0]*x+[1] : [2]*x*x+[3]*x+[4]",0,1000)
    calibModel.SetLineWidth(1)
    for wType in weightTitles :
        calibGr[wType].Fit(calibModel,'MER+')
        calibGr[wType].GetFunction(calibModel.GetName()).SetLineColor(calibGr[wType].GetListOfGraphs().At(0).GetLineColor())

    #show results
    resCorrectionGr=showCalibrationCurves(calibGr=calibGr,calibRanges=etaRanges,outDir=outDir,calibPostFix=calibPostFix)
    showResolutionCurves(resGr=resGr,outDir=outDir,calibPostFix=calibPostFix,model=0)

    #save all to file
    calibModelRes=ROOT.TF1('calibmodelres',"[0]*x*x+[1]*x+[2]",1.45,3.1)
    calibF=ROOT.TFile.Open('%s/calib_%s.root'%(outDir,calibPostFix),'RECREATE')
    for wType in weightTitles :
        calibGr[wType].Write()
        calibGr[wType].GetFunction(calibModel.GetName()).Write('%s_calib'%wType)
        resCorrectionGr[wType].Write()
        resCorrectionGr[wType].Fit(calibModelRes,'WMR+')
        resCorrectionGr[wType].GetFunction(calibModelRes.GetName()).Write('%s_calib_res'%wType)
    calibF.Close()

    #compute compensation weights
    if len(swCompParams)==0:
        computeCompensationWeights(enRanges,etaRanges,ws,outDir)

"""
steer 
"""
def main():
    
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i',      '--in' ,      dest='input',         help='Input file',                                                     default=None)
    parser.add_option('-w',      '--ws' ,      dest='wsUrl',         help='Workspace file',                                                 default=None)
    parser.add_option('--emCalib' ,            dest='emCalibUrl',    help='em calibration files (e.g. EE:calib_ee.root,HEF:calib_hef.root', default=None)
    parser.add_option('--calib' ,              dest='calibUrl',      help='pion calibration file',                                          default=None)
    parser.add_option('--compWeights' ,        dest='compWeights',   help='file with software compensation weights',                        default=None)
    parser.add_option('--compCalib' ,          dest='compCalib',     help='calibration after compensation weights applied',                 default=None)
    parser.add_option('--vetoTrackInt',        dest='vetoTrackInt',  help='flag if tracker interactions should be removed',                 default=False, action="store_true")
    parser.add_option('--vetoHEBLeaks',        dest='vetoHEBLeaks',  help='flag if HEB leaks are allowed',                                  default=False, action='store_true')
    parser.add_option('--hefhebComb',          dest='hefhebCombUrl', help='Location of the parameterization for HEB/HEF combination',       default=None)
    parser.add_option('--noEE',                dest='noEE',          help='Assign weight 0 to EE',                                          default=False, action='store_true')
    parser.add_option('--ehComb',              dest='ehCombUrl',     help='Location of the parameterization for EE+HE(F+B) combination',    default=None)
    parser.add_option('-v',      '--var' ,     dest='treeVarName',   help='Variable to use as energy estimator',                            default='edep_sim')
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

    runCalibrationStudy(opt)
    
if __name__ == "__main__":
    sys.exit(main())
