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
creates an histogram and fits a gaussian 
"""
def fitGaussianToPeak(data,hrange=None,nbins=None):
    meanX,meanXerr=0,0
    try:
        data.InheritsFrom('TH1')
        maxBin=data.GetMaximumBin()
        fitRangeMin=data.GetXaxis().GetBinCenter(maxBin-7)
        fitRangeMax=data.GetXaxis().GetBinCenter(maxBin+7)
        data.Fit('gaus','LMRQ0+','',fitRangeMin,fitRangeMax)
        meanX,meanXerr=data.GetFunction('gaus').GetParameter(1),data.GetFunction('gaus').GetParError(1)
    except:
        if len(data)<5: return 0,0
        h=ROOT.TH1F('datah','',nbins,hrange[0],hrange[1])
        for d in data: h.Fill(d)
        maxBin=h.GetMaximumBin()
        fitRangeMin=h.GetXaxis().GetBinCenter(maxBin-7)
        fitRangeMax=h.GetXaxis().GetBinCenter(maxBin+7)
        h.Fit('gaus','LMRQ0+','',fitRangeMin,fitRangeMax)
        meanX,meanXerr=h.GetFunction('gaus').GetParameter(1),h.GetFunction('gaus').GetParError(1)
        h.Delete()

    return meanX,meanXerr

"""
wraps the computation of pi/e
"""
def computePiOverE(x,xresp,gr):

    #require minimum 5 events for pi/e...
    if len(x)>5:

        #response by mean or mode
        mean,     meanErr     = numpy.mean(x),     numpy.std(x)/numpy.sqrt(len(x))
        meanResp, meanRespErr = numpy.mean(xresp), numpy.std(xresp)/numpy.sqrt(len(xresp))
        mode,     modeErr     = fitGaussianToPeak(data=x,     hrange=(0,800), nbins=50)
        modeResp, modeRespErr = fitGaussianToPeak(data=xresp, hrange=(0,2),   nbins=50)
    
        np=gr.GetN()
        if mean<5 or mode<=0 or modeResp<0.1:
            gr.SetPoint(np,mean,meanResp)
            gr.SetPointError(np,meanErr,meanRespErr)
            return [mean,meanErr,meanResp,meanRespErr]
        else:
            gr.SetPoint(np,mode,modeResp)
            gr.SetPointError(np,modeErr,modeRespErr)
            return [mode,modeErr,modeResp,modeRespErr]

    #nothing done
    return []

"""
draws the correlation for different sub-dets
"""
def computeSubdetectorResponse(enRanges,etaRanges,xaxis,yaxis,ws,outDir):

    #auxiliary, for plotting purposes
    canvas   = ROOT.TCanvas('c','c',1000,500)
    allLegs  = []
    allProjs = []

    #pi/e model to adjust
    #responseFunc=ROOT.TF1('responseFunc','TMath::Min(TMath::Max([0]*(x-[1])/(x+[2])+[3],0.001),10.0)',0,1000)
    #responseFunc=ROOT.TF1('responseFunc','TMath::Min(TMath::Max([0]*(x-[1])/(x+[2])*exp(-[3]*x)+[4],0.001),10.0)',0,1000)
    responseFunc=ROOT.TF1('responseFunc','[0]*(1-exp(-[1]*x))/(1+exp(-[2]*x))*exp(-[3]*x)+[4]',0,1000)

    #determine the names for each axis
    xaxisName,hasCorrX='',False
    for var,corrFunc in xaxis: 
        xaxisName+=var
        if not (corrFunc is None) :  hasCorrX=True
    yaxisName,hasCorrY='',False
    for var,corrFunc in yaxis: 
        yaxisName+=var
        if not (corrFunc is None) : hasCorrY=True 
    postfix=''
    if hasCorrX : postfix += '_corrx'
    if hasCorrY : postfix += '_corry'

    #prepare the to derive the responses from the projections of the scatter plots
    all2DScatters={
        xaxisName : ROOT.TH2F('hscatrecx',';%s energy;E(rec)/E(gen);Events'%(xaxisName), 50,0,400,100,0,2.5),
        yaxisName : ROOT.TH2F('hscatrecy',';%s energy;E(rec)/E(gen);Events'%(yaxisName), 50,0,400,100,0,2.5)
        }
    for key in all2DScatters:
        all2DScatters[key].SetDirectory(0)
        all2DScatters[key].Sumw2()

    incResolutionScatter={}

    responseProfiles={xaxisName:ROOT.TGraphErrors(),
                      yaxisName:ROOT.TGraphErrors(),
                      xaxisName+'_gen':ROOT.TGraphErrors(),
                      yaxisName+'_gen':ROOT.TGraphErrors(),
                      xaxisName+yaxisName+'_comb':ROOT.TGraphErrors(),
                      xaxisName+yaxisName+'_combgen':ROOT.TGraphErrors(),
                      xaxisName+yaxisName+'_combdiag':ROOT.TGraphErrors()}
    for key in responseProfiles:
        marker=20
        if yaxisName in key: marker=24
        if yaxisName in key and xaxisName in key : marker=22
        responseProfiles[key].SetTitle(key)
        responseProfiles[key].SetMarkerStyle(marker)
        responseProfiles[key].SetFillStyle(0)
        responseProfiles[key].SetName('%s_prof'%key)

    #loop over available energies
    for ien in xrange(0,len(enRanges)):
        genEn_min=enRanges[ien][0]
        genEn_max=enRanges[ien][1]
        genEn_mean=0.5*(genEn_max+genEn_min)

        genEnKey='%3.0f'%genEn_mean
        all2DScatters[genEnKey]=ROOT.TH2F('hprof%d'%(ien),';%s energy/E_{beam};%s energy/E_{beam};Events'%(xaxisName,yaxisName),100,0,2.5,100,0,2.5)
        all2DScatters[genEnKey].SetDirectory(0)
        all2DScatters[genEnKey].Sumw2()
        incResolutionScatter[genEnKey]=ROOT.TH1F('hincprof%d'%(ien),';Total energy/E_{beam};Events',100,0,2.5)
        incResolutionScatter[genEnKey].SetDirectory(0)
        incResolutionScatter[genEnKey].Sumw2()
        incResolutionScatter[genEnKey+'_diag']=incResolutionScatter[genEnKey].Clone('hincprof%d_diag'%(ien))
        incResolutionScatter[genEnKey+'_diag'].SetDirectory(0)

        redData=ws.data('data').reduce('en>=%f && en<=%f && eta>=1.6 && eta<=2.8'%(genEn_min,genEn_max))
        if redData.numEntries()<10 : continue

        xvaluesAtY0,     yvaluesAtX0      = [], []
        xrespvaluesAtY0, yrespvaluesAtX0  = [], []
        xyvalues,        xyrespvalues     = [], []
        xydiagvalues,    xydiagrespvalues = [], []

        for ientry in xrange(0,redData.numEntries()):

            entryVars=redData.get(ientry)
            
            #raw estimate
            totalEnHEF=entryVars.find('en_HEF').getVal()
            totalEnHEB=entryVars.find('en_HEB').getVal()

            #corrected estimate
            xval=0
            for var,corrFunc in xaxis:
                ienVal=entryVars.find('en_%s'%var).getVal()
                if ienVal<0.1 : continue
                if not (corrFunc is None) :
                    ienVal /= corrFunc.Eval(ienVal)
                xval+=ienVal

            yval=0
            for var,corrFunc in yaxis:
                ienVal=entryVars.find('en_%s'%var).getVal()
                if ienVal<0.1 : continue
                if not (corrFunc is None):
                    ienVal /= corrFunc.Eval(ienVal)
                yval+=ienVal
            if xval+yval<0.1 : continue   


            xresp, yresp, combresp =xval/genEn_mean, yval/genEn_mean, (xval+yval)/genEn_mean

            all2DScatters[genEnKey].Fill(xresp,yresp)
            all2DScatters[yaxisName].Fill(yval,combresp)
            all2DScatters[xaxisName].Fill(xval,combresp)

            incResolutionScatter[genEnKey].Fill(combresp)
            xyvalues.append(xval+yval)
            xyrespvalues.append(combresp)

            if xresp/combresp>0.2 or yresp/combresp>0.2 :  
                incResolutionScatter[genEnKey+'_diag'].Fill(combresp)
                xydiagvalues.append(xval+yval)
                xydiagrespvalues.append(combresp)

            if xval<1 and yval>0.1:
                yvaluesAtX0.append(yval)
                yrespvaluesAtX0.append(yresp)

            if yval<1 and xval>0.1:
                xvaluesAtY0.append(xval)
                xrespvaluesAtY0.append(xresp)
        
        #pi/e 
        xPiOverE    = computePiOverE(x=xvaluesAtY0,  xresp=xrespvaluesAtY0,  gr=responseProfiles[xaxisName])
        if len(xPiOverE)==4:
            np=responseProfiles[xaxisName+'_gen'].GetN()
            responseProfiles[xaxisName+'_gen'].SetPoint(np,genEn_mean,xPiOverE[2])
            responseProfiles[xaxisName+'_gen'].SetPointError(np,0,xPiOverE[3])
        yPiOverE    = computePiOverE(x=yvaluesAtX0,  xresp=yrespvaluesAtX0,  gr=responseProfiles[yaxisName])
        if len(yPiOverE)==4:
            np=responseProfiles[yaxisName+'_gen'].GetN()
            responseProfiles[yaxisName+'_gen'].SetPoint(np,genEn_mean,yPiOverE[2])
            responseProfiles[yaxisName+'_gen'].SetPointError(np,0,yPiOverE[3])
        combPiOverE = computePiOverE(x=xyvalues,     xresp=xyrespvalues,     gr=responseProfiles[xaxisName+yaxisName+'_comb'])
        if len(combPiOverE)==4:
            np=responseProfiles[xaxisName+yaxisName+'_combgen'].GetN()
            responseProfiles[xaxisName+yaxisName+'_combgen'].SetPoint(np,genEn_mean,combPiOverE[2])
            responseProfiles[xaxisName+yaxisName+'_combgen'].SetPointError(np,0,combPiOverE[3])
        combdiagPiOverE = computePiOverE(x=xydiagvalues, xresp=xydiagrespvalues, gr=responseProfiles[xaxisName+yaxisName+'_combdiag'])
       
        #
        # SHOW ENERGY SLICE
        #
        canvas.Clear()
        canvas.Divide(2,1)            

        #scatter
        p=canvas.cd(1)
        p.SetLogz()
        p.SetRightMargin(0.1)
        all2DScatters[genEnKey].Draw('colz')
        all2DScatters[genEnKey].GetZaxis().SetTitleOffset(-0.5)            
        line=ROOT.TLine(0,1,1,0)
        line.SetLineStyle(7)
        line.Draw('same')

        #projections
        p=canvas.cd(2)
        nLegs=len(allLegs)
        allLegs.append( ROOT.TLegend(0.2,0.7,0.9,0.94) )
        allLegs[nLegs].SetFillStyle(0)
        allLegs[nLegs].SetBorderSize(0)
        allLegs[nLegs].SetTextFont(42)
        allLegs[nLegs].SetTextSize(0.035)        
        drawOpt='hist'
        for profKey in [genEnKey,genEnKey+'_diag']:
            totalFound=incResolutionScatter[profKey].Integral()
            if totalFound<0: continue
            profTitle, color, fill = 'inc', ROOT.kGray, 1001
            if 'diag' in profKey : profTitle, color, fill = 'inc, >0.2E_{rec}', 38, 3344
            incResolutionScatter[profKey].Rebin()
            incResolutionScatter[profKey].SetLineColor(color)
            incResolutionScatter[profKey].SetFillStyle(fill)
            incResolutionScatter[profKey].SetFillColor(color)
            incResolutionScatter[profKey].SetMarkerColor(color)
            incResolutionScatter[profKey].SetLineWidth(1)
            incResolutionScatter[profKey].Scale(1./totalFound)
            fixExtremities(incResolutionScatter[profKey])            
            incResolutionScatter[profKey].GetXaxis().SetTitle('E(rec)/E(gen)')
            incResolutionScatter[profKey].GetYaxis().SetRangeUser(0,0.3)
            title='<#pi/e>(%s) = '%profTitle
            if 'diag' in profKey :
                if len(combdiagPiOverE)==4 : title += '%3.2f'%combdiagPiOverE[2]
                else                       : title += 'n/a'
            else:
                if len(combPiOverE)==4     : title += '%3.2f'%combPiOverE[2]
                else                       : title += 'n/a'
            incResolutionScatter[profKey].SetTitle(title)
            incResolutionScatter[profKey].Draw(drawOpt)
            allLegs[nLegs].AddEntry( incResolutionScatter[profKey],title,'fp')
            drawOpt='histsame'

        for iaxis in xrange(0,2):
            nProjs=len(allProjs)
            title=''
            if iaxis==1:
                if len(xaxisName)==0 : continue
                allProjs.append( all2DScatters[genEnKey].ProjectionX('projx_%d'%nProjs,1,1) )
                allProjs[nProjs].SetTitle( xaxisName )
                title='<#pi/e>(%s) = '%xaxisName
                if len(xPiOverE)==4 : title += '%3.2f'%xPiOverE[2]
                else                : title += 'n/a'
            else:
                if len(yaxisName)==0 : continue
                allProjs.append( all2DScatters[genEnKey].ProjectionY('projy_%d'%nProjs,1,1) )
                allProjs[nProjs].SetTitle( yaxisName )
                title='<#pi/e>(%s) = '%yaxisName
                if len(yPiOverE)==4 : title += '%3.2f'%yPiOverE[2]
                else                : title += 'n/a'

            totalEvts=allProjs[nProjs].Integral()
            if totalEvts<1: continue
            allProjs[nProjs].Scale(1./totalEvts)
            fixExtremities(allProjs[nProjs])
            allProjs[nProjs].Rebin()
            allProjs[nProjs].SetLineColor(1+iaxis)
            allProjs[nProjs].SetMarkerColor(1+iaxis)
            allProjs[nProjs].SetMarkerStyle(20+4*iaxis)
            allProjs[nProjs].SetTitle(title)
            allProjs[nProjs].SetDirectory(0)
            allProjs[nProjs].Draw(drawOpt)
            allLegs[nLegs].AddEntry( allProjs[nProjs],title,'p')
            drawOpt='histsame'

        #all done
        allLegs[nLegs].Draw()    
        line1d=ROOT.TLine(1,0,1,0.3)
        line1d.SetLineStyle(7)
        line1d.Draw('same')
    
        p=canvas.cd(1)
        MyPaveText('#bf{CMS} #it{simulation}  Energy=%s GeV'%genEnKey)
        canvas.Modified()
        canvas.Update()
        canvas.SaveAs('%s/profile_%s%s.png'%(outDir,genEnKey,postfix))


    #final pi/e for x/y-axis separate
    canvas.Clear()
    canvas.SetWindowSize(500,500)
    canvas.SetCanvasSize(500,500)
    drawOpt='ap'
    nLegs=len(allLegs)
    allLegs.append( ROOT.TLegend(0.2,0.85,0.9,0.94) )
    allLegs[nLegs].SetFillStyle(0)
    allLegs[nLegs].SetBorderSize(0)
    allLegs[nLegs].SetTextFont(42)
    allLegs[nLegs].SetTextSize(0.035)   
    canvas.SetLogx()
    responseCtr=0
    for key in [xaxisName,yaxisName]:
        if responseProfiles[key].GetN()==0: continue
        responseProfiles[key].Draw(drawOpt)
        drawOpt='p'
        responseProfiles[key].GetYaxis().SetRangeUser(0.6,1.2)
        responseProfiles[key].GetYaxis().SetTitle('#pi/e mean or mode')
        responseProfiles[key].GetXaxis().SetTitle('E(rec) mean or mode [GeV]')
        allLegs[nLegs].AddEntry(responseProfiles[key],responseProfiles[key].GetTitle(),'p')
        #responseFunc.SetParameter(0,0.95)
        #responseFunc.SetParameter(1,0)
        #responseFunc.SetParLimits(1,0,100)
        #responseFunc.SetParameter(2,0)
        #responseFunc.SetParLimits(2,0,100)
        #responseFunc.SetParameter(3,0.001)
        #responseFunc.SetParLimits(3,-1,1)
        #responseFunc.SetParameter(4,1)
        if key=='HEB':
            responseFunc.SetParameter(0,1.0)
            responseFunc.SetParLimits(0,0,2)
            responseFunc.SetParameter(1,0.6)
            responseFunc.SetParLimits(1,0,10)
            responseFunc.SetParameter(2,0.12)
            responseFunc.SetParLimits(2,0,10)
            responseFunc.SetParameter(3,0.0002)
            responseFunc.SetParLimits(3,0,0.01)
            responseFunc.SetParameter(4,0)
            responseFunc.SetParLimits(4,0,2)
            responseProfiles[key].Fit(responseFunc,'MRQ+','',0,500)
            responseProfiles[key].GetFunction(responseFunc.GetName()).SetRange(0,1000)
        else:
            responseFunc.SetParameter(0,1.0)
            responseFunc.SetParLimits(0,0,2)
            responseFunc.SetParameter(1,0.6)
            responseFunc.SetParLimits(1,0,10)
            responseFunc.SetParameter(2,0.12)
            responseFunc.SetParLimits(2,0,10)
            responseFunc.SetParameter(3,0.0002)
            responseFunc.SetParLimits(3,0,0.01)
            responseFunc.SetParameter(4,0)
            responseFunc.SetParLimits(4,0,2)
            responseProfiles[key].Fit(responseFunc,'MR+','')
        responseProfiles[key].GetFunction(responseFunc.GetName()).SetLineStyle(responseCtr+1)
        responseCtr+=1

    allLegs[nLegs].SetNColumns(2)
    allLegs[nLegs].Draw()
    MyPaveText('#bf{CMS} #it{simulation}')
    canvas.Modified()
    canvas.Update()
    canvas.SaveAs('%s/piovereprofiles%s.png'%(outDir,postfix))
          
    #save response parameterisations
    fOut=ROOT.TFile.Open('%s/%s%s_response%s.root'%(outDir,xaxisName,yaxisName,postfix),'RECREATE')
    toReturn=[]
    try:
        toReturn.append( responseProfiles[xaxisName].GetFunction('responseFunc').Clone('%s_responseFunc'%xaxisName) )
        responseProfiles[xaxisName].Write()
        toReturn[0].Clone('%s_responseFunc'%xaxisName).Write()
    except:
        pass
    try:
        toReturn.append( responseProfiles[yaxisName].GetFunction('responseFunc').Clone('%s_responseFunc'%yaxisName) )
        responseProfiles[yaxisName].Write()
        toReturn[1].Clone('%s_responseFunc'%yaxisName).Write()
    except:
        pass
    print 'pi/e written for %s, %s in %s'%(xaxisName,yaxisName,fOut.GetName())
    fOut.Close()

    #all done here
    return toReturn


"""
draws the correlation between showerVolume and energy estimator
"""
#def computeCompensationWeights(enRanges,etaRanges,ws,outDir):
#
#    print '[computeCompensationWeights] will look at correlations between reconstructed energy and hit fraction (global) or shower density (local)'
#
#    #wgtfunc       = ROOT.TF1('wgtfunc','[0]+[1]*pow(x,[2])',0,2)
#    wgtfunc       = ROOT.TF1('wgtfunc','[0]+[1]*x',0,2)
#    weightOptimGr = {'rho':[],'c':[]}
#    aGr           = {'rho':[],'c':[]}
#    for key in aGr:
#        for ip in xrange(0,wgtfunc.GetNpar()):
#            aGr[key].append( ROOT.TGraphErrors() )
#            aGr[key][ip].SetMarkerStyle(20)
#            aGr[key][ip].SetName('%s_swweights_%d'%(key,ip))
#
#
#    canvas=ROOT.TCanvas('c','c',1200,800)
#    all2DScatters={'c_Si':[],'c_Sci':[],'rho_Si':[],'rho_Sci':[]}
#    all2DProfiles={'c_Si':[],'c_Sci':[],'rho_Si':[],'rho_Sci':[]}
#    all1DProfiles={'c_Si':[],'c_Sci':[],'rho_Si':[],'rho_Sci':[]}
#    for ien in xrange(0,len(enRanges)):
#
#        genEn_min=enRanges[ien][0]
#        genEn_max=enRanges[ien][1]
#        genEn_mean=0.5*(genEn_max+genEn_min)
#        iprof1d=len(all1DProfiles['rho_Si'])
#
#        #optimize in energy density ranges, based on the inclusive profile
#        all1DProfiles['rho_Si'].append( ROOT.TH1F('rho_Si_hprof1d_%d'%ien,'%d GeV;Shower density [GeV/Volume];PDF'%genEn_mean,25,0,2) )
#        all1DProfiles['rho_Si'][iprof1d].SetDirectory(0)
#        all1DProfiles['rho_Si'][iprof1d].Sumw2()
#        all1DProfiles['rho_Sci'].append( all1DProfiles['rho_Si'][iprof1d].Clone('rho_Sci_hprof1d_%d'%ien) )
#        all1DProfiles['rho_Sci'][iprof1d].SetDirectory(0)
#
#        all1DProfiles['c_Si'].append( ROOT.TH1F('c_Si_hprof1d_%d'%ien,'%d GeV;C(10 MIP);PDF'%genEn_mean,25,0.5,2.0) )
#        all1DProfiles['c_Si'][iprof1d].SetDirectory(0)
#        all1DProfiles['c_Si'][iprof1d].Sumw2()
#        all1DProfiles['c_Sci'].append( all1DProfiles['c_Si'][iprof1d].Clone( 'c_Sci_hprof1d_%d'%ien ) )
#        all1DProfiles['c_Sci'][iprof1d].SetDirectory(0)
#
#        #fill histograms and prepare for quantile computation
#        rho_values,c_values=[],[]
#        redIncData=ws.data('data').reduce('en>=%f && en<=%f'%(genEn_min,genEn_max))
#        nEntries=redIncData.numEntries()
#        if nEntries<10 : continue
#        for ientry in xrange(0,nEntries):
#            entryVars=redIncData.get(ientry)
#            for subDet in ['EE','HEF','HEB']:
#                rhoval=entryVars.find('rho_%s'%subDet).getVal()
#                cval=entryVars.find('c_%s'%subDet).getVal()
#                if subDet=='HEF':
#                    rho_values.append(rhoval)
#                    c_values.append(cval)
#                if rhoval>0 : all1DProfiles['rho_%s'%subDet][iprof1d].Fill(rhoval,1./nEntries)
#                if cval>0 : all1DProfiles['c_%s'%subDet][iprof1d].Fill(cval,1./nEntries)
#
#        #sums for weight optimization
#        rhoRanges,       cRanges       = [], []
#        prevRhoQuantile, prevCQuantile = 0,  0 
#        for q in [5,25,50,75]:
#            rhoQuantile, cQuantile = numpy.percentile(rho_values,q), numpy.percentile(c_values,q)
#            if prevRhoQuantile!=rhoQuantile : rhoRanges.append([prevRhoQuantile,rhoQuantile])
#            if prevCQuantile!=cQuantile     : cRanges.append([prevCQuantile,cQuantile])
#            prevRhoQuantile, prevCQuantile = rhoQuantile, cQuantile 
#        rhoRanges.append([prevRhoQuantile,ROOT.TMath.Min(numpy.percentile(rho_values,99),2)])
#        cRanges.append([prevCQuantile,ROOT.TMath.Min(numpy.percentile(c_values,99),2)])
#        print 'For E=%d GeV compensation weights will be optimised in the following en. density ranges'%genEn_mean
#        print 'Density : ',rhoRanges
#        print 'C       : ',cRanges
#
#        #arrays for optimisation:  delta^2 = [ w Erec / Egen - 1]^2 
#        # => [ Erec/Egen ] [ w Erec / Egen - 1 ] = 0
#        # ~  a(w.b +c ) = 0 <=> w = - ca / ab
#        # with ca = - [Erec/Egen] 
#        #      ab =  [Erec/Egen]^2
#        ca_rho_Sum, ab_rho_Sum = [0]*len(rhoRanges), [0]*len(rhoRanges)
#        ca_c_Sum,   ab_c_Sum   = [0]*len(cRanges), [0]*len(cRanges)
#        rho_values, c_values   = [], []
#        for irho in xrange(0,len(rhoRanges)) : rho_values.append( [] )
#        for ic in xrange(0,len(cRanges))     : c_values.append( [] )
#
#        iprof=len(all2DScatters['rho_EE'])
#        all2DScatters['rho_EE'].append( ROOT.TH2F('rho_EE_hprof2d_%d'%ien,';E_{rec} [GeV];Shower density [GeV/Volume];Events',30,0.5*genEn_mean,3*genEn_mean,25,0,2) )
#        all2DScatters['rho_EE'][iprof].SetDirectory(0)
#        all2DScatters['rho_EE'][iprof].Sumw2()
#        all2DScatters['rho_HEF'].append( all2DScatters['rho_EE'][iprof].Clone('rho_HEF_hprof2d_%d'%ien) )
#        all2DScatters['rho_HEF'][iprof].SetDirectory(0)
#        all2DScatters['rho_HEB'].append( all2DScatters['rho_EE'][iprof].Clone('rho_HEB_hprof2d_%d'%ien) )
#        all2DScatters['rho_HEB'][iprof].SetDirectory(0)
#        all2DScatters['c_EE'].append( ROOT.TH2F('c_EE_hprof2d_%d'%ien,';E_{rec} [GeV];C(10 MIP);Events',30,0.5*genEn_mean,3*genEn_mean,25,0.5,2) )
#        all2DScatters['c_EE'][iprof].SetDirectory(0)
#        all2DScatters['c_EE'][iprof].Sumw2()
#        all2DScatters['c_HEF'].append( all2DScatters['c_EE'][iprof].Clone('c_HEF_hprof2d_%d'%ien) )
#        all2DScatters['c_HEF'][iprof].SetDirectory(0)
#        all2DScatters['c_HEB'].append( all2DScatters['c_EE'][iprof].Clone('c_HEB_hprof2d_%d'%ien) )
#        all2DScatters['c_HEB'][iprof].SetDirectory(0)
#
#        for ieta in xrange(0,len(etaRanges)):
#            genEta_min=etaRanges[ieta][0]
#            genEta_max=etaRanges[ieta][1]
#            
#            #fill a new profile
#            redData=ws.data('data').reduce('en>=%f && en<=%f && eta>=%f && eta<=%f'%(genEn_min,genEn_max,genEta_min,genEta_max))
#            if redData.numEntries()<10 : continue
#
#            for ientry in xrange(0,redData.numEntries()):
#                entryVars=redData.get(ientry)
#
#                egen=entryVars.find('en').getVal()         
#                #FIXME use pi/e corrections
#                e_EE  = entryVars.find('en_EE').getVal() #*ws.var('k_EE').getVal()
#                e_HEF = entryVars.find('en_HEF').getVal() #*ws.var('k_HEF').getVal()
#                e_HEB = entryVars.find('en_HEB').getVal() #*ws.var('k_HEB').getVal()
#                e_tot = e_EE+e_HEF+e_HEB
#
#                for subDet in ['EE','HEF','HEB']:
#
#                    rhoval=entryVars.find('rho_%s'%subDet).getVal()
#                    all2DScatters['rho_%s'%subDet][iprof].Fill(e_tot,rhoval)
#
#                    cval=entryVars.find('c_%s'%subDet).getVal()
#                    all2DScatters['c_%s'%subDet][iprof].Fill(e_tot,cval)
#
#                    if subDet!='HEF': continue
#
#                    for irho in xrange(0,len(rhoRanges)):
#                        if rhoval<rhoRanges[irho][0] or rhoval>rhoRanges[irho][1] : continue
#                        rho_values[ irho ].append( rhoval )
#                        ca_rho_Sum[ irho ] += (e_tot/egen)
#                        ab_rho_Sum[ irho ] += (e_tot/egen)*(e_tot/egen)
#                        
#                    for ic in xrange(0,len(cRanges)):
#                        if cval<cRanges[ic][0] or cval>cRanges[ic][1] : continue
#                        c_values[ ic ].append( cval )
#                        ca_c_Sum[ ic ] += (e_tot/egen)
#                        ab_c_Sum[ ic ] += (e_tot/egen)*(e_tot/egen)
#
#        #now optimize weights for this energy
#        for key in weightOptimGr:
#            weightOptimGr[key].append( ROOT.TGraph() )
#            weightOptimGr[key][ien].SetName('%s_sweights_%d'%(key,ien))
#            weightOptimGr[key][ien].SetTitle('%d GeV'%genEn_mean)
#            weightOptimGr[key][ien].SetMarkerStyle(20)
#            if key=='rho':
#                for irho in xrange(0,len(rhoRanges)):
#                    if ab_rho_Sum[irho]==0: continue
#                    np=weightOptimGr[key][ien].GetN()
#                    weightOptimGr[key][ien].SetPoint(np,numpy.mean(rho_values[irho]),ca_rho_Sum[irho]/ab_rho_Sum[irho])
#            else:
#                for ic in xrange(0,len(cRanges)):
#                    if ab_c_Sum[ic]==0: continue
#                    np=weightOptimGr[key][ien].GetN()
#                    weightOptimGr[key][ien].SetPoint(np,numpy.mean(c_values[ic]),ca_c_Sum[ic]/ab_c_Sum[ic])
#                    
#            #save evolution of the parameters for this energy       
#            weightOptimGr[key][ien].Fit(wgtfunc,'MQR+')
#            for ip in xrange(0,wgtfunc.GetNpar()):
#                np=aGr[key][ip].GetN()
#                aGr[key][ip].SetPoint(np,genEn_mean,wgtfunc.GetParameter(ip))
#                aGr[key][ip].SetPointError(np,0,wgtfunc.GetParError(ip))
#        
#        #show profile
#        canvas.Clear()
#        canvas.Divide(3,2)   
#        ipad=0
#        for var in ['rho','c']:
#            for subDet in ['EE','HEF','HEB']: 
#                ipad+=1
#                p=canvas.cd(ipad)
#                p.SetLogz()
#                p.SetRightMargin(0.1)
#                key='%s_%s'%(var,subDet)
#                all2DProfiles[key].append( all2DScatters[key][iprof].ProfileY('%s_prof'%all2DScatters[key][iprof].GetName()) )
#                all2DProfiles[key][iprof].SetMarkerStyle(20)
#                all2DScatters[key][iprof].Draw('colz')
#                all2DScatters[key][iprof].GetZaxis().SetTitleOffset(-0.5)            
#                all2DScatters[key][iprof].GetYaxis().SetTitleOffset(1.2)            
#                #all2DProfiles[key][iprof].Draw('e1same')
#                MyPaveText('[ %s ]'%subDet,0.8,0.95,0.95,0.99)
#                if ipad>1: continue
#                MyPaveText('#bf{CMS} #it{simulation} Energy=%d'%genEn_mean)
#        canvas.Modified()
#        canvas.Update()
#        canvas.SaveAs('%s/showerdens_%d_profile.png'%(outDir,ien))
#       
#    #compare shower densities for different events
#    canvas.Clear()
#    canvas.Divide(3,2)   
#    ipad=0
#    for var in ['rho','c']:
#        for subDet in ['EE','HEF','HEB']:
#            key='%s_%s'%(var,subDet)
#            ipad+=1
#            p=canvas.cd(ipad)
#            p.Clear()
#            p.SetLogy(True)
#            for iprof1d in xrange(0,len(all1DProfiles[key])):
#                all1DProfiles[key][iprof1d].SetLineColor(45-2*iprof1d)
#                all1DProfiles[key][iprof1d].SetMarkerColor(45-2*iprof1d)
#                all1DProfiles[key][iprof1d].SetMarkerStyle(1)
#                all1DProfiles[key][iprof1d].SetLineWidth(2)
#                if iprof1d==0 : 
#                    all1DProfiles[key][iprof1d].Draw('hist')
#                    all1DProfiles[key][iprof1d].GetYaxis().SetRangeUser(1e-3,1.0)
#                    all1DProfiles[key][iprof1d].GetYaxis().SetTitleOffset(1.2)
#                else:
#                    all1DProfiles[key][iprof1d].Draw('histsame')
#
#            MyPaveText('[ %s ]'%subDet,0.8,0.95,0.95,0.99)
#
#            if ipad>1: continue
#            MyPaveText('#bf{CMS} #it{simulation}')
#            leg=p.BuildLegend(0.6,0.6,0.9,0.94)
#            leg.SetFillStyle(0)
#            leg.SetBorderSize(0)
#            leg.SetTextFont(42)
#            leg.SetTextSize(0.03)
#
#    canvas.Modified()
#    canvas.Update()
#    canvas.SaveAs('%s/showerdens.png'%outDir)
#
# 
#    #compare evolution
#    npar=wgtfunc.GetNpar()
#    canvas.Clear()
#    canvas.SetWindowSize(400*npar,400)
#    canvas.SetCanvasSize(400*npar,400)
#    canvas.Divide(npar,1)
#    for ip in xrange(0,wgtfunc.GetNpar()):
#        p=canvas.cd(ip+1)
#        p.SetLogy(False)
#        key='rho'
#        aGr[key][ip].Draw('ap')
#        aGr[key][ip].GetXaxis().SetTitle('Energy [GeV]')
#        aGr[key][ip].GetYaxis().SetTitle('Weight function parameter # %d'%ip)
#        aGr[key][ip].GetYaxis().SetTitleOffset(1.4)
#        aFunc=ROOT.TF1('%s_evfunc_%d'%(key,ip),'[0]+[1]*(1-exp(-x*[2]))',0,1000)
#        aFunc.SetParLimits(0,-100,100)
#        aFunc.SetParLimits(1,-50,50)
#        aFunc.SetParLimits(2,0.01,100)
#        aGr[key][ip].Fit(aFunc,'MR+')
#        if ip==0 : MyPaveText('#bf{CMS} #it{simulation}')
#        MyPaveText('Parameter evolution\\p_{%d}(E)=q_{1}+q_{2}(1-e^{-q_{3}E})\\q_{1}=%3.3f#pm%3.3f\\q_{2}=%3.3f#pm%3.3f\\q_{3}=%3.3f#pm%3.3f'
#                   %(ip,aFunc.GetParameter(0),aFunc.GetParError(0),
#                     aFunc.GetParameter(1),aFunc.GetParError(1),
#                     aFunc.GetParameter(2),aFunc.GetParError(2)),
#                   0.4,0.3,0.9,0.6).SetTextSize(0.035)
#    canvas.Modified()
#    canvas.Update()
#    canvas.SaveAs('%s/swcompweights.png'%(outDir))
#    canvas.SaveAs('%s/swcompweights.C'%(outDir))
#
#
#    #save weights and weight function to file
#    swcompUrl='%s/swcompweights.root'%outDir
#    fOut=ROOT.TFile(swcompUrl,'RECREATE')
#    fOut.cd()
#    for key in aGr : 
#        for gr in aGr[key] : gr.Write()
#        for w in weightOptimGr[key] : w.Write()
#    wgtfunc.Write()
#    fOut.Close()
#    
#    #all done
#    return swcompUrl
#

"""
Adapts the workspace for pion calibration
"""
def adaptWorkspaceForPionCalibration(opt,outDir):

    wsUrl             = opt.wsUrl

    #prepare workspace (if needed) and output
    if wsUrl is None :

        #electromagnetic energy scale of each section
        emCalibMap={'EE':None,'HEF':None,'HEB':None}
        if opt.emCalibUrl :
            for iurl in opt.emCalibUrl.split(','):
                subDet,url=iurl.split(':')
                print 'Replacing default energy scale in %s with calibration from %s'%(subDet,url)
                calibF=ROOT.TFile.Open(url)
                #emCalibMap[subDet]=(calibF.Get('simple_calib'),calibF.Get('calib_3_simple_res'))
                #don't apply residuals! dangerous for low energy!!! 
                emCalibMap[subDet]=(calibF.Get('simple_calib'),None)
                calibF.Close()

        #readout material overburden file
        matParamBeforeHGCMap={}
        matFurl='%s/src/UserCode/HGCanalysis/data/HGCMaterialOverburden.root'%os.environ['CMSSW_BASE']
        matF=ROOT.TFile.Open(matFurl)
        matParamBeforeHGCGr=matF.Get("lambdaOverburden")
        matF.Close()
        if not( matParamBeforeHGCGr is None) : print 'Material overburden has been read from %s'%matFurl

        #init weighting scheme
        weightingScheme={
            'EE': [([1, 1 ], matParamBeforeHGCGr, 0.010, emCalibMap['EE']),
                   ([2, 11], None,                0.036, emCalibMap['EE']),
                   ([12,21], None,                0.043, emCalibMap['EE']),
                   ([22,30], None,                0.056, emCalibMap['EE'])],
            'HEF':[([31,31], None,                0.338, emCalibMap['HEF']),
                   ([32,42], None,                0.273, emCalibMap['HEF'])],
            'HEB':[([43,54], None,                0.475, emCalibMap['HEB'])]
            }
        
        #prepare the workspace and get new url
        wsUrl=prepareWorkspace(url=opt.input,weightingScheme=weightingScheme,vetoTrackInt=opt.vetoTrackInt,vetoHEBLeaks=opt.vetoHEBLeaks,treeVarName=opt.treeVarName)
    
    #all done here
    print 'Workspace ready for pion calibration, stored @ %s'%wsUrl

    #return the workspace
    wsOutF=ROOT.TFile.Open(wsUrl)
    ws=wsOutF.Get('w')
    wsOutF.Close()
    return ws


"""
Steer the calibration study
"""
def runCalibrationStudy(opt):

    #prepare the output
    outDir="./"
    if opt.wsUrl is None:
        outDir=os.path.basename(opt.input).replace('.root','')
        os.system('mkdir -p '+outDir)
    else:
        outDir=os.path.dirname(opt.wsUrl)
    
    #get the workspace
    ws=adaptWorkspaceForPionCalibration(opt,outDir=outDir)
    
    #init phase space regions of interest
    etaRanges = [[1.55,1.75],[1.75,2.0],[2.0,2.25],[2.25,2.5],[2.5,2.7],[2.7,2.9]]
    enRanges  = [[1.8,2.2],[2.8,3.2],[4.5,5.5],[7.5,8.5],[9,11],[19,21],[39,41],[49,51],[74,76],[99,101],[124,126],[174,176],[249,251],[399,401],[499,501]]

    pioverE_EE, pioverE_HEF, pioverE_HEB = None, None,None

    #get pi/e for HEB
    try:
        hebFin=ROOT.TFile.Open(opt.hebRespUrl)
        pioverE_HEB  = hebFin.Get('HEB_responseFunc')
        print 'Readout pi/e response for HEB from %s'%hebFin.GetName()
        if opt.noEE and opt.noHEF:
            computeSubdetectorResponse(enRanges=enRanges,etaRanges=etaRanges,xaxis=[('HEB',pioverE_HEB)],yaxis=[],ws=ws,outDir=outDir)
        hebFin.Close()
    except:
        if opt.noEE and opt.noHEF:
            print 'Will compute pi/e response for HEB sub-detector'
            pioverE_HEB = computeSubdetectorResponse(enRanges=enRanges,etaRanges=etaRanges,xaxis=[('HEB',pioverE_HEB)],yaxis=[],ws=ws,outDir=outDir)[0]

    #get pi/e for HEF
    if opt.noHEF==False:
        try:
            hefhebFin=ROOT.TFile.Open(opt.hefRespUrl)       
            pioverE_HEF  = hefhebFin.Get('HEF_responseFunc')
            print 'Readout pi/e responses for HEF from %s'%hefhebFin.GetName()
            if opt.noEE :
                computeSubdetectorResponse(enRanges=enRanges,           etaRanges=etaRanges,
                                           xaxis=[('HEF',pioverE_HEF)], yaxis=[('HEB',pioverE_HEB)],
                                           ws=ws,                       outDir=outDir)
            hefhebFin.Close()
        except:
            if opt.noEE :
                print 'Will compute pi/e response for HEF  sub-detector'
                pioverE_HEF,_ = computeSubdetectorResponse(enRanges=enRanges,etaRanges=etaRanges,
                                                           xaxis=[('HEF',None)],yaxis=[('HEB',None)],
                                                           ws=ws,outDir=outDir)

    #get pi/e for EE
    if opt.noEE==False:
        try:
            ehFin=ROOT.TFile.Open(opt.eeRespUrl)
            pioverE_EE = ehFin.Get('EE_responseFunc')
            print 'Readout pi/e response for EE from %s'%ehFin.GetName()
            computeSubdetectorResponse(enRanges=enRanges,etaRanges=etaRanges,
                                       xaxis=[('EE',pioverE_EE)],yaxis=[('HEF',pioverE_HEF),('HEB',pioverE_HEB)],
                                       ws=ws,outDir=outDir)
            ehFin.Close()
        except:
            print 'Will compute pi/e for EE'
            pioverE_EE,_ =computeSubdetectorResponse(enRanges=enRanges,etaRanges=etaRanges,
                                                     xaxis=[('EE',None)],yaxis=[('HEF',pioverE_HEF),('HEB',pioverE_HEB)],
                                                     ws=ws,outDir=outDir)

    #read sw compensation weights
    #    swCompParams=[]
    #    if opt.compWeights is None:
    #        opt.compWeights=computeCompensationWeights(enRanges,etaRanges,ws,outDir)
    #    swF=ROOT.TFile.Open(opt.compWeights)
    #    print 'Reading out compensation weights from %s'%(opt.compWeights)
    #    swwgtfunc=swF.Get('wgtfunc')
    #    for ip in xrange(0,swwgtfunc.GetNpar()):
    #        swCompParams.append( swF.Get('rho_swweights_%d'%ip).GetFunction('rho_evfunc_%d'%ip) )
    #        #swCompParams.append( swF.Get('rho_swweights_%d'%ip) )
    #    swF.Close()


    #nothing else to be done
    if opt.noResCalib==True : 
        print 'Bailing out, no residual calibration needs to be run'
        return
    
    #read calibrations from file, if available
    weightTitles={'simple':'Simple sum ','gc':'Global compensation'} #,'lc':'Local compensation'}
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
            calibMapRes[wType]=calibF.Get('calib_3_%s_res'%wType).Clone()  #use the 2.0-2.25 eta range for residuals
        calibF.Close()
    except:
        print 'No calibration will be applied'

    #create the dataset for calibration
    uncalibDataVars=ROOT.RooArgSet(ws.var('en'), ws.var('eta'), ws.var('phi'))
    for wType in weightTitles: 
        ws.factory('%sEn[0,0,999999999]'%wType)
        ws.var('%sEn'%wType).SetTitle( '%s'%weightTitles[wType] )
        uncalibDataVars.add( ws.var('%sEn'%wType) )
    getattr(ws,'import')( ROOT.RooDataSet('data_uncalib_final','data_uncalib_final',uncalibDataVars) )
    for ientry in xrange(0,ws.data('data').numEntries()):

        entryVars=ws.data('data').get(ientry)
        newEntry=ROOT.RooArgSet()

        for baseVar in ['en','eta','phi']:
            ws.var(baseVar).setVal( entryVars.find(baseVar).getVal() )
            newEntry.add( ws.var(baseVar) )

        enEstimators={}

        #simple sum
        e_EE    = entryVars.find('en_EE').getVal()
        if opt.noEE : e_EE=0
        e_HEF   = entryVars.find('en_HEF').getVal()
        if opt.noHEF : e_HEF=0
        e_HEB   = entryVars.find('en_HEB').getVal()

        if not (pioverE_HEB is None) and not opt.noComp and e_HEB>0.1: 
            e_HEB /= pioverE_HEB.Eval(e_HEB) 
        if not (pioverE_HEF is None) and not opt.noComp and e_HEF>0.1: 
            e_HEF /= pioverE_HEF.Eval(e_HEF)
        if not (pioverE_EE is None) and not opt.noComp and e_EE>0.1: 
            e_EE /= pioverE_EE.Eval(e_EE)

        e_tot   = e_EE + e_HEF + e_HEB
        enEstimators['simple']=e_tot

        #global compensation weights
        c_EE    = entryVars.find('c_EE').getVal()
        c_HEF   = entryVars.find('c_HEF').getVal()
        c_HEB   = entryVars.find('c_HEB').getVal()
        e_c_tot = e_EE + c_HEF*e_HEF + e_HEB
        enEstimators['gc']=e_c_tot

        #software compensation weight
        #        e_rho_tot = e_tot
        #        rho_HEF = entryVars.find('rho_HEF').getVal()
        #        if e_HEF > 0 and rho_HEF>0:
        #            for ip in xrange(0,swwgtfunc.GetNpar()):
        #                swwgtfunc.SetParameter(ip,swCompParams[ip].Eval(e_tot))
        #            swWgt=swwgtfunc.Eval(ROOT.TMath.Min(rho_HEF,2.0))
        #            e_rho_tot = e_EE + swWgt*e_HEF + e_HEB
        #        enEstimators['lc']=e_rho_tot

        #now apply calibration
        for wType in weightTitles :

            ienVal=enEstimators[wType]

            if wType in calibMap:
                ienVal=calibMap[wType].GetX(ienVal)

            if wType in calibMapRes:
                resCorr=calibMapRes[wType].Eval(ienVal)
                ienVal*=(1-resCorr/100.)
                
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
                v_mean, v_sigma, v_skew    = redData.mean(ws.function(vName)), redData.sigma(ws.function(vName)), redData.skewness(ws.function(vName))
                #v_min, v_max       = ROOT.TMath.Max(v_mean*0.5,v_mean-5*v_sigma), v_mean+5*v_sigma
                #v_fitMin, v_fitMax = ROOT.TMath.Max(1,v_mean-nSigmasToFit*v_sigma), v_mean+nSigmasToFit*v_sigma
                v_min, v_max       = ROOT.TMath.Max(1,v_mean-3*v_sigma), v_mean+3*v_sigma
                v_fitMin, v_fitMax = ROOT.TMath.Max(1,v_mean-nSigmasToFit*v_sigma), v_mean+nSigmasToFit*v_sigma

                #define PDF
                fitName          = 'range%d%d_%s'%(iEtaRange,iEnRange,vName)            
                ws.var(vName).setRange(fitName,v_min,v_max)
                ws.var(vName).setRange('fit_%s'%fitName,v_fitMin, v_fitMax)
                iniAlphaValue=0
                if v_skew<0 : iniAlphaValue=2
                if v_skew>0 : iniAlphaValue=-2
                ws.factory('RooCBShape::resol_%s(%s,mean_%s[%f,%f,%f],sigma_%s[%f,%f,%f],alpha_%s[%f,-10.0,10.0],n_%s[2,1,3])'%
                           (fitName,vName,
                            fitName,v_mean,v_min,v_max,
                            fitName,v_sigma,v_sigma*0.001, v_sigma*1.5,
                            fitName,iniAlphaValue,
                            fitName)
                           ) 


                #fit
                theVar=ws.var(vName)
                thePDF=ws.pdf('resol_%s'%fitName)
                fres = thePDF.fitTo( redData, ROOT.RooFit.Range('fit_%s'%fitName), ROOT.RooFit.Save(True) )
                meanFit, meanFit_error   = ws.var('mean_%s'%fitName).getVal(), ws.var('mean_%s'%fitName).getError()
                sigmaFit, sigmaFit_error = ws.var('sigma_%s'%fitName).getVal(), ws.var('sigma_%s'%fitName).getError()
                scanStep=(v_max-v_min)*1e-4
                tol=scanStep/4
                effSigma = ROOT.getEffSigma(theVar,thePDF,v_min,v_max,scanStep,tol)
                sigmaEffVal=0.5*(effSigma.second-effSigma.first)
                
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
                                           theLabel='#it{Energy=%d GeV, %3.1f<#eta<%3.1f}\\#mu=%3.2f#pm%3.2f\\#sigma=%3.2f#pm%3.2f\\#sigma_{eff}=%3.2f'%(genEn_mean,genEta_min,genEta_max,meanFit,meanFit_error,sigmaFit,sigmaFit_error,sigmaEffVal),
                                           fitName=fitName,
                                           outDir=outDir)

        for wType in weightTitles: 
            calibGr[wType].Add(etaSliceCalibGr[wType],'p')
            resGr[wType].Add(etaSliceResGr[wType],'p')

    #derive calibration
    calibModel=ROOT.TF1('calibmodel',"[0]*x",0,1000)
    calibModel.SetLineWidth(1)
    for wType in weightTitles :
        calibGr[wType].Fit(calibModel,'MER+','',5,1000)
        calibGr[wType].GetFunction(calibModel.GetName()).SetRange(0,1000)
        calibGr[wType].GetFunction(calibModel.GetName()).SetLineColor(calibGr[wType].GetListOfGraphs().At(0).GetLineColor())

    #show results
    resCorrectionGr,resCalibGr=showCalibrationCurves(calibGr=calibGr,calibRanges=etaRanges,outDir=outDir,calibPostFix=calibPostFix)
    showResolutionCurves(resGr=resGr,outDir=outDir,calibPostFix=calibPostFix,model=0)

    #save all to file
    calibModelRes=ROOT.TF1('calibmodelres',"[0]*x*x+[1]*x+[2]",1.45,3.1)
    calibF=ROOT.TFile.Open('%s/calib_%s.root'%(outDir,calibPostFix),'RECREATE')
    for wType in weightTitles :
        calibGr[wType].Write()
        calibGr[wType].GetFunction(calibModel.GetName()).Write('%s_calib'%wType)
        for gr in resCalibGr[wType].GetListOfGraphs():
            gr.Fit(calibModelRes,'WMR+')
            gr.Write()            
        resCorrectionGr[wType].Write()
        resCorrectionGr[wType].Fit(calibModelRes,'WMR+')
        resCorrectionGr[wType].GetFunction(calibModelRes.GetName()).Write('%s_calib_res'%wType)
        for gr in resGr[wType].GetListOfGraphs():
            gr.Write()
    calibF.Close()


"""
steer 
"""
def main():
    
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i',      '--in' ,      dest='input',         help='Input file',                                                     default=None)
    parser.add_option('-w',      '--ws' ,      dest='wsUrl',         help='Workspace file',                                                 default=None)
    parser.add_option('--emCalib' ,            dest='emCalibUrl',    help='em calibration files (e.g. EE:calib_ee.root,HEF:calib_hef.root', default=None)
    parser.add_option('--noComp' ,             dest='noComp',        help='no attempt to correct for compensation',                         default=False, action="store_true")
    parser.add_option('--calib' ,              dest='calibUrl',      help='pion calibration file',                                          default=None)
    parser.add_option('--compWeights' ,        dest='compWeights',   help='file with software compensation weights',                        default=None)
    parser.add_option('--vetoTrackInt',        dest='vetoTrackInt',  help='flag if tracker interactions should be removed',                 default=False, action="store_true")
    parser.add_option('--vetoHEBLeaks',        dest='vetoHEBLeaks',  help='flag if HEB leaks are allowed',                                  default=False, action='store_true')
    parser.add_option('--noResCalib',          dest='noResCalib',    help='Don\'t run calibration fits in E/eta slices',                    default=False, action='store_true')
    parser.add_option('--noEE',                dest='noEE',          help='Assign weight 0 to EE',                                          default=False, action='store_true')
    parser.add_option('--noHEF',               dest='noHEF',         help='Assign weight 0 to HEF',                                         default=False, action='store_true')
    parser.add_option('--hebResp',             dest='hebRespUrl',    help='Location of the parameterization for HEB pi/e',                  default=None)
    parser.add_option('--hefResp',             dest='hefRespUrl',    help='Location of the parameterization for HEF and HEF pi/e',         default=None)
    parser.add_option('--eeResp',              dest='eeRespUrl',     help='Location of the parameterization for EE pi/e',                   default=None)
    parser.add_option('-v',      '--var' ,     dest='treeVarName',   help='Variable to use as energy estimator',                            default='edep_rec')
    (opt, args) = parser.parse_args()

    #check inputs    
    if opt.input is None and opt.wsUrl is None:
        parser.print_help()
        sys.exit(1)

    #basic ROOT customization
    customROOTstyle()
    ROOT.gSystem.Load( "libUserCodeHGCanalysis")

    #ROOT.gROOT.SetBatch(False)
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.RooMsgService.instance().setSilentMode(True);
    ROOT.RooMsgService.instance().getStream(0).removeTopic(ROOT.RooFit.Minimization);
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Minimization);
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.ObjectHandling);
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.DataHandling);
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Fitting);
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Plotting);
    ROOT.RooMsgService.instance().getStream(0).removeTopic(ROOT.RooFit.InputArguments);
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.InputArguments);
    ROOT.RooMsgService.instance().getStream(0).removeTopic(ROOT.RooFit.Eval);
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Eval);
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Integration);
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.NumIntegration);
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.NumIntegration);

    runCalibrationStudy(opt)
    
if __name__ == "__main__":
    sys.exit(main())
