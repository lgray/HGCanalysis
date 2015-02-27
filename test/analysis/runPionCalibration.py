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
            redData=ws.data('data').reduce('en>=%f && en<=%f && eta>=%f && eta<=%f'%(genEn_min,genEn_max,genEta_min,genEta_max))
            if redData.numEntries()<10 : continue

            iprof=len(all2DScatters)
            all2DScatters.append( ROOT.TH2F('hprof%d%d'%(ieta,ien),';%s energy/E_{beam};%s energy/E_{beam};Events'%(xaxis,yaxis),75,0,1.5,75,0,1.5) )
            all2DScatters[iprof].SetDirectory(0)
            all2DScatters[iprof].Sumw2()
            for ientry in xrange(0,redData.numEntries()):
                entryVars=redData.get(ientry)
                xval=entryVars.find('en_%s'%xaxis).getVal()/genEn_mean
                yval=0
                for iyvar in xrange(0,len(yAxisVars)):
                    yval+=entryVars.find('en_%s'%yAxisVars[iyvar]).getVal()/yAxisWeights[iyvar]
                yval/=genEn_mean
                all2DScatters[iprof].Fill(xval,yval)
    
            #show profile and fit
            p=canvas.cd(ieta+1)
            p.SetLogz()
            p.SetRightMargin(0.1)
            
            #find fit range from projection quantiles to avoid biases
            projY=all2DScatters[iprof].ProjectionX('tmpproj',2)
            probSum = array.array('d', [0.05,0.4])
            if len(yAxisVars)==2 : probSum = array.array('d', [0.1,0.7])
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
            redData=ws.data('data').reduce('en>=%f && en<=%f && eta>=%f && eta<=%f'%(genEn_min,genEn_max,genEta_min,genEta_max))
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
                xval=entryVars.find('en_%s'%xaxis).getVal()/genEn_mean
                yval=0
                for iyvar in xrange(0,len(yAxisVars)):
                    yval+=entryVars.find('en_%s'%yAxisVars[iyvar]).getVal()/yAxisWeights[iyvar]
                yval/=genEn_mean
                allSums[iprof].Fill(xval+yval)
                allWeightedSums[iprof].Fill(xval+yval/finalCombSlope)

            #show raw sum and after weighted combination
            pad=canvas.cd(ieta+1)
            allSums[iprof].Draw('hist')
            allSums[iprof].GetYaxis().SetRangeUser(0,allSums[iprof].GetMaximum()*2)
            allSums[iprof].SetTitle("sum")
            allWeightedSums[iprof].Draw('histsame')
            allWeightedSums[iprof].SetTitle("weighted sum")
            MyPaveText('#it{Energy=%d GeV, %3.1f<#eta<%3.1f}\\RMS/mean sum : %3.3f\\RMS/mean weighted sum : %3.3f'%
                       (genEn_mean,genEta_min,genEta_max,
                        allSums[iprof].GetRMS()/allSums[iprof].GetMean(),
                        allWeightedSums[iprof].GetRMS()/allWeightedSums[iprof].GetMean()),
                       0.15,0.95,0.5,0.7).SetTextSize(0.04)
            if ieta==0 :
                MyPaveText('#bf{CMS} #it{simulation}')
                leg=ROOT.TLegend(0.15,0.68,0.5,0.5)
                leg.SetFillStyle(0)
                leg.SetBorderSize(0)
                leg.SetTextFont(42)
                leg.SetTextSize(0.04)
                leg.AddEntry(allSums[iprof],allSums[iprof].GetTitle(),'f')
                leg.AddEntry(allWeightedSums[iprof],allWeightedSums[iprof].GetTitle(),'f')
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

    print '[computeCompensationWeights] will look at correlations between reconstructed energy and hit fraction (global) or shower density (local)'

    #wgtfunc       = ROOT.TF1('wgtfunc','[0]+[1]*pow(x,[2])',0,2)
    wgtfunc       = ROOT.TF1('wgtfunc','[0]+[1]*x',0,2)
    weightOptimGr = {'rho':[],'c':[]}
    aGr           = {'rho':[],'c':[]}
    for key in aGr:
        for ip in xrange(0,wgtfunc.GetNpar()):
            aGr[key].append( ROOT.TGraphErrors() )
            aGr[key][ip].SetMarkerStyle(20)
            aGr[key][ip].SetName('%s_swweights_%d'%(key,ip))


    canvas=ROOT.TCanvas('c','c',1200,800)
    all2DScatters={'c_EE':[],'c_HEF':[],'c_HEB':[],'rho_EE':[],'rho_HEF':[],'rho_HEB':[]}
    all2DProfiles={'c_EE':[],'c_HEF':[],'c_HEB':[],'rho_EE':[],'rho_HEF':[],'rho_HEB':[]}
    all1DProfiles={'c_EE':[],'c_HEF':[],'c_HEB':[],'rho_EE':[],'rho_HEF':[],'rho_HEB':[]}
    for ien in xrange(0,len(enRanges)):

        genEn_min=enRanges[ien][0]
        genEn_max=enRanges[ien][1]
        genEn_mean=0.5*(genEn_max+genEn_min)
        iprof1d=len(all1DProfiles['rho_EE'])

        #optimize in energy density ranges, based on the inclusive profile
        all1DProfiles['rho_EE'].append( ROOT.TH1F('rho_EE_hprof1d_%d'%ien,'%d GeV;Shower density [GeV/Volume];PDF'%genEn_mean,25,0,2) )
        all1DProfiles['rho_EE'][iprof1d].SetDirectory(0)
        all1DProfiles['rho_EE'][iprof1d].Sumw2()
        all1DProfiles['rho_HEF'].append( all1DProfiles['rho_EE'][iprof1d].Clone('rho_HEF_hprof1d_%d'%ien) )
        all1DProfiles['rho_HEF'][iprof1d].SetDirectory(0)
        all1DProfiles['rho_HEB'].append( all1DProfiles['rho_EE'][iprof1d].Clone('rho_HEB_hprof1d_%d'%ien) )
        all1DProfiles['rho_HEB'][iprof1d].SetDirectory(0)

        all1DProfiles['c_EE'].append( ROOT.TH1F('c_EE_hprof1d_%d'%ien,'%d GeV;C(10 MIP);PDF'%genEn_mean,25,0.5,2.0) )
        all1DProfiles['c_EE'][iprof1d].SetDirectory(0)
        all1DProfiles['c_EE'][iprof1d].Sumw2()
        all1DProfiles['c_HEF'].append( all1DProfiles['c_EE'][iprof1d].Clone( 'c_HEF_hprof1d_%d'%ien ) )
        all1DProfiles['c_HEF'][iprof1d].SetDirectory(0)
        all1DProfiles['c_HEB'].append( all1DProfiles['c_EE'][iprof1d].Clone( 'c_HEB_hprof1d_%d'%ien ) )
        all1DProfiles['c_HEB'][iprof1d].SetDirectory(0)

        #fill histograms and prepare for quantile computation
        rho_values,c_values=[],[]
        redIncData=ws.data('data').reduce('en>=%f && en<=%f'%(genEn_min,genEn_max))
        nEntries=redIncData.numEntries()
        for ientry in xrange(0,nEntries):
            entryVars=redIncData.get(ientry)
            for subDet in ['EE','HEF','HEB']:
                rhoval=entryVars.find('rho_%s'%subDet).getVal()
                cval=entryVars.find('c_%s'%subDet).getVal()
                if subDet=='HEF':
                    rho_values.append(rhoval)
                    c_values.append(cval)
                if rhoval>0 : all1DProfiles['rho_%s'%subDet][iprof1d].Fill(rhoval,1./nEntries)
                if cval>0 : all1DProfiles['c_%s'%subDet][iprof1d].Fill(cval,1./nEntries)

        #sums for weight optimization
        rhoRanges,       cRanges       = [], []
        prevRhoQuantile, prevCQuantile = 0,  0 
        for q in [5,25,50,75]:
            rhoQuantile, cQuantile = numpy.percentile(rho_values,q), numpy.percentile(c_values,q)
            if prevRhoQuantile!=rhoQuantile : rhoRanges.append([prevRhoQuantile,rhoQuantile])
            if prevCQuantile!=cQuantile     : cRanges.append([prevCQuantile,cQuantile])
            prevRhoQuantile, prevCQuantile = rhoQuantile, cQuantile 
        rhoRanges.append([prevRhoQuantile,ROOT.TMath.Min(numpy.percentile(rho_values,99),2)])
        cRanges.append([prevCQuantile,ROOT.TMath.Min(numpy.percentile(c_values,99),2)])
        print 'For E=%d GeV compensation weights will be optimised in the following en. density ranges'%genEn_mean
        print 'Density : ',rhoRanges
        print 'C       : ',cRanges

        #arrays for optimisation:  delta^2 = [ w Erec / Egen - 1]^2 
        # => [ Erec/Egen ] [ w Erec / Egen - 1 ] = 0
        # ~  a(w.b +c ) = 0 <=> w = - ca / ab
        # with ca = - [Erec/Egen] 
        #      ab =  [Erec/Egen]^2
        ca_rho_Sum, ab_rho_Sum = [0]*len(rhoRanges), [0]*len(rhoRanges)
        ca_c_Sum,   ab_c_Sum   = [0]*len(cRanges), [0]*len(cRanges)
        rho_values, c_values   = [], []
        for irho in xrange(0,len(rhoRanges)) : rho_values.append( [] )
        for ic in xrange(0,len(cRanges))     : c_values.append( [] )

        iprof=len(all2DScatters['rho_EE'])
        all2DScatters['rho_EE'].append( ROOT.TH2F('rho_EE_hprof2d_%d'%ien,';E_{rec} [GeV];Shower density [GeV/Volume];Events',30,0.5*genEn_mean,3*genEn_mean,25,0,2) )
        all2DScatters['rho_EE'][iprof].SetDirectory(0)
        all2DScatters['rho_EE'][iprof].Sumw2()
        all2DScatters['rho_HEF'].append( all2DScatters['rho_EE'][iprof].Clone('rho_HEF_hprof2d_%d'%ien) )
        all2DScatters['rho_HEF'][iprof].SetDirectory(0)
        all2DScatters['rho_HEB'].append( all2DScatters['rho_EE'][iprof].Clone('rho_HEB_hprof2d_%d'%ien) )
        all2DScatters['rho_HEB'][iprof].SetDirectory(0)
        all2DScatters['c_EE'].append( ROOT.TH2F('c_EE_hprof2d_%d'%ien,';E_{rec} [GeV];C(10 MIP);Events',30,0.5*genEn_mean,3*genEn_mean,25,0.5,2) )
        all2DScatters['c_EE'][iprof].SetDirectory(0)
        all2DScatters['c_EE'][iprof].Sumw2()
        all2DScatters['c_HEF'].append( all2DScatters['c_EE'][iprof].Clone('c_HEF_hprof2d_%d'%ien) )
        all2DScatters['c_HEF'][iprof].SetDirectory(0)
        all2DScatters['c_HEB'].append( all2DScatters['c_EE'][iprof].Clone('c_HEB_hprof2d_%d'%ien) )
        all2DScatters['c_HEB'][iprof].SetDirectory(0)

        for ieta in xrange(0,len(etaRanges)):
            genEta_min=etaRanges[ieta][0]
            genEta_max=etaRanges[ieta][1]
            
            #fill a new profile
            redData=ws.data('data').reduce('en>=%f && en<=%f && eta>=%f && eta<=%f'%(genEn_min,genEn_max,genEta_min,genEta_max))
            if redData.numEntries()<10 : continue

            for ientry in xrange(0,redData.numEntries()):
                entryVars=redData.get(ientry)

                egen=entryVars.find('en').getVal()                
                e_EE  = entryVars.find('en_EE').getVal()*ws.var('k_EE').getVal()
                e_HEF = entryVars.find('en_HEF').getVal()*ws.var('k_HEF').getVal()
                e_HEB = entryVars.find('en_HEB').getVal()*ws.var('k_HEB').getVal()
                e_tot = e_EE+e_HEF+e_HEB

                for subDet in ['EE','HEF','HEB']:

                    rhoval=entryVars.find('rho_%s'%subDet).getVal()
                    all2DScatters['rho_%s'%subDet][iprof].Fill(e_tot,rhoval)

                    cval=entryVars.find('c_%s'%subDet).getVal()
                    all2DScatters['c_%s'%subDet][iprof].Fill(e_tot,cval)

                    if subDet!='HEF': continue

                    for irho in xrange(0,len(rhoRanges)):
                        if rhoval<rhoRanges[irho][0] or rhoval>rhoRanges[irho][1] : continue
                        rho_values[ irho ].append( rhoval )
                        ca_rho_Sum[ irho ] += (e_tot/egen)
                        ab_rho_Sum[ irho ] += (e_tot/egen)*(e_tot/egen)
                        
                    for ic in xrange(0,len(cRanges)):
                        if cval<cRanges[ic][0] or cval>cRanges[ic][1] : continue
                        c_values[ ic ].append( cval )
                        ca_c_Sum[ ic ] += (e_tot/egen)
                        ab_c_Sum[ ic ] += (e_tot/egen)*(e_tot/egen)

        #now optimize weights for this energy
        for key in weightOptimGr:
            weightOptimGr[key].append( ROOT.TGraph() )
            weightOptimGr[key][ien].SetName('%s_sweights_%d'%(key,ien))
            weightOptimGr[key][ien].SetTitle('%d GeV'%genEn_mean)
            weightOptimGr[key][ien].SetMarkerStyle(20)
            if key=='rho':
                for irho in xrange(0,len(rhoRanges)):
                    if ab_rho_Sum[irho]==0: continue
                    np=weightOptimGr[key][ien].GetN()
                    weightOptimGr[key][ien].SetPoint(np,numpy.mean(rho_values[irho]),ca_rho_Sum[irho]/ab_rho_Sum[irho])
            else:
                for ic in xrange(0,len(cRanges)):
                    if ab_c_Sum[ic]==0: continue
                    np=weightOptimGr[key][ien].GetN()
                    weightOptimGr[key][ien].SetPoint(np,numpy.mean(c_values[ic]),ca_c_Sum[ic]/ab_c_Sum[ic])
                    
            #save evolution of the parameters for this energy       
            weightOptimGr[key][ien].Fit(wgtfunc,'MQR+')
            for ip in xrange(0,wgtfunc.GetNpar()):
                np=aGr[key][ip].GetN()
                aGr[key][ip].SetPoint(np,genEn_mean,wgtfunc.GetParameter(ip))
                aGr[key][ip].SetPointError(np,0,wgtfunc.GetParError(ip))
        
        #show profile
        canvas.Clear()
        canvas.Divide(3,2)   
        ipad=0
        for var in ['rho','c']:
            for subDet in ['EE','HEF','HEB']: 
                ipad+=1
                p=canvas.cd(ipad)
                p.SetLogz()
                p.SetRightMargin(0.1)
                key='%s_%s'%(var,subDet)
                all2DProfiles[key].append( all2DScatters[key][iprof].ProfileX('%s_prof'%all2DScatters[key][iprof].GetName()) )
                all2DProfiles[key][iprof].SetMarkerStyle(20)
                all2DScatters[key][iprof].Draw('colz')
                all2DScatters[key][iprof].GetZaxis().SetTitleOffset(-0.5)            
                all2DScatters[key][iprof].GetYaxis().SetTitleOffset(1.2)            
                all2DProfiles[key][iprof].Draw('e1same')
                MyPaveText('[ %s ]'%subDet,0.8,0.95,0.95,0.99)
                if ipad>1: continue
                MyPaveText('#bf{CMS} #it{simulation} Energy=%d'%genEn_mean)
        canvas.Modified()
        canvas.Update()
        canvas.SaveAs('%s/showerdens_%d_profile.png'%(outDir,ien))
       
    #compare shower densities for different events
    canvas.Clear()
    canvas.Divide(3,2)   
    ipad=0
    for var in ['rho','c']:
        for subDet in ['EE','HEF','HEB']:
            key='%s_%s'%(var,subDet)
            ipad+=1
            p=canvas.cd(ipad)
            p.Clear()
            p.SetLogy(True)
            for iprof1d in xrange(0,len(all1DProfiles[key])):
                all1DProfiles[key][iprof1d].SetLineColor(45-2*iprof1d)
                all1DProfiles[key][iprof1d].SetMarkerColor(45-2*iprof1d)
                all1DProfiles[key][iprof1d].SetMarkerStyle(1)
                all1DProfiles[key][iprof1d].SetLineWidth(2)
                if iprof1d==0 : 
                    all1DProfiles[key][iprof1d].Draw('hist')
                    all1DProfiles[key][iprof1d].GetYaxis().SetRangeUser(1e-3,1.0)
                    all1DProfiles[key][iprof1d].GetYaxis().SetTitleOffset(1.2)
                else:
                    all1DProfiles[key][iprof1d].Draw('histsame')

            MyPaveText('[ %s ]'%subDet,0.8,0.95,0.95,0.99)

            if ipad>1: continue
            MyPaveText('#bf{CMS} #it{simulation}')
            leg=p.BuildLegend(0.6,0.6,0.9,0.94)
            leg.SetFillStyle(0)
            leg.SetBorderSize(0)
            leg.SetTextFont(42)
            leg.SetTextSize(0.03)

    canvas.Modified()
    canvas.Update()
    canvas.SaveAs('%s/showerdens.png'%outDir)

 
    #compare evolution
    npar=wgtfunc.GetNpar()
    canvas.Clear()
    canvas.SetWindowSize(400*npar,400)
    canvas.SetCanvasSize(400*npar,400)
    canvas.Divide(npar,1)
    for ip in xrange(0,wgtfunc.GetNpar()):
        p=canvas.cd(ip+1)
        p.SetLogy(False)
        key='rho'
        aGr[key][ip].Draw('ap')
        aGr[key][ip].GetXaxis().SetTitle('Energy [GeV]')
        aGr[key][ip].GetYaxis().SetTitle('Weight function parameter # %d'%ip)
        aGr[key][ip].GetYaxis().SetTitleOffset(1.4)
        aFunc=ROOT.TF1('%s_evfunc_%d'%(key,ip),'[0]+[1]*(1-exp(-x*[2]))',0,1000)
        aFunc.SetParLimits(0,-100,100)
        aFunc.SetParLimits(1,-50,50)
        aFunc.SetParLimits(2,0.01,100)
        aGr[key][ip].Fit(aFunc,'MR+')
        if ip==0 : MyPaveText('#bf{CMS} #it{simulation}')
        MyPaveText('Parameter evolution\\p_{%d}(E)=q_{1}+q_{2}(1-e^{-q_{3}E})\\q_{1}=%3.3f#pm%3.3f\\q_{2}=%3.3f#pm%3.3f\\q_{3}=%3.3f#pm%3.3f'
                   %(ip,aFunc.GetParameter(0),aFunc.GetParError(0),
                     aFunc.GetParameter(1),aFunc.GetParError(1),
                     aFunc.GetParameter(2),aFunc.GetParError(2)),
                   0.4,0.3,0.9,0.6).SetTextSize(0.035)
    canvas.Modified()
    canvas.Update()
    canvas.SaveAs('%s/swcompweights.png'%(outDir))
    canvas.SaveAs('%s/swcompweights.C'%(outDir))


    #save weights and weight function to file
    swcompUrl='%s/swcompweights.root'%outDir
    fOut=ROOT.TFile(swcompUrl,'RECREATE')
    fOut.cd()
    for key in aGr : 
        for gr in aGr[key] : gr.Write()
        for w in weightOptimGr[key] : w.Write()
    wgtfunc.Write()
    fOut.Close()
    
    #all done
    return swcompUrl


"""
Adapts the workspace for pion calibration
"""
def adaptWorkspaceForPionCalibration(opt,outDir):

    wsUrl             = opt.wsUrl

    #prepare workspace (if needed) and output
    if wsUrl is None :

        #electromagnetic energy scale of each section
        emCalibMap={'EE':1.0,'HEF':1.0,'HEB':1.0}
        if opt.emCalibUrl :
            for iurl in opt.emCalibUrl.split(','):
                subDet,url=iurl.split(':')
                print 'Replacing default energy scale in %s with calibration from %s'%(subDet,url)

                #discard any offset and only use the slope
                calibF=ROOT.TFile.Open(url)
                emCalibMap[subDet]=1./calibF.Get('lambda_calib').GetParameter(0)
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
    enRanges  = [[9,11],[19,21],[39,41],[49,51],[74,75],[99,101],[124,126],[174,176],[249,251],[399,401]]

    #determine how to combine HEF and HEB
    hefhebCombSlope, hefhebCombSlope_err = 1.0, 0.0
    try:
        hefhebFin=ROOT.TFile.Open(opt.hefhebCombUrl)
        hefhebCombSlope, hefhebCombSlope_err = hefhebFin.Get('HEFHEB_combfunc').GetParameter(0), hefhebFin.Get('HEFHEB_combfunc').GetParError(0)
        hefhebFin.Close()
    except:
        print 'Will compute HEF vs HEB combination slope'
        hefhebCombSlope, hefhebCombSlope_err = computeSubdetectorWeights(enRanges=enRanges,etaRanges=etaRanges,ws=ws,xaxis='HEF',yaxis='HEB',outDir=outDir)
    if opt.vetoHEB:
        print 'vetoeing HEB in the end...'
        hefhebCombSlope=99999999.
    ws.factory('k_HEB[%f]'%(1.0/hefhebCombSlope))
    #ws.factory('k_HEB[1.0]')

    #determine how to combine EE and HE(F+B)
    ehCombSlope, ehCombSlope_err = 1.0, 0.0
    ws.factory('k_EE[%f]'%(1.0-int(opt.noEE)))
    if opt.noEE==False:
        try:
            ehFin=ROOT.TFile.Open(opt.ehCombUrl)
            ehCombSlope, ehCombSlope_err = ehFin.Get('EEHEFHEB_combfunc').GetParameter(0), ehFin.Get('EEHEFHEB_combfunc').GetParError(0)
            ehFin.Close()
        except:
            print 'Will compute EE vs HE(F+B) combination slope using %f coefficient for HEB'%ws.var('k_HEB').getVal()
            ehCombSlope, ehCombSlope_err = computeSubdetectorWeights(enRanges=enRanges,etaRanges=etaRanges,ws=ws,xaxis='EE',yaxis='HEF+HEB/%3.4f'%hefhebCombSlope,outDir=outDir)
    ws.factory('k_HEF[%f]'%(1.0/ehCombSlope))
    ws.var('k_HEB').setVal(ws.var('k_HEB').getVal()/ehCombSlope)

    if opt.noComp:
        ws.var('k_HEF').setVal(1.0)
        ws.var('k_HEB').setVal(1.0)

    print 'Will use the following combination of sub-detectors %3.4f x EE + %3.4f x HEF + %3.4f x HEB'%(ws.var('k_EE').getVal(),ws.var('k_HEF').getVal(),ws.var('k_HEB').getVal())
 
    #read sw compensation weights
    swCompParams=[]
    if opt.compWeights is None:
        opt.compWeights=computeCompensationWeights(enRanges,etaRanges,ws,outDir)
    swF=ROOT.TFile.Open(opt.compWeights)
    print 'Reading out compensation weights from %s'%(opt.compWeights)
    swwgtfunc=swF.Get('wgtfunc')
    for ip in xrange(0,swwgtfunc.GetNpar()):
        swCompParams.append( swF.Get('rho_swweights_%d'%ip).GetFunction('rho_evfunc_%d'%ip) )
        #swCompParams.append( swF.Get('rho_swweights_%d'%ip) )
    swF.Close()

    #read calibrations from file, if available
    weightTitles={'simple':'Simple sum ','gc':'Global compensation','lc':'Local compensation'}
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
            #calibMapRes[wType]=[]
            #for ietaRange in xrange(0,len(etaRanges)):
            #    calibMapRes[wType].append( calibF.Get('calib_%d_%s_res'%(ietaRange,wType)).Clone() )
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

        #iEtaRange=-1
        #for etaRange in etaRanges:
        #    iEtaRange+=1
        #    if ws.var('eta')<etaRange[0] or ws.var('eta')>etaRange[1]: continue 
        #    break

        enEstimators={}

        #simple sum
        e_EE    = entryVars.find('en_EE').getVal()*ws.var('k_EE').getVal()
        e_HEF   = entryVars.find('en_HEF').getVal()*ws.var('k_HEF').getVal()
        e_HEB   = entryVars.find('en_HEB').getVal()*ws.var('k_HEB').getVal()
        e_tot   = e_EE + e_HEF + e_HEB
        enEstimators['simple']=e_tot

        #global compensation weights
        c_EE    = entryVars.find('c_EE').getVal()
        c_HEF   = entryVars.find('c_HEF').getVal()
        c_HEB   = entryVars.find('c_HEB').getVal()
        e_c_tot = e_EE + c_HEF*e_HEF + e_HEB
        enEstimators['gc']=e_c_tot

        #software compensation weight
        e_rho_tot = e_tot
        rho_HEF = entryVars.find('rho_HEF').getVal()
        if e_HEF > 0 and rho_HEF>0:
            for ip in xrange(0,swwgtfunc.GetNpar()):
                swwgtfunc.SetParameter(ip,swCompParams[ip].Eval(e_tot))
            swWgt=swwgtfunc.Eval(ROOT.TMath.Min(rho_HEF,2.0))
            e_rho_tot = e_EE + swWgt*e_HEF + e_HEB
        enEstimators['lc']=e_rho_tot

        #now apply calibration
        for wType in weightTitles :

            ienVal=enEstimators[wType]

            calib_offset, calib_slope = 0.0, 1.0
            if wType in calibMap:
                calib_offset=calibMap[wType].GetParameter(1)
                calib_slope=calibMap[wType].GetParameter(0)
            
            calib_residual = 0.0
            if wType in calibMapRes:
                calib_residual=calibMapRes[wType].Eval(ROOT.TMath.Abs(ws.var('eta').getVal()))
                #calib_residual=calibMapRes[wType][iEtaRange].Eval(ienVal)/100
                
            #calibrated energy estimator
            ienVal=((ienVal-calib_offset)/calib_slope)*(1-calib_residual)
                
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
                ws.factory('RooCBShape::resol_%s(%s,mean_%s[%f,%f,%f],sigma_%s[%f,%f,%f],alpha_%s[-3.0,-10.0,10.0],n_%s[2,1,3])'%
                           (fitName,vName,
                            fitName,v_mean,v_min,v_max,
                            fitName,v_sigma,v_sigma*0.001, v_sigma*1.5,
                            fitName,
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
    calibModel=ROOT.TF1('calibmodel',"[0]*x+[1]",0,800)
    #calibModel=ROOT.TF1('calibmodel',"x>100 ? [0]*x+[1] : [2]*x*x+[3]*x+[4]",0,1000)
    calibModel.SetLineWidth(1)
    for wType in weightTitles :
        calibGr[wType].Fit(calibModel,'MER+')
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
    parser.add_option('--hefhebComb',          dest='hefhebCombUrl', help='Location of the parameterization for HEB/HEF combination',       default=None)
    parser.add_option('--noEE',                dest='noEE',          help='Assign weight 0 to EE',                                          default=False, action='store_true')
    parser.add_option('--vetoHEB',             dest='vetoHEB',       help='Veto HEB',                                                       default=False, action='store_true')
    parser.add_option('--ehComb',              dest='ehCombUrl',     help='Location of the parameterization for EE+HE(F+B) combination',    default=None)
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
