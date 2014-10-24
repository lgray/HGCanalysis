#!/usr/bin/env python

import ROOT
import numpy as numpy
from array import array
import io,os,sys
import json
import optparse
import commands
from UserCode.HGCanalysis.PlotUtils import *

"""
Converts all to a workspace and returns optimized weights
"""
def prepareWorkspace(url,integRanges,treeVarName,deriveOffset,outUrl):
    
    #prepare the workspace
    ws=ROOT.RooWorkspace("w")
    dsVars=ROOT.RooArgSet( ws.factory('eta[1.5,1.45,3.1]'), ws.factory('en[0,0,9999999999]'), ws.factory('phi[0,-3.2,3.2]') )
    for ireg in xrange(0,len(integRanges)): dsVars.add( ws.factory('edep%d[0,0,99999999]'%ireg) )
    getattr(ws,'import')( ROOT.RooDataSet('data','data',dsVars) )

    #optimization with linear regression
    optimVec    = numpy.zeros( len(integRanges)+int(deriveOffset) )
    optimMatrix = numpy.zeros( (len(integRanges)+int(deriveOffset), len(integRanges)+int(deriveOffset) ) )

    #read all to a RooDataSet
    fin=ROOT.TFile.Open(url)
    HGC=fin.Get('analysis/HGC')
    for entry in xrange(0,HGC.GetEntriesFast()+1):
        HGC.GetEntry(entry)
        genEn=HGC.genEn
        genEta=ROOT.TMath.Abs(HGC.genEta)
        genPhi=HGC.genPhi
        if genEta<1.4 or genEta>3.0 : continue
        ws.var('eta').setVal(genEta)
        ws.var('en').setVal(genEn)
        ws.var('phi').setVal(genPhi)
        newEntry=ROOT.RooArgSet( ws.var('eta'), ws.var('en'), ws.var('phi') )

        #get the relevant energy deposits and add new row
        edeps=[]
        for ireg in xrange(0,len(integRanges)):
            totalEnInIntegRegion=0
            for ilayer in xrange(integRanges[ireg][0],integRanges[ireg][1]+1):
                totalEnInIntegRegion=totalEnInIntegRegion+(getattr(HGC,treeVarName))[ilayer-1]
            edeps.append(totalEnInIntegRegion)
            ws.var('edep%d'%ireg).setVal(totalEnInIntegRegion)
            newEntry.add(ws.var('edep%d'%(ireg)))

        #conversion veto
        #if edeps[0]>20 : continue
        ws.data('data').add( newEntry )

        #fill the optmization matrix and vector for a subset of events (low energy, low eta)

        #include weights and offset
        if deriveOffset:
            for ie in xrange(0,len(edeps)):
                optimVec[ie]=optimVec[ie]+edeps[ie]/genEn
                for je in xrange(0,len(edeps)): 
                    optimMatrix[ie][je]=optimMatrix[ie][je]+edeps[ie]*edeps[je]/(genEn*genEn)
                optimMatrix[len(edeps)][ie]=optimMatrix[len(edeps)][ie]+edeps[ie]/(genEn*genEn)
                optimMatrix[ie][len(edeps)]=optimMatrix[ie][len(edeps)]+edeps[ie]/(genEn*genEn)
            optimVec[len(edeps)]=optimVec[len(edeps)]+1.0/genEn        
            optimMatrix[len(edeps)][len(edeps)]=optimMatrix[len(edeps)][len(edeps)]+1.0/(genEn*genEn)
            
        #include only weights
        else:
            for ie in xrange(0,len(edeps)):
                optimVec[ie]=optimVec[ie]+edeps[ie]/genEn
                for je in xrange(0,len(edeps)): 
                    optimMatrix[ie][je]=optimMatrix[ie][je]+edeps[ie]*edeps[je]/(genEn*genEn)


    fin.Close()

    #all done, write to file
    ws.writeToFile(outUrl,True)
    print 'Created the analysis RooDataSet with %d events, stored @ %s'%(ws.data('data').numEntries(),outUrl)

    #finalize weight optimization (solve linear equation system)
    optimWeights=numpy.linalg.solve(optimMatrix,optimVec)
    return optimWeights.tolist()

    

"""
runs the fits to the calibration and resolution
"""
def showEMcalibrationResults(methodNames,calibFunc,resFunc,outDir,isHadronic=False):

    nMethods=len(methodNames)
    nEta=len(resFunc)/nMethods

    #resolSummary=[]

    #
    #resolution
    #
    resolModel=ROOT.TF1('resolmodel',"sqrt([0]*[0]/x+[1]*[1])",0,1000)
    #resolModel=ROOT.TF1('resolmodel',"sqrt([0]*x*x+[1]+[2]*x*x*x*x)",0,1)
    resolModel.SetParameter(0,0.2);
    resolModel.SetParLimits(0,0,2);
    resolModel.SetParameter(1,0);
    resolModel.SetParLimits(1,0,1.0);
    resolModel.SetLineWidth(1)

    c=ROOT.TCanvas('cresol','cresol',500,500)
    c.SetTopMargin(0.05)
    c.SetRightMargin(0.02)
    c.SetLeftMargin(0.15)
    leg=ROOT.TLegend(0.6,0.6,0.9,0.95)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.03)
    c.SetLogx()
    resolgr=[]
    for imethod in xrange(0,nMethods):
        resolgr.append(ROOT.TMultiGraph())
        for ieta in xrange(0,nEta): resolgr[imethod].Add( resFunc[ieta+imethod*nEta], 'p' )
        if imethod==0 : 
            resolgr[imethod].Draw('a')
            resolgr[imethod].GetXaxis().SetTitle("E [GeV]")
            resolgr[imethod].GetYaxis().SetTitle("Relative energy resolution") 
            resolgr[imethod].GetXaxis().SetLabelSize(0.035)
            resolgr[imethod].GetYaxis().SetLabelSize(0.035)
            resolgr[imethod].GetXaxis().SetTitleSize(0.045)
            resolgr[imethod].GetYaxis().SetTitleSize(0.045)
            resolgr[imethod].GetYaxis().SetTitleOffset(1.2)
            resolgr[imethod].GetXaxis().SetMoreLogLabels(True)
            if isHadronic : resolgr[imethod].GetYaxis().SetRangeUser(0,1.0)
            else          : resolgr[imethod].GetYaxis().SetRangeUser(0,0.2)
        else : resolgr[imethod].Draw('p')

        #fit and add to legend
        resolgr[imethod].Fit(resolModel,'MER+')
        resolgr[imethod].GetFunction(resolModel.GetName()).SetLineColor(resFunc[imethod*nEta].GetLineColor())
        sigmaStoch=resolModel.GetParameter(0)
        sigmaStochErr=resolModel.GetParError(0)
        sigmaConst=resolModel.GetParameter(1)
        sigmaConstErr=resolModel.GetParError(1)
        leg.AddEntry(resFunc[imethod*nEta],"#splitline{#bf{[ %s ]}}{#frac{#sigma}{E} #propto #frac{%3.4f}{#sqrt{E}} #oplus %3.4f}"
                     %(methodNames[imethod],sigmaStoch,sigmaConst),
                     "l")
        #resolSummary.append( [sigmaStoch,sigmaStochErr,sigmaConst,sigmaConstErr] )        
    leg.Draw()
    MyPaveText('#bf{CMS} #it{simulation}')
    c.Modified()
    c.Update()
    for ext in ['png','pdf','C'] : c.SaveAs('%s/resolution.%s'%(outDir,ext))
    
    #
    #calibration
    #
    calibModel=ROOT.TF1('calibmodel',"[0]*x+[1]",0,800)
    
    c=ROOT.TCanvas('ccalib','ccalib',1200,500)
    c.SetTopMargin(0.05)
    c.SetRightMargin(0.5)
    leg=ROOT.TLegend(0.52,0.1,0.99,0.2+0.085*len(calibFunc)/3)
    leg.SetNColumns(3)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.025)
    for ir in xrange(0,len(calibFunc)):
        
        if ir==0:
            calibFunc[ir].Draw('ap')
            calibFunc[ir].GetXaxis().SetTitle("Energy [GeV]")
            calibFunc[ir].GetYaxis().SetTitle("Reconstructed energy")
            calibFunc[ir].GetXaxis().SetLabelSize(0.04)
            calibFunc[ir].GetYaxis().SetLabelSize(0.04)
            calibFunc[ir].GetXaxis().SetTitleSize(0.05)
            calibFunc[ir].GetYaxis().SetTitleSize(0.05)
            calibFunc[ir].GetYaxis().SetTitleOffset(1.3)
        else:
            calibFunc[ir].Draw('p')
        
        if ir%nEta==0 :
            calibFunc[ir].Fit(calibModel,'MER+')
            calibFunc[ir].GetFunction(calibModel.GetName()).SetLineColor(calibFunc[ir].GetLineColor())
        else     :
            calibFunc[ir].Fit(calibModel,'MER0+')
        
        slope=calibModel.GetParameter(0)
        slopeErr=calibModel.GetParError(0)
        offset=calibModel.GetParameter(1)
        offsetErr=calibModel.GetParError(1)
        leg.AddEntry(calibFunc[ir],
                     "#splitline{[#scale[0.8]{#bf{#it{%s}}}]}{#hat{E} = %3.1f#timesE_{beam} + %3.1f}"
                     %(calibFunc[ir].GetTitle(),slope,offset),
                     "fp")

        #resolSummary[ir].append(slope)
        #resolSummary[ir].append(slopeErr)
        #resolSummary[ir].append(offset)
        #resolSummary[ir].append(offsetErr)
        
    leg.Draw()
    MyPaveText('#bf{CMS} #it{simulation}')

    c.Modified()
    c.Update()
    for ext in ['png','pdf','C'] : c.SaveAs('%s/calibration.%s'%(outDir,ext))

    #return resolSummary


"""
"""
def runResolutionStudy(url,treeVarName='edep_thr0',isHadronic=False,deriveOffset=False,wsOutUrl=None):
    
    #init weights
    integRanges = [[1,1],[2,11],[12,21],[22,30],[31,31],[32,42],[43,54]]
    defWeights  = [0.08,  0.62,   0.81, 1.19,   3.58,   2.72,   2.31 ]
    if isHadronic : defWeights = [0.01, 0.036, 0.042, 0.055, 0.34, 0.25, 0.21 ] 

    #integRanges  = [[1,1]]
    #defWeights   = [0.08]
    #integRanges += [[2,2],[3,3],[4,4],[5,5],[6,6],[7,7],[8,8],[9,9],[10,10],[11,11]]
    #defWeights  += [0.92, 0.60,  0.57, 0.60,  0.57, 0.60,  0.57, 0.60,  0.57,   0.60]
    #integRanges += [[12,12],[13,13],[14,14],[15,15],[16,16],[17,17],[18,18],[19,19],[20,20],[21,21]]
    #defWeights  += [0.87,   0.79,   0.87,   0.79,   0.87,   0.79,   0.87,   0.79,   0.87,   0.79]
    #integRanges += [[22,22],[23,23],[24,24],[25,25],[26,26],[27,27],[28,28],[29,29],[30,30]]
    #defWeights  += [1.27,   1.20,   1.27,   1.20,   1.27,   1.20,   1.27,   1.20,   1.27]
    #integRanges += [[31,31], [32,42], [43,54]]
    #defWeights  += [3.58,    2.72,    2.31]
    #if isHadronic:
    #    defWeights   = [0.3]
    #    defWeights  += [1.0,   1.0,    0.5,    1.0,    0.5,    1.0,    0.5,    1.0,    0.5,   1.0]
    #    defWeights  += [0.8,   1.1,    0.8,    1.1,    0.8,    1.1,    0.8,    1.1,    0.8,   1.1]
    #    defWeights  += [1.1,   1.4,    1.1,    1.4,    1.1,    1.4,    1.1,    1.4,    1.1]
    #    defWeights  += [7.6,   5.6,    4.8]
    
    #calibration constant 
    defWeights += [0]

    #dummy initialization
    optimWeights = defWeights

    #prepare output
    outDir=url.replace('.root','')
    os.system('mkdir -p '+outDir)

    #get workspace
    if wsOutUrl is None:
        wsOutUrl='%s/workspace.root'%outDir
        optimWeights=prepareWorkspace(url=url,integRanges=integRanges,treeVarName=treeVarName,outUrl=wsOutUrl,deriveOffset=deriveOffset)
        if not deriveOffset : optimWeights.append(0)
    wsOutF=ROOT.TFile.Open(wsOutUrl)
    ws=wsOutF.Get('w')
    wsOutF.Close()

    #output weights to a file
    calibrationData={}
    calibrationData['IntegrationRanges'] = [ {'first':fLayer, 'last':lLayer} for fLayer,lLayer in integRanges ]
    calibrationData['DefaultWeights']    = [ item for item in defWeights ]
    calibrationData['OptimWeights']      = [ item for item in optimWeights ]
    with io.open('%s/weights.dat'%outDir, 'w', encoding='utf-8') as f: f.write(unicode(json.dumps(calibrationData, sort_keys = True, ensure_ascii=False, indent=4)))

    #prepare energy estimators
    funcArgs='{edep0'
    rawEnFunc,weightEnFunc,optimWeightEnFunc='edep0','%f*edep0'%defWeights[0],'%f*edep0'%optimWeights[0]
    for ireg in xrange(1,len(integRanges)) :
        funcArgs          += ',edep%d'%(ireg)
        rawEnFunc         += '+edep%d'%(ireg)
        weightEnFunc      += '+%f*edep%d'%(defWeights[ireg],ireg)
        optimWeightEnFunc += '+%f*edep%d'%(optimWeights[ireg],ireg)
    weightEnFunc      += '+%f'%defWeights[len(defWeights)-1]
    optimWeightEnFunc += '+%f'%optimWeights[len(optimWeights)-1]
    funcArgs=funcArgs+'}'

    baseVar='E'
    if treeVarName.find('nhits')>=0 : baseVar='N'
    vars=[
        [ws.data('data').addColumn(ws.factory("RooFormulaVar::rawEnFunc('"+rawEnFunc+"',"+funcArgs+")" )),                '#Sigma %s_{i}'%baseVar],
        [ws.data('data').addColumn(ws.factory("RooFormulaVar::weightEnFunc('"+weightEnFunc+"',"+funcArgs+")")),           '#Sigma w_{i}%s_{i}'%baseVar ],
        [ws.data('data').addColumn(ws.factory("RooFormulaVar::optimWeightEnFunc('"+optimWeightEnFunc+"',"+funcArgs+")")), '#Sigma w^{optim}_{i}%s_{i}'%baseVar]
        ]
    enEstimatorsSet=ROOT.RooArgSet(ws.var('eta'),ws.var('en'), ws.var('phi'))
    for v in vars:
        vName=v[0].GetName().replace('Func','')
        enEstimatorsSet.add( ws.factory('%s[0,0,999999999]'%vName) )
    getattr(ws,'import')( ROOT.RooDataSet('fitdata','fitdata',enEstimatorsSet) )
    methodNames=[]
    for i in xrange(0,len(vars)): methodNames.append( vars[i][1] )

    #create the fit dataset (skip direct usage of RooFormulaVar in a fit) and store final value
    for ientry in xrange(0,ws.data('data').numEntries()):
        entryVars=ws.data('data').get(ientry)
        for baseVar in ['eta','en','phi']: ws.var(baseVar).setVal( entryVars.find(baseVar).getVal() )
        newEntry=ROOT.RooArgSet( ws.var('eta'), ws.var('en'), ws.var('phi') )
        for v in vars:
            vName=v[0].GetName().replace('Func','')
            ws.var(vName).setVal( v[0].getVal() )
            newEntry.add( ws.var(vName) )
        ws.data('fitdata').add( newEntry )
    
    #first run the calibration
    varCtr=0
    calibFunc=[]
    resFunc=[]
    calibModel=ROOT.TF1('calibmodel',"[0]*x+[1]",0,800)
    c=ROOT.TCanvas('c','c',1400,800)
    nSigmasToFit=2.8

    etaRanges=[[1.5,2.0],[2.0,2.5],[2.5,2.9]]
    enRanges  =[[4,6],[19,21],[29,31],[49,51],[74,76],[99,101],[149,151],[249,251],[499,501]]
    #etaRanges=[[1.5,2.0]]
    #enRanges  =[[4,6]]
    for var,varTitle in vars:
        varCtr+=1
        vName=var.GetName().replace('Func','')
        outliersPhi=ROOT.TH1F('outliersphi',';#phi [rad];Events',100,0,3.2)

        #run calibration and resolution fits in different rapidity ranges
        for iEtaRange in xrange(0,len(etaRanges)):
            etaMin=etaRanges[iEtaRange][0]
            etaMax=etaRanges[iEtaRange][1]
            avgEta=0.5*(etaMax+etaMin)
            nv=len(calibFunc)
            calibFunc.append(ROOT.TGraphErrors())
            calibFunc[nv].SetMarkerStyle(20+iEtaRange)
            calibFunc[nv].SetTitle('%s %3.1f<#eta<%3.1f'%(varTitle,etaMin,etaMax))
            calibFunc[nv].SetFillStyle(0)
            calibFunc[nv].SetMarkerColor(varCtr+1)
            calibFunc[nv].SetLineColor(varCtr+1)
            calibFunc[nv].SetName('calib_%s_%d'%(vName,iEtaRange))
            resFunc.append(calibFunc[nv].Clone('resol_%s_%d'%(vName,iEtaRange)))

            fitCtr=0
            c.Clear()
            c.Divide(len(enRanges)/3+2,2)
            ipad=0
            for enMin,enMax in enRanges:
                fitCtr=fitCtr+1
                postfix='fit%d%d%d'%(fitCtr,nv,iEtaRange)
                print enMin,enMax,etaMin,etaMax
                redData=ws.data('fitdata').reduce('en>=%f && en<=%f && eta>=%f && eta<=%f'%(enMin, enMax, etaMin, etaMax))
                if redData.numEntries()<10 : continue
                
                #generator level information
                fitDataMean,   fitDataSigma   = redData.mean(ws.var(vName)), redData.sigma(ws.var(vName))
                fitDataEnMean, fitDataEnSigma = redData.mean(ws.var('en')),  redData.sigma(ws.var('en'))

                #define PDF and ranges to fit/show
                xaxisMin,xaxisMax=fitDataMean-5*fitDataSigma,fitDataMean+5*fitDataSigma
                xfitMin, xfitMax=fitDataMean-nSigmasToFit*fitDataSigma,fitDataMean+nSigmasToFit*fitDataSigma
                if fitDataSigma/fitDataMean > 0.2:
                    xaxisMin, xaxisMax = fitDataMean*0.5, fitDataMean*2
                    xfitMin,  xfitMax = fitDataMean*0.9, fitDataMean*1.1
                ws.var(vName).setRange(postfix,xaxisMin,xaxisMax)
                ws.var(vName).setRange('fit'+postfix,xfitMin,xfitMax)
                #ws.factory('Gaussian::g_%s(%s,mean_%s[%f,%f,%f],sigma_%s[%f,%f,%f])'%
                #           (postfix,vName,
                #            postfix,fitDataMean,fitDataMean-25*fitDataSigma,fitDataMean+25*fitDataSigma,
                #            postfix,fitDataSigma,fitDataSigma*0.001,fitDataSigma*2)
                #    )
                ws.factory('RooCBShape::g_%s(%s,mean_%s[%f,%f,%f],sigma_%s[%f,%f,%f],alpha_%s[0.001.,0.0001,20.0],n_%s[1])'%
                           (postfix,vName,
                            postfix,fitDataMean,fitDataMean-25*fitDataSigma,fitDataMean+25*fitDataSigma,
                            postfix,fitDataSigma,fitDataSigma*0.001,fitDataSigma*2,
                            postfix,postfix)
                           )

                #fit
                fres=ws.pdf('g_%s'%postfix).fitTo( redData, ROOT.RooFit.Range('fit'+postfix), ROOT.RooFit.Save(True) )
                meanFit       = ws.var('mean_%s'%postfix).getVal()
                sigmaFit       = ws.var('sigma_%s'%postfix).getVal()
                meanFitError  = ws.var('mean_%s'%postfix).getError()
                if meanFit<0: continue

                #save results
                np=calibFunc[nv].GetN()
                calibFunc[nv].SetPoint(np,fitDataEnMean,meanFit)
                calibFunc[nv].SetPointError(np,fitDataEnSigma,meanFitError)
                
                #show results
                ipad=ipad+1
                c.cd(ipad)
                pframe=ws.var(vName).frame(ROOT.RooFit.Range(postfix))
                redData.plotOn(pframe)
                ws.pdf('g_%s'%postfix).plotOn(pframe,ROOT.RooFit.Range(postfix))
                pframe.Draw()
                pframe.GetXaxis().SetTitle(varTitle )
                pframe.GetYaxis().SetTitle('Events')
                pframe.GetYaxis().SetTitleOffset(1.2)
                pframe.GetYaxis().SetRangeUser(0.01,1.8*pframe.GetMaximum())
                pframe.GetXaxis().SetNdivisions(5)
                pt=MyPaveText('[Energy=%d GeV, %3.1f<#eta<%3.1f]\\<E_{raw}>=%3.2f RMS=%3.2f\\#mu=%3.2f #sigma=%3.2f'%
                              (fitDataEnMean,etaMin,etaMax,fitDataMean,fitDataSigma,meanFit,sigmaFit),
                              0.18,0.9,0.5,0.6)
                pt.SetTextFont(42)
                pt.SetTextSize(0.05)
                if ipad==1:
                    simpt=MyPaveText('#bf{CMS} #it{simulation}')
                    simpt.SetTextSize(0.06)
                    simpt.SetTextFont(42)
                               

            #save the canvas with the fitted values
            for ext in ['png','pdf','C'] : c.SaveAs(outDir+'/fits_%s_eta%3.1f.%s'%(var.GetName(),avgEta,ext))

            #add the calibrated variable bias for the resolution fit
            calibFunc[nv].Fit(calibModel,'MER+')
            calibFunc[nv].GetFunction(calibModel.GetName()).SetLineColor(calibFunc[nv].GetLineColor())
            calib_offset=calibFunc[nv].GetFunction(calibModel.GetName()).GetParameter(1)
            calib_slope=calibFunc[nv].GetFunction(calibModel.GetName()).GetParameter(0)
            #calibFname='%sCalibBias_%d'%(var.GetName(),nv)
            #varValName=var.GetName().replace('Func','')
            #ws.factory("RooFormulaVar::%s('(@0-%f)/%f-@1',{%s,en})"%(calibFname,calib_offset,calib_slope,varValName))
            calibFname='%sCalib_%d'%(var.GetName(),nv)
            varValName=var.GetName().replace('Func','')
            ws.factory("RooFormulaVar::%s('(@0-%f)/%f',{%s})"%(calibFname,calib_offset,calib_slope,varValName))

            
            #use the values
            calibVname=calibFname.replace('Func','')
            calibEnEstimatorsSet=ROOT.RooArgSet(ws.var('eta'),ws.var('en'),ws.var('phi'),ws.factory('%s[0,-9999999,9999999]'%calibVname))
            getattr(ws,'import')( ROOT.RooDataSet('fitdata'+calibVname,'fitdata'+calibVname,calibEnEstimatorsSet) )
            for ientry in xrange(0,ws.data('fitdata').numEntries()):
                entryVars=ws.data('fitdata').get(ientry)
                for baseVar in ['eta','en','phi',varValName]: ws.var(baseVar).setVal( entryVars.find(baseVar).getVal() )
                newEntry=ROOT.RooArgSet( ws.var('eta'), ws.var('en'), ws.var('phi') )
                ws.var(calibVname).setVal( ws.function(calibFname).getVal() )
                newEntry.add( ws.var(calibVname) )
                ws.data('fitdata'+calibVname).add( newEntry )
                if ROOT.TMath.Abs(ws.function(calibFname).getVal()/ws.var('en').getVal())>0.2:
                    outliersPhi.Fill(ws.var('phi').getVal())
                
            
            #repeat for resolution fit but with calibrated energy
            c.Clear()
            c.Divide(len(enRanges)/3+2,2)
            ipad=0
            for enMin,enMax in enRanges:
                postfix='fit%d%d%d%d'%(fitCtr,nv,iEtaRange,ipad)
                ipad=ipad+1

                #cut data
                redData=ws.data('fitdata'+calibVname).reduce('en>=%f && en<=%f && eta>=%f && eta<=%f'%(enMin,enMax,etaMin,etaMax))
                if redData.numEntries()<10 : 
                    ipad=ipad-1
                    continue
                
                #calibrate and fit the residual bias
                calibFitDataMean,   calibFitDataSigma   = redData.mean(ws.var(calibVname)), redData.sigma(ws.var(calibVname))
                fitDataEnMean, fitDataEnSigma = redData.mean(ws.var('en')),  redData.sigma(ws.var('en'))
                
                xaxisMin,xaxisMax=calibFitDataMean-5*calibFitDataSigma,calibFitDataMean+5*calibFitDataSigma
                xfitMin, xfitMax=calibFitDataMean-nSigmasToFit*calibFitDataSigma,calibFitDataMean+nSigmasToFit*calibFitDataSigma
                #if calibFitDataSigma/calibFitDataMean > 0.2:
                #    xaxisMin, xaxisMax = calibFitDataMean*0.25, calibFitDataMean*4
                #    xfitMin,  xfitMax = calibFitDataMean*0.9, calibFitDataMean*1.1
                ws.var(calibVname).setRange(postfix+calibVname,xaxisMin,xaxisMax)
                ws.var(calibVname).setRange('fit'+postfix+calibVname,xfitMin,xfitMax)
                #ws.factory('Gaussian::gcalib_%s_%s(%s,calibmean_%s_%s[%f,%f,%f],calibsigma_%s_%s[%f,%f,%f])'%
                #           (postfix, calibVname,  calibVname,
                #            postfix, calibVname,  0,  -5*calibFitDataSigma, +5*calibFitDataSigma,
                #            postfix, calibVname,  calibFitDataSigma, calibFitDataSigma*0.001,              calibFitDataSigma*2)
                #    )
                ws.factory('RooCBShape::gcalib_%s_%s(%s,calibmean_%s_%s[%f,%f,%f],calibsigma_%s_%s[%f,%f,%f],calibalpha_%s_%s[0.001.,0.0001,20.0],calibn_%s_%s[1])'%
                           (postfix,calibVname, calibVname,
                            postfix,calibVname,calibFitDataMean,calibFitDataMean-calibFitDataSigma,calibFitDataMean+calibFitDataSigma,
                            postfix,calibVname,calibFitDataSigma,calibFitDataSigma*0.001,calibFitDataSigma*2,
                            postfix,calibVname,postfix,calibVname)
                           )

                fres=ws.pdf('gcalib_%s_%s'%(postfix,calibVname)).fitTo( redData, ROOT.RooFit.Range('fit'+postfix+calibVname), ROOT.RooFit.Save(True) )
                meanFit       = ws.var('calibmean_%s_%s'%(postfix,calibVname)).getVal()
                sigmaFit      = ws.var('calibsigma_%s_%s'%(postfix,calibVname)).getVal()
                sigmaFitError = ws.var('calibsigma_%s_%s'%(postfix,calibVname)).getError()

                #save results
                np=resFunc[nv].GetN()
                resFunc[nv].SetPoint(np,fitDataEnMean,sigmaFit/fitDataEnMean)
                resFunc[nv].SetPointError(np,fitDataEnSigma,sigmaFitError/fitDataEnMean)

                #show results
                c.cd(ipad)
                pframe=ws.var(calibVname).frame(ROOT.RooFit.Range(postfix+calibVname))
                redData.plotOn(pframe)
                ws.pdf('gcalib_%s_%s'%(postfix,calibVname)).plotOn(pframe,ROOT.RooFit.Range(postfix+calibVname))
                pframe.Draw()
                pframe.GetXaxis().SetTitle( 'Calibrated ' + varTitle + ' [GeV]' )
                pframe.GetYaxis().SetTitle( 'Events' )
                pframe.GetYaxis().SetTitleOffset(1.2)
                pframe.GetYaxis().SetRangeUser(0.01,1.8*pframe.GetMaximum())
                pframe.GetXaxis().SetNdivisions(5)
                pt=MyPaveText('[Energy=%d GeV, |#eta|=%3.1f]\\<E>=%3.2f RMS=%3.2f\\#mu=%3.2f #sigma=%3.2f'%
                              (fitDataEnMean,avgEta,calibFitDataMean,calibFitDataSigma,meanFit,sigmaFit),
                    0.18,0.9,0.5,0.6)
                pt.SetTextFont(42)
                pt.SetTextSize(0.05)
                if ipad==1:
                    simpt=MyPaveText('#bf{CMS} #it{simulation}')
                    simpt.SetTextSize(0.06)
                    simpt.SetTextFont(42)
                
            c.Modified()
            c.Update()
            for ext in ['png','pdf','C'] : c.SaveAs(outDir+'/biasfits_%s_eta%3.1f.%s'%(vName,avgEta,ext))

        #show the outliers phi distribution
        c.Clear()
        outliersPhi.Draw()
        simpt=MyPaveText('#bf{CMS} #it{simulation}')
        simpt.SetTextSize(0.06)
        simpt.SetTextFont(42)
        for ext in ['png','pdf','C'] : c.SaveAs(outDir+'/%s_outliers.%s'%(vName,ext))

    showEMcalibrationResults(methodNames=methodNames,calibFunc=calibFunc,resFunc=resFunc,outDir=outDir,isHadronic=isHadronic)

"""
steer 
"""
def main():
    
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i',      '--in' ,      dest='input',        help='Input file',                              default=None)
    parser.add_option('--isHad' ,              dest='isHad',        help='Flag if it\'s an hadronic calibration',   default=False, action="store_true")
    parser.add_option('--offset' ,             dest='deriveOffset', help='Flag if calibration offset is to be fit', default=False, action="store_true")
    parser.add_option('-w',      '--ws' ,      dest='ws',           help='ROOT file with previous workspace',       default=None)
    parser.add_option('-v',      '--var' ,     dest='treeVarName',  help='Variable to use as energy estimotor',     default='edep_thr0')
    (opt, args) = parser.parse_args()

     #check inputs                                                                                                                                                                                                  
    if opt.input is None :
        parser.print_help()
        sys.exit(1)

    #basic ROOT customization
    customROOTstyle()
    #ROOT.gROOT.SetBatch(False)
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    
    runResolutionStudy(url=opt.input,treeVarName=opt.treeVarName,isHadronic=opt.isHad,deriveOffset=opt.deriveOffset,wsOutUrl=opt.ws)
    
if __name__ == "__main__":
    sys.exit(main())
