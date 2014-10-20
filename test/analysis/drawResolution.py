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
def prepareWorkspace(url,integRanges,mipEn,outUrl):
    
    #prepare the workspace
    ws=ROOT.RooWorkspace("w")
    dsVars=ROOT.RooArgSet( ws.factory('eta[1.5,1.45,3.1]'), ws.factory('en[0,0,9999999999]') )
    for ireg in xrange(0,len(integRanges)): dsVars.add( ws.factory('edep%d[0,0,99999999]'%ireg) )
    getattr(ws,'import')( ROOT.RooDataSet('data','data',dsVars) )

    #optimization with linear regression
    optimVec    = numpy.zeros( len(integRanges)+1 )
    optimMatrix = numpy.zeros( (len(integRanges)+1, len(integRanges)+1 ) )

    #read all to a RooDataSet
    fin=ROOT.TFile.Open(url)
    HGC=fin.Get('analysis/HGC')
    for entry in xrange(0,HGC.GetEntriesFast()+1):
        HGC.GetEntry(entry)
        genEn=HGC.genEn
        genEta=ROOT.TMath.Abs(HGC.genEta)
        if genEta<1.4 or genEta>3.0 : continue
        ws.var('eta').setVal(genEta)
        ws.var('en').setVal(genEn)
        newEntry=ROOT.RooArgSet( ws.var('eta'), ws.var('en') )

        #get the relevant energy deposits and add new row
        edeps=[]
        for ireg in xrange(0,len(integRanges)):
            totalEnInIntegRegion=0
            for ilayer in xrange(integRanges[ireg][0],integRanges[ireg][1]+1):
                totalEnInIntegRegion=totalEnInIntegRegion+HGC.edeps[ilayer-1]
            edeps.append(totalEnInIntegRegion*1e6 / mipEn[ireg] )
            ws.var('edep%d'%ireg).setVal(totalEnInIntegRegion)
            newEntry.add(ws.var('edep%d'%(ireg)))
        ws.data('data').add( newEntry )

        #fill the optmization matrix and vector
        for ie in xrange(0,len(edeps)):
            optimVec[ie]=optimVec[ie]+edeps[ie]*genEn
            for je in xrange(0,len(edeps)):
                optimMatrix[ie][je]=optimMatrix[ie][je]+edeps[ie]*edeps[je]
            optimMatrix[len(edeps)][ie]=optimMatrix[len(edeps)][ie]+edeps[ie]
            optimMatrix[ie][len(edeps)]=optimMatrix[ie][len(edeps)]+edeps[ie]
        optimMatrix[len(edeps)][len(edeps)]=optimMatrix[len(edeps)][len(edeps)]+1
        optimVec[len(edeps)]=optimVec[len(edeps)]+genEn        
    fin.Close()

    #all done, write to file
    ws.writeToFile(outUrl,True)
    print 'Created the analysis RooDataSet with %d events, stored @ %s'%(ws.data('data').numEntries(),outUrl)
    
    #finalize weight optimization
    optimWeights=numpy.linalg.solve(optimMatrix,optimVec)
    return optimWeights

    

"""
runs the fits to the calibration and resolution
"""
def showEMcalibrationResults(calibFunc,resFunc,outDir,isHadronic=False):

    nEta=len(resFunc)/3

    resolSummary=[]

    #
    #resolution
    #
    resolModel=ROOT.TF1('resolmodel',"sqrt([0]*[0]/x+[1]*[1])",0,1000)
    #resolModel=ROOT.TF1('resolmodel',"sqrt([0]*x*x+[1])",0,1)
    #resolModel=ROOT.TF1('resolmodel',"sqrt([0]*x*x+[1]+[2]*x*x*x*x)",0,1)
    resolModel.SetParameter(0,0.2);
    resolModel.SetParLimits(0,0,2);
    resolModel.SetParameter(1,0);
    resolModel.SetParLimits(1,0,1.0);
    resolModel.SetLineWidth(1)

    c=ROOT.TCanvas('cresol','cresol',1200,500)
    c.SetTopMargin(0.05)
    c.SetRightMargin(0.5)
    c.SetLogx()
    leg=ROOT.TLegend(0.52,0.1,0.99,0.2+0.085*len(calibFunc)/3)
    leg.SetNColumns(3)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.025)
    for ir in xrange(0,len(resFunc)):
        if ir==0:
            resFunc[ir].Draw('ap')
            resFunc[ir].GetXaxis().SetTitle("E [GeV]")
            resFunc[ir].GetYaxis().SetTitle("Relative energy resolution") 
            resFunc[ir].GetXaxis().SetLabelSize(0.04)
            resFunc[ir].GetYaxis().SetLabelSize(0.04)
            resFunc[ir].GetXaxis().SetTitleSize(0.05)
            resFunc[ir].GetYaxis().SetTitleSize(0.05)
            resFunc[ir].GetYaxis().SetTitleOffset(1.3)
            resFunc[ir].GetXaxis().SetMoreLogLabels(True)
            maxY=resFunc[ir].GetMaximum()
            if isHadronic: resFunc[ir].GetYaxis().SetRangeUser(0,1.0)
            else: resFunc[ir].GetYaxis().SetRangeUser(0,0.1)
        else:
            resFunc[ir].Draw('p')

        if ir%nEta==0 :
            resFunc[ir].Fit(resolModel,'MER+')
            resFunc[ir].GetFunction(resolModel.GetName()).SetLineColor(resFunc[ir].GetLineColor())
        else :
            resFunc[ir].Fit(resolModel,'MER0+')
        
        sigmaStoch=resolModel.GetParameter(0)
        sigmaStochErr=resolModel.GetParError(0)
        sigmaConst=resolModel.GetParameter(1)
        sigmaConstErr=resolModel.GetParError(1)

        #sigmaStoch=ROOT.TMath.Sqrt(resolModel.GetParameter(0));
        #if sigmaStoch>0 : sigmaStochErr=resolModel.GetParError(0)/(2*ROOT.TMath.Sqrt(sigmaStoch))
        #sigmaConstErr=0
        #sigmaConst=ROOT.TMath.Sqrt(resolModel.GetParameter(1))
        #if sigmaConst>0 : sigmaConstErr=resolModel.GetParError(1)/(2*ROOT.TMath.Sqrt(sigmaConst))
        #sigmaNoise=ROOT.TMath.Sqrt(resolModel.GetParameter(2))
        #sigmaNoiseErr=resolModel.GetParError(2)/(2*ROOT.TMath.Sqrt(sigmaNoise))

        leg.AddEntry(resFunc[ir],
                     #"#splitline{[#bf{#it{%s}}]}{#frac{#sigma}{E} #propto #frac{%3.3f #pm %3.3f}{#sqrt{E}} #oplus %3.4f#pm%3.4f}"
                     #"[#bf{#it{%s}}] #frac{#sigma}{E} #propto #frac{%3.4f #pm %3.4f}{#sqrt{E}} #oplus %3.5f#pm%3.5f"
                     "#splitline{[#scale[0.8]{#bf{#it{%s}}}]}{#frac{#sigma}{E} #propto #frac{%3.4f}{#sqrt{E}} #oplus %3.5f}"
                     %(resFunc[ir].GetTitle(),sigmaStoch,sigmaConst),
                     "fp")
        
        resolSummary.append( [sigmaStoch,sigmaStochErr,sigmaConst,sigmaConstErr] )

    leg.SetTextSize(0.03)
    leg.Draw()
    MyPaveText('#bf{CMS} #it{simulation}')

    #
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
                     #"#splitline{[#bf{#it{%s}}]}{E #propto (%3.1f#pm%3.1f)#timesE_{rec} +%3.0f#pm%3.0f}"
                     "#splitline{[#scale[0.8]{#bf{#it{%s}}}]}{#hat{E} = %3.1f#timesE_{beam} + %3.0f}"
                     %(calibFunc[ir].GetTitle(),slope,offset),
                     "fp")

        resolSummary[ir].append(slope)
        resolSummary[ir].append(slopeErr)
        resolSummary[ir].append(offset)
        resolSummary[ir].append(offsetErr)
        
    leg.Draw()
    MyPaveText('#bf{CMS} #it{simulation}')

    c.Modified()
    c.Update()
    for ext in ['png','pdf','C'] : c.SaveAs('%s/calibration.%s'%(outDir,ext))

    return resolSummary


"""
"""
def runResolutionStudy(url,isHadronic=False,wsOutUrl=None):
    
    #init weights
    mipEn       = [55.1,        55.1,        55.1,        55.1,        85.0,     1498.4]   
    integRanges = [[1,1],       [2,11],      [12,21],     [22,30],     [31,42],  [43,54]]
    defWeights  = [0.209/0.494, 0.494/0.494, 0.797/0.494, 1.207/0.494, 0.,       0.,      0.]
    if isHadronic :
        #defWeights = [0.028,       0.065,        0.105,       0.160,      1.0,       1.667,   0.]        
        defWeights = [0.028,       0.065,        0.105,       0.160, 0.]
    optimWeights=defWeights
    
    #prepare output
    outDir=url.replace('.root','')
    os.system('mkdir -p '+outDir)

    #get workspace
    if wsOutUrl is None:
        wsOutUrl='%s/workspace.root'%outDir
        optimWeights=prepareWorkspace(url=url,integRanges=integRanges,mipEn=mipEn,outUrl=wsOutUrl)
    wsOutF=ROOT.TFile.Open(wsOutUrl)
    ws=wsOutF.Get('w')
    wsOutF.Close()

    #output weights to a file
    calibrationData={}
    calibrationData['IntegrationRanges'] = [ {'first':fLayer, 'last':lLayer} for fLayer,lLayer in integRanges ]
    calibrationData['MIPinKeV']          = [ item for item in mipEn ]
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
    vars=[
        [ws.data('data').addColumn(ws.factory("RooFormulaVar::rawEnFunc('"+rawEnFunc+"',"+funcArgs+")" )),                '#Sigma E_{i}'],
        [ws.data('data').addColumn(ws.factory("RooFormulaVar::weightEnFunc('"+weightEnFunc+"',"+funcArgs+")")),           '#Sigma w_{i}E_{i}' ],
        [ws.data('data').addColumn(ws.factory("RooFormulaVar::optimWeightEnFunc('"+optimWeightEnFunc+"',"+funcArgs+")")), '#Sigma w^{optim}_{i}E_{i}']
        ]
    enEstimatorsSet=ROOT.RooArgSet(ws.var('eta'),ws.var('en'))
    for v in vars:
        vName=v[0].GetName().replace('Func','')
        enEstimatorsSet.add( ws.factory('%s[0,0,999999999]'%vName) )
    getattr(ws,'import')( ROOT.RooDataSet('fitdata','fitdata',enEstimatorsSet) )

    #create the fit dataset (skip direct usage of RooFormulaVar in a fit) and store final value
    for ientry in xrange(0,ws.data('data').numEntries()):
        entryVars=ws.data('data').get(ientry)
        for baseVar in ['eta','en']: ws.var(baseVar).setVal( entryVars.find(baseVar).getVal() )
        newEntry=ROOT.RooArgSet( ws.var('eta'), ws.var('en') )
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
    nSigmasToFit=2.0

    etaRanges=[[1.5,2.8]]
    enRanges=[[99,101]]
    for var,varTitle in vars:
        varCtr+=1
        vName=var.GetName().replace('Func','')

        #run calibration and resolution fits in different rapidity ranges
        c.Clear()
        fitCtr=0
        for iEtaRange in xrange(0,len(etaRanges)):
            etaMin=etaRanges[iEtaRange][0]
            etaMax=etaRanges[iEtaRange][1]
            avgEta=0.5*(etaMax+etaMin)
            nv=len(calibFunc)
            calibFunc.append(ROOT.TGraphErrors())
            calibFunc[nv].SetMarkerStyle(20+iEtaRange) #nv+4*nv%2)
            calibFunc[nv].SetTitle('%s |#eta|=%3.1f'%(varTitle,avgEta))
            calibFunc[nv].SetFillStyle(0)
            calibFunc[nv].SetMarkerColor(varCtr+1)
            calibFunc[nv].SetLineColor(varCtr+1)
            calibFunc[nv].SetName('calib_%s_%d'%(vName,iEtaRange))
            resFunc.append(calibFunc[nv].Clone('resol_%s_%d'%(vName,iEtaRange)))

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
                ws.var(vName).setRange(postfix,fitDataMean-4*fitDataSigma,fitDataMean+4*fitDataSigma)
                ws.var(vName).setRange('fit'+postfix,fitDataMean*0.9,fitDataMean*1.1)
                ws.factory('Gaussian::g_%s(%s,mean_%s[%f,%f,%f],sigma_%s[%f,%f,%f])'%
                           (postfix,vName,
                            postfix,fitDataMean,fitDataMean-25*fitDataSigma,fitDataMean+25*fitDataSigma,
                            postfix,fitDataSigma,fitDataSigma*0.001,fitDataSigma*2)
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

            calibFunc[nv].Fit(calibModel,'MER+')
            calibFunc[nv].GetFunction(calibModel.GetName()).SetLineColor(calibFunc[nv].GetLineColor())
            calib_offset=calibFunc[nv].GetFunction(calibModel.GetName()).GetParameter(1)
            calib_slope=calibFunc[nv].GetFunction(calibModel.GetName()).GetParameter(0)

            #add the calibrated variable bias for the resolution fit
            calibFname='%sCalibBias_%d'%(var.GetName(),nv)
            varValName=var.GetName().replace('Func','')
            ws.factory("RooFormulaVar::%s('(@0-%f)/%f-@1',{%s,en})"%(calibFname,calib_offset,calib_slope,varValName))

            #use the values
            calibVname=calibFname.replace('Func','')
            calibEnEstimatorsSet=ROOT.RooArgSet(ws.var('eta'),ws.var('en'),ws.factory('%s[0,-9999999,9999999]'%calibVname))
            getattr(ws,'import')( ROOT.RooDataSet('fitdata'+calibVname,'fitdata'+calibVname,calibEnEstimatorsSet) )
            for ientry in xrange(0,ws.data('fitdata').numEntries()):
                entryVars=ws.data('fitdata').get(ientry)
                for baseVar in ['eta','en',varValName]: ws.var(baseVar).setVal( entryVars.find(baseVar).getVal() )
                newEntry=ROOT.RooArgSet( ws.var('eta'), ws.var('en') )
                ws.var(calibVname).setVal( ws.function(calibFname).getVal() )
                newEntry.add( ws.var(calibVname) )
                ws.data('fitdata'+calibVname).add( newEntry )
            
            #repeat for resolution fit but with calibrated energy
            c.Clear()
            c.Divide(len(enRanges)/3+2,2)
            ipad=0
            for ptRang in enRanges:
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
                
                if (not isHadronic) and calibFitDataSigma/calibFitDataMean>0.1 :
                    ws.var(calibVname).setRange(postfix+calibVname,-0.3*fitDataEnMean,0.3*fitDataEnMean)
                    ws.var(calibVname).setRange('fit'+postfix+calibVname,-0.1*fitDataEnMean,0.1*fitDataEnMean)
                else:
                    ws.var(calibVname).setRange(postfix+calibVname,calibFitDataMean-4*calibFitDataSigma,calibFitDataMean+4*calibFitDataSigma)
                    ws.var(calibVname).setRange('fit'+postfix+calibVname,calibFitDataMean-nSigmasToFit*calibFitDataSigma,calibFitDataMean+nSigmasToFit*calibFitDataSigma)
                ws.factory('Gaussian::gcalib_%s_%s(%s,calibmean_%s_%s[%f,%f,%f],calibsigma_%s_%s[%f,%f,%f])'%
                           (postfix, calibVname,  calibVname,
                            postfix, calibVname,  0,  -5*calibFitDataSigma, +5*calibFitDataSigma,
                            postfix, calibVname,  calibFitDataSigma, calibFitDataSigma*0.001,              calibFitDataSigma*2)
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
                pframe.GetXaxis().SetTitle( 'Calibrated ' + varTitle )
                pframe.GetYaxis().SetTitle( 'Events')
                pframe.GetYaxis().SetTitleOffset(1.2)
                pframe.GetYaxis().SetRangeUser(0.01,1.8*pframe.GetMaximum())
                pframe.GetXaxis().SetNdivisions(5)
                pt=MyPaveText('[Energy=%d GeV, |#eta|=%3.1f]\\<E>=%3.1f RMS=%3.1f GeV\\<bias>=%3.1f #sigma=%3.1f GeV'%
                              (fitDataEnMean,avgEta,calibFitDataMean+fitDataEnMean,calibFitDataSigma,meanFit,sigmaFit),
                    0.18,0.9,0.5,0.6)
                pt.SetTextFont(42)
                pt.SetTextSize(0.06)
                if ipad==1:
                    simpt=MyPaveText('#bf{CMS} #it{simulation}')
                    simpt.SetTextSize(0.06)
                    simpt.SetTextFont(42)
                
            c.Modified()
            c.Update()
            raw_input()
            for ext in ['png','pdf','C'] : c.SaveAs(outDir+'/efits_%s_eta%3.1f.%s'%(vName,avgEta,ext))

    return showEMcalibrationResults(calibFunc=calibFunc,resFunc=resFunc,outDir=outDir,isHadronic=isHadronic)

"""
steer 
"""
def main():
    
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i',      '--in' ,      dest='input',    help='Input file',                                     default='Single11_0.root')
    parser.add_option('--isHad' ,              dest='isHad',    help='Flag if it\'s an hadronic calibration',          default=False, action="store_true")
    parser.add_option('-o',      '--out' ,     dest='output',   help='Output file',                                    default='HitIntegrationAnalysis.root')
    parser.add_option('-w',      '--ws' ,      dest='ws',       help='ROOT file with previous workspace',              default=None)
    (opt, args) = parser.parse_args()

    #basic ROOT customization
    customROOTstyle()
    ROOT.gROOT.SetBatch(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    
    runResolutionStudy(url=opt.input,isHadronic=opt.isHad,wsOutUrl=opt.ws)
    
if __name__ == "__main__":
    sys.exit(main())
