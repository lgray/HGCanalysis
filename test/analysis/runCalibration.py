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
def prepareWorkspace(url,integRanges,vetoTrackInt,treeVarName):

    #optimization with linear regression
    optimVec = numpy.zeros( len(integRanges) )
    optimMatrix = numpy.zeros( (len(integRanges), len(integRanges) ) )
 
    #prepare the workspace
    outUrl=url.replace('.root','')
    os.system('mkdir -p %s'%outUrl)
    ws=ROOT.RooWorkspace("w")
    dsVars=ROOT.RooArgSet( ws.factory('eta[1.5,1.45,3.1]'), ws.factory('en[0,0,9999999999]'), ws.factory('phi[0,-3.2,3.2]') )
    for ireg in xrange(0,len(integRanges)): dsVars.add( ws.factory('edep%d[0,0,99999999]'%ireg) )
    getattr(ws,'import')( ROOT.RooDataSet('data','data',dsVars) )

    #read all to a RooDataSet
    fin=ROOT.TFile.Open(url)
    HGC=fin.Get('analysis/HGC')
    for entry in xrange(0,HGC.GetEntriesFast()+1):
        HGC.GetEntry(entry)

        #veto interactions in the tracker
        if vetoTrackInt and HGC.hasInteractionBeforeHGC : continue

        #require generated particle to be in the endcap
        genEn=HGC.genEn
        genEta=ROOT.TMath.Abs(HGC.genEta)
        genPhi=HGC.genPhi
        if genEta<1.4 or genEta>3.0 : continue
        ws.var('en').setVal(genEn)
        ws.var('eta').setVal(genEta)
        ws.var('phi').setVal(genPhi)
        newEntry=ROOT.RooArgSet(ws.var('en'), ws.var('eta'),  ws.var('phi') )

        #get the relevant energy deposits and add new row
        for ireg in xrange(0,len(integRanges)):
            totalEnInIntegRegion=0
            for ilayer in xrange(integRanges[ireg][0],integRanges[ireg][1]+1):
                totalEnInIntegRegion=totalEnInIntegRegion+(getattr(HGC,treeVarName))[ilayer-1]
            geomCorrection=ROOT.TMath.TanH(genEta)
            ws.var('edep%d'%ireg).setVal(totalEnInIntegRegion*geomCorrection)
            newEntry.add(ws.var('edep%d'%(ireg)))
        ws.data('data').add( newEntry )

        #for optimization
        for ie in xrange(0,len(integRanges)):
            optimVec[ie]=optimVec[ie]+ws.var('edep%d'%ie).getVal()/genEn
            for je in xrange(0,len(integRanges)):
                optimMatrix[ie][je]=optimMatrix[ie][je]+ws.var('edep%d'%ie).getVal()*ws.var('edep%d'%je).getVal()/(genEn*genEn)


    fin.Close()

    #finalize optimization
    try:
        optimWeights=numpy.linalg.solve(optimMatrix,optimVec)
        optimData={}
        optimData['IntegrationRanges'] = [ {'first':fLayer, 'last':lLayer} for fLayer,lLayer in integRanges ]
        optimData['OptimWeights'] = [ item for item in optimWeights ]
        with io.open('%s/optim_weights.dat'%outUrl, 'w', encoding='utf-8') as f: f.write(unicode(json.dumps(optimData, sort_keys = True, ensure_ascii=False, indent=4)))
    except:
        print 'Failed to optimize - singular matrix?'

    #all done, write to file
    wsFileUrl='%s/workspace.root'%outUrl
    ws.writeToFile(wsFileUrl,True)
    print 'Created the analysis RooDataSet with %d events, stored @ %s'%(ws.data('data').numEntries(),wsFileUrl)
    return wsFileUrl


"""
Helper function to show and save the results of the fit to a slice
"""
def showCalibrationFitResults(theVar,theData,thePDF,theLabel,fitName,outDir) :
    canvas=ROOT.TCanvas('c','c',500,500)
    
    pframe=theVar.frame(ROOT.RooFit.Range(fitName))
    theData.plotOn(pframe)
    thePDF.plotOn(pframe,ROOT.RooFit.Range(fitName))
    pframe.Draw()
    pframe.GetXaxis().SetTitle(theVar.GetTitle())
    pframe.GetYaxis().SetTitle('Events')
    pframe.GetYaxis().SetTitleOffset(1.2)
    pframe.GetYaxis().SetRangeUser(0.01,1.8*pframe.GetMaximum())
    pframe.GetXaxis().SetNdivisions(5)
    MyPaveText(theLabel,0.15,0.95,0.5,0.7).SetTextSize(0.035)
    MyPaveText('#bf{CMS} #it{simulation}')
    canvas.SaveAs('%s/%s.png'%(outDir,fitName))

"""
shows a set of calibration curves
"""
def showCalibrationCurves(calibGr,outDir,calibPostFix) :
    canvas=ROOT.TCanvas('c','c',500,500)
    canvas.cd()
    up=ROOT.TPad('up','up',0.0,0.4,1.0,1.0)
    up.SetBottomMargin(0.01)
    up.SetTopMargin(0.08)
    up.Draw()
    up.cd()
    leg=ROOT.TLegend(0.8,0.6,0.9,0.9)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.045)
    resCalibGr={}
    igr=0
    for wType in calibGr:
        if igr==0:
            calibGr[wType].Draw('a')
            calibGr[wType].GetXaxis().SetTitleSize(0)
            calibGr[wType].GetXaxis().SetLabelSize(0)
            calibGr[wType].GetYaxis().SetTitle('Reconstructed energy')
            calibGr[wType].GetYaxis().SetTitleOffset(0.9)
            calibGr[wType].GetYaxis().SetTitleSize(0.07)
            calibGr[wType].GetYaxis().SetLabelSize(0.05)
            calibGr[wType].GetYaxis().SetRangeUser(1,calibGr[wType].GetYaxis().GetXmax()*2)
            for gr in calibGr[wType].GetListOfGraphs():
                leg.AddEntry(gr,gr.GetTitle(),"p")
        else:
            calibGr[wType].Draw()
        igr+=1

        lcol=calibGr[wType].GetListOfGraphs().At(0).GetLineColor()
        ffunc=calibGr[wType].GetListOfFunctions().At(0)
        ffunc.SetLineColor(lcol)
        ffunc.SetLineStyle(9)
        calib_offset=ffunc.GetParameter(1)
        calib_slope=ffunc.GetParameter(0)
        MyPaveText('#it{%s} : %3.4f E_{rec} + %3.4f'%(calibGr[wType].GetTitle(),1./calib_slope,-calib_offset/calib_slope),
                   0.15,0.95-igr*0.05,0.4,0.92-igr*0.05).SetTextColor(lcol)

        #compute the residuals
        resCalibGr[wType]=ROOT.TMultiGraph()
        for gr in calibGr[wType].GetListOfGraphs() :
            newGr=gr.Clone('%s_res'%gr.GetName())
            newGr.Set(0)
            xval, yval = ROOT.Double(0), ROOT.Double(0)
            for ip in xrange(0,gr.GetN()):
                gr.GetPoint(ip,xval,yval)
                xval_error=gr.GetErrorX(ip)
                yval_error=gr.GetErrorY(ip)
                xrec=(yval-calib_offset)/calib_slope
                xrec_error=yval_error/calib_slope
                newGr.SetPoint(ip,xval,100*(xrec/xval-1))
                newGr.SetPointError(ip,xval_error,100*xrec_error/xval)
            resCalibGr[wType].Add(newGr,'p')
    leg.Draw()
    MyPaveText('#bf{CMS} #it{simulation}').SetTextSize(0.06)

    canvas.cd()
    dp=ROOT.TPad('dp','dp',0.0,0.0,1.0,0.4)
    dp.SetTopMargin(0.01)
    dp.SetBottomMargin(0.2)
    dp.Draw()
    dp.cd()
    igr=0
    for wType in resCalibGr:
        if igr==0:
            resCalibGr[wType].Draw('a')
            resCalibGr[wType].GetYaxis().SetRangeUser(-6.5,6.5)
            resCalibGr[wType].GetYaxis().SetNdivisions(5)
            resCalibGr[wType].GetXaxis().SetTitle('Generated energy [GeV]')
            resCalibGr[wType].GetYaxis().SetTitle('<E_{rec}-E_{gen}>/E_{gen} [%]')
            resCalibGr[wType].GetYaxis().SetTitleOffset(0.7)
            resCalibGr[wType].GetYaxis().SetTitleSize(0.09)
            resCalibGr[wType].GetYaxis().SetLabelSize(0.07)
            resCalibGr[wType].GetXaxis().SetTitleSize(0.09)
            resCalibGr[wType].GetXaxis().SetLabelSize(0.08)
        else:
            resCalibGr[wType].Draw()
        igr+=1

    canvas.cd()
    canvas.Modified()
    canvas.Update()
    canvas.SaveAs('%s/calib%s.png'%(outDir,calibPostFix))
    up.Delete()
    dp.Delete()

"""
shows a set of resolution curves
"""
def showResolutionCurves(resGr,outDir,calibPostFix) :

    resolModel=ROOT.TF1('resolmodel',"sqrt([0]*[0]/x+[1]*[1])",0,1000)
    resolModel.SetParameter(0,0.2);
    resolModel.SetParLimits(0,0,2);
    resolModel.SetParameter(1,0);
    resolModel.SetParLimits(1,0,1.0);

    canvas=ROOT.TCanvas('c','c',500,500)
    canvas.cd()
    leg=ROOT.TLegend(0.75,0.6,0.9,0.93)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.03)

    igr=0
    pt=[]
    for wType in resGr:
        if igr==0:
            resGr[wType].Draw('a')
            resGr[wType].GetXaxis().SetTitleSize(0)
            resGr[wType].GetXaxis().SetLabelSize(0)
            resGr[wType].GetXaxis().SetTitle('Generated energy [GeV]')
            resGr[wType].GetYaxis().SetTitle('#sigma_{E} / E')
            resGr[wType].GetYaxis().SetTitleOffset(1.3)
            resGr[wType].GetYaxis().SetTitleSize(0.07)
            resGr[wType].GetYaxis().SetLabelSize(0.05)
            resGr[wType].GetYaxis().SetRangeUser(1,resGr[wType].GetYaxis().GetXmax()*2)
            for gr in resGr[wType].GetListOfGraphs():
                leg.AddEntry(gr,gr.GetTitle(),"p")
        else:
            resGr[wType].Draw()
        igr+=1

        lcol=resGr[wType].GetListOfGraphs().At(0).GetLineColor()
        resGr[wType].Fit(resolModel,'MER+')
        ffunc=resGr[wType].GetListOfFunctions().At(0)
        ffunc.SetLineStyle(9)
        ffunc.SetLineColor(lcol)
        sigmaStoch    = ffunc.GetParameter(0)
        sigmaStochErr = ffunc.GetParError(0)
        sigmaConst    = ffunc.GetParameter(1)
        sigmaConstErr = ffunc.GetParError(1)
        pt.append( MyPaveText('#it{%s} :  #frac{%3.4f}{#sqrt{E}} #oplus %3.4f'%(resGr[wType].GetTitle(),sigmaStoch,sigmaConst),
                              0.2,0.93-igr*0.05,0.4,0.90-igr*0.05) )
        pt[igr-1].SetTextColor(lcol)
        pt[igr-1].SetTextSize(0.03)

    leg.Draw()
    MyPaveText('#bf{CMS} #it{simulation}').SetTextSize(0.04)

    canvas.Modified()
    canvas.Update()
    canvas.SaveAs('%s/resol%s.png'%(outDir,calibPostFix))




"""
Steer the calibration study
"""
def runCalibrationStudy(opt):

    url               = opt.input
    wsUrl             = opt.wsUrl
    calibUrl          = opt.calibUrl
    treeVarName       = opt.treeVarName
    vetoTrackInt      = opt.vetoTrackInt
    
    #init phase space regions of interest
    #etaRanges = [[1.5,2.0],[2.0,2.5],[2.5,2.9]]
    etaRanges = [[1.6,1.75],[1.75,2.0],[2.0,2.25],[2.25,2.5],[2.5,2.75],[2.75,2.9]]
    enRanges  = [[9,11],[19,21],[39,41],[49,51],[74,75],[99,101],[249,251]] #[4,6],[19,21],[29,31],[49,51],[74,76],[99,101],[149,151],[249,251],[499,501]]
    
    #init weights and integration ranges
    integRanges       = [[1,1],[2,11],[12,21],[22,30],[31,31],[32,42],[43,54]]
    weights={}
    #weights["Trivial"] = [1.0,  1.0,   1.0,    1.0,    1.0,    1.0,    1.0]
    weights["x0"]      = [0.08, 0.620, 0.809,  1.239,  3.580,  3.103,  5.228]
    weights["lambda"]  = [0.01, 0.036, 0.043,  0.056,  0.338,  0.273,  0.476]
    #weights["had_linreg"] = [0.060478641299467166,  0.010821968277230637, 0.0096620531479292698, 0.017915306115493686, 0.045039116195834318, 0.05297744952272318, 0.11989320894062623]
    #weights["em_linreg"]  = [ 0.041309493816764332,  0.008462210377895853, 0.009629802090212889, 0.016843239031382701, 0.1264599727885049, 0.15573504738067073, 0.2407144864311227]
    weightTitles={}
    weightTitles["Trivial"] = "Trivial weights"
    weightTitles["x0"]      = "X_{0}-based weights"
    weightTitles["lambda"]  = "#lambda-based weights"
    weightTitles["had_linreg"] = "had lin. regression weights"
    weightTitles["em_linreg"]  = "e.m. lin. regression weights"

    #prepare workspace (if needed) and output
    outDir="./"
    if wsUrl is None :
        outDir=url.replace('.root','')
        os.system('mkdir -p '+outDir)
        wsUrl=prepareWorkspace(url=url,integRanges=integRanges,vetoTrackInt=vetoTrackInt,treeVarName=treeVarName)
    else:
        outDir=os.path.dirname(wsUrl)

    #get the workspace from the file
    wsOutF=ROOT.TFile.Open(wsUrl)
    ws=wsOutF.Get('w')
    wsOutF.Close()

    #readout material overburden file
    matParamBeforeHGCMap={}
    matF=ROOT.TFile.Open("/afs/cern.ch/user/p/psilva/work/HGCal/Calib/CMSSW_6_2_0_SLHC20/src/UserCode/HGCanalysis/data/HGCMaterialOverburden.root")
    for wType in weights:
        matParamBeforeHGCMap[wType]=matF.Get("%sOverburden"%wType)
        matF.Close()

    #readout calibration
    calibPostFix='uncalib'
    calibMap={}
    if calibUrl:
        calibPostFix=os.path.basename(calibUrl)
        calibPostFix='_'+calibPostFix.replace('.root','')
        calibF=ROOT.TFile.Open(calibUrl)
        print 'Reading out calibrations from %s'%calibUrl
        for wType in weights:
            calibMap[wType]=calibF.Get('%s_calib'%wType).Clone()
        calibF.Close()

    #prepare energy estimators
    funcArgs=''
    funcFormula={}
    for wType in weights: funcFormula[wType]=''
    for ireg in xrange(0,len(integRanges)) :
        #list of function arguments
        if ireg==0:
            funcArgs += '{edep%d'%ireg
        else:
            funcArgs += ',edep%d'%ireg
        if ireg==len(integRanges)-1:
            funcArgs += '}'
        #functions
        for wType in weights:
            if ireg==0:
                funcFormula[wType] += '%f*edep%d'%(weights[wType][ireg],ireg)
            else:
                funcFormula[wType] += '+%f*edep%d'%(weights[wType][ireg],ireg)
            if ireg==len(integRanges)-1:
                funcFormula[wType] += ''
    for wType in weights: 
        calib_offset, calib_slope = 0.0, 1.0
        if wType in calibMap:
            calib_offset=calibMap[wType].GetParameter(1)
            calib_slope=calibMap[wType].GetParameter(0)
        ws.factory('%sEn_%s[0,0,999999999]'%(wType,calibPostFix))
        ws.var('%sEn_%s'%(wType,calibPostFix)).SetTitle( weightTitles[wType] )
        theFunc=ws.factory("RooFormulaVar::%sEn_func('(%s-%f)/%f',%s)"%(wType,funcFormula[wType],calib_offset,calib_slope,funcArgs))
        ws.data('data').addColumn( theFunc ) 

        
    #create dataset for calibration
    uncalibDataVars=ROOT.RooArgSet(ws.var('en'), ws.var('eta'), ws.var('phi'))
    for wType in weights: uncalibDataVars.add( ws.var('%sEn_%s'%(wType,calibPostFix) ) )
    getattr(ws,'import')( ROOT.RooDataSet('data_%s'%calibPostFix,'data_%s'%calibPostFix,uncalibDataVars) )
    for ientry in xrange(0,ws.data('data').numEntries()):
        entryVars=ws.data('data').get(ientry)
        newEntry=ROOT.RooArgSet()
        for baseVar in ['en','eta','phi']: 
            ws.var(baseVar).setVal( entryVars.find(baseVar).getVal() )
            newEntry.add( ws.var(baseVar) )

        for wType in weights :
            ienVal=entryVars.find('%sEn_func'%wType).getVal()

            #check if a material correction is available
            try:
                genEta = ROOT.TMath.Abs(ws.var('eta').getVal())
                edep0  = ws.var('edep0').getVal()/ROOT.TMath.TanH(genEta) #remove previously applied eta correction
                ienVal+= edep0*matParamBeforeHGCMap[wType].Eval(genEta)*int(vetoTrackInt) 
            except:
                pass
            ws.var('%sEn_%s'%(wType,calibPostFix)).setVal(ienVal)
            newEntry.add(ws.var('%sEn_%s'%(wType,calibPostFix)))
        ws.data('data_%s'%calibPostFix).add( newEntry )



    #calibrate the energy estimators (split up in different energies and pseudo-rapidity ranges)
    nSigmasToFit=2 #2.8
    calibGr={}
    resGr={}
    for wType in weights:
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
        for wType in weights:
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
            redData=ws.data('data_%s'%calibPostFix).reduce('en>=%f && en<=%f && eta>=%f && eta<=%f'%(genEn_min,genEn_max,genEta_min,genEta_max))
            if redData.numEntries()<10 : continue
            genEn_mean, genEn_sigma = redData.mean(ws.var('en')),  redData.sigma(ws.var('en'))

            for wType in weights :

                #prepare the fit to this slice
                vName = '%sEn_%s'%(wType,calibPostFix)
                v_mean, v_sigma    = redData.mean(ws.var(vName)), redData.sigma(ws.var(vName))
                v_min, v_max       = v_mean-5*v_sigma, v_mean+5*v_sigma
                v_fitMin, v_fitMax = v_mean-nSigmasToFit*v_sigma, v_mean+nSigmasToFit*v_sigma

                #define PDF
                fitName          = 'range%d%d_%s'%(iEtaRange,iEnRange,vName)            
                ws.var(vName).setRange(fitName,v_min,v_max)
                ws.var(vName).setRange('fit_%s'%fitName,v_fitMin, v_fitMax)
                ws.factory('RooCBShape::resol_%s(%s,mean_%s[%f,%f,%f],sigma_%s[%f,%f,%f],alpha_%s[0.001.,0.0001,20.0],n_%s[1])'%
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

                #print genEn_min,genEn_max,genEn_mean,wType
                #raw_input()
                
                #for debug purposes
                showCalibrationFitResults( theVar=ws.var(vName),
                                           theData=redData,
                                           thePDF=ws.pdf('resol_%s'%fitName),
                                           theLabel='#it{Energy=%d GeV, %3.1f<#eta<%3.1f}\\#mu=%3.2f#pm%3.2f\\#sigma=%3.2f#pm%3.2f'%(genEn_mean,genEta_min,genEta_max,meanFit,meanFit_error,sigmaFit,sigmaFit_error),
                                           fitName=fitName,
                                           outDir=outDir)

        for wType in weights: 
            calibGr[wType].Add(etaSliceCalibGr[wType],'p')
            resGr[wType].Add(etaSliceResGr[wType],'p')

    #derive calibration
    calibModel=ROOT.TF1('calibmodel',"[0]*x+[1]",0,800)
    calibModel.SetLineWidth(1)
    calibF=ROOT.TFile.Open('%s/calib.root'%outDir,'RECREATE')
    for wType in weights :
        calibGr[wType].Fit(calibModel,'MER+')
        calibGr[wType].GetFunction(calibModel.GetName()).SetLineColor(calibGr[wType].GetListOfGraphs().At(0).GetLineColor())
        calibGr[wType].GetFunction(calibModel.GetName()).Clone('%s_calib'%wType).Write()
    calibF.Close()

    #show results
    showCalibrationCurves(calibGr=calibGr,outDir=outDir,calibPostFix=calibPostFix)
    showResolutionCurves(resGr=resGr,outDir=outDir,calibPostFix=calibPostFix)

"""
steer 
"""
def main():
    
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i',      '--in' ,      dest='input',        help='Input file',                                     default=None)
    parser.add_option('-w',      '--ws' ,      dest='wsUrl',        help='Workspace file',                                 default=None)
    parser.add_option('-c',      '--calib' ,   dest='calibUrl',     help='Calibration file',                               default=None)
    parser.add_option('--vetoTrackInt',        dest='vetoTrackInt', help='flag if tracker interactions should be removed', default=False, action="store_true")
    parser.add_option('-v',      '--var' ,     dest='treeVarName',  help='Variable to use as energy estimotor',            default='edep_sim')
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
