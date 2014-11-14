#!/usr/bin/env python

import ROOT
import io,os,sys
import optparse
import commands
from UserCode.HGCanalysis.PlotUtils import *
from UserCode.HGCanalysis.HGCTree2Workspace import *


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
    #etaRanges = [[1.6,1.75],[1.75,2.0],[2.0,2.25],[2.25,2.5],[2.5,2.75],[2.75,2.9]]
    #enRanges  = [[9,11],[19,21],[39,41],[49,51],[74,75],[99,101],[249,251]] 
    
    etaRanges = [[1.6,2.9]]
    enRanges  = [[9,11],[29,31],[49,51]]
    
    #init weights and integration ranges
    integRanges       = [[1,1],[2,11],[12,21],[22,30],[31,31],[32,42],[43,54]]
    weights={}
    #weights["Trivial"] = [1.0,  1.0,   1.0,    1.0,    1.0,    1.0,    1.0]
    #weights["x0"]      = [0.08, 0.620, 0.809,  1.239,  3.580,  3.103,  5.228]
    weights["lambda"]  = [0.01, 0.036, 0.043,  0.056,  0.338,  0.273,  0.476]
    #weights["had_linreg"] = [0.060478641299467166,  0.010821968277230637, 0.0096620531479292698, 0.017915306115493686, 0.045039116195834318, 0.05297744952272318, 0.11989320894062623]
    #weights["em_linreg"]  = [ 0.041309493816764332,  0.008462210377895853, 0.009629802090212889, 0.016843239031382701, 0.1264599727885049, 0.15573504738067073, 0.2407144864311227]
    weightTitles={}
    weightTitles["Trivial"] = "Trivial weights"
    weightTitles["x0"]      = "X_{0}-based weights"
    weightTitles["lambda"]  = "#lambda-based weights"
    #weightTitles["had_linreg"] = "had lin. regression weights"
    #weightTitles["em_linreg"]  = "e.m. lin. regression weights"

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
    calibF=ROOT.TFile.Open('%s/calib_%s.root'%(outDir,calibPostFix),'RECREATE')
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
