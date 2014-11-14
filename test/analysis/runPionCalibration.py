#!/usr/bin/env python

import ROOT
import io, os, sys
import optparse
import commands
from UserCode.HGCanalysis.PlotUtils import *
from UserCode.HGCanalysis.HGCTree2Workspace import *


"""
Adapts the workspace for pion calibration
"""
def adaptWorkspaceForPionCalibration(opt):
    url               = opt.input
    wsUrl             = opt.wsUrl
    treeVarName       = opt.treeVarName
    vetoTrackInt      = opt.vetoTrackInt
    vetoHEBLeaks      = opt.vetoHEBLeaks
    emCalibUrls       = {}
    if opt.emCalibUrl :
        for url in opt.emCalibUrl.split(','):
            key_url=url.split(':')
            emCalibUrls[ key_url[0] ] = key_url[1]
    
    #init phase space regions of interest
    etaRanges = [[1.6,1.75],[1.75,2.0],[2.0,2.25],[2.25,2.5],[2.5,2.75],[2.75,2.9]]
    enRanges  = [[9,11],[19,21],[39,41],[49,51],[74,75],[99,101],[249,251]]
    
    #init weights and integration ranges
    subDets            = ['EE','HEF','HEB']
    integRanges        = [[1,1],[2,11],[12,21],[22,30],[31,31],[32,42],[43,54]]
    subDetRanges       = ['EE', 'EE',  'EE',   'EE',   'HEF',  'HEF',  'HEB']
    weights            = {}
    weights["lambda"]    = [0.01, 0.036, 0.043,  0.056,  0.338,  0.273,  0.476]
    weights["lambda_em"] = [0.01, 0.036, 0.043,  0.056,  0.338,  0.273,  0.476]
    weightTitles={}
    weightTitles["lambda"]  = "#lambda-based weights"
    weightTitles["lambda_em"]  = "#lambda-based + e.m. scale weights"

    #prepare workspace (if needed) and output
    outDir="./"
    if wsUrl is None :
        outDir=url.replace('.root','')
        os.system('mkdir -p '+outDir)
        wsUrl=prepareWorkspace(url=url,integRanges=integRanges,vetoTrackInt=vetoTrackInt,vetoHEBLeaks=vetoHEBLeaks,treeVarName=treeVarName)
    else:
        outDir=os.path.dirname(wsUrl)

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
    for key in emCalibUrls:
        if key in subDets:
            calibF=ROOT.TFile.Open(emCalibUrls[key])
            if calibF:
                print 'Reading out e.m. calibration for %s from %s'%(key,emCalibUrls[key])
                emCalibMap[key]=calibF.Get('lambda_calib').Clone('%s_lambda_calib'%key)
                calibF.Close()

    #global calibration (tbd)
    calibPostFix=''

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

            ws.factory('%sEn_%s[0,0,999999999]'%(wType,subDet))
            ws.var('%sEn_%s'%(wType,subDet)).SetTitle( '%s, %s'%(subDet,weightTitles[wType] ) )
            theFunc=ws.factory("RooFormulaVar::%sEn_%s_func('%s',%s)"%(wType,subDet,funcFormula[wType],funcArgs))
            ws.data('data').addColumn( theFunc ) 

    #now create the sum of the sub-detectors
    for wType in weights:

        funcArgs='{'
        funcFormula='('
        varCtr=0
        for subDet in subDets:

            if varCtr>0 :
                funcArgs += ','
                funcFormula+='+'
            funcArgs+='%sEn_%s_func'%(wType,subDet)
            calib_slope, calib_offset = 1.0, 0.0
            funcFormula += '(@%d-%f)/%f'%(varCtr,calib_offset,calib_slope)
            varCtr+=1
        funcArgs += '}'
        funcFormula+=')'
        ws.factory('%sEn[0,0,999999999]'%(wType))
        ws.var('%sEn'%wType).SetTitle( '%s'%weightTitles[wType])
        theFunc=ws.factory("RooFormulaVar::%sEn_func('%s',%s)"%(wType,funcFormula,funcArgs))
        
    ws.Print('v')
    raw_input()

    #finally create dataset for calibration
    uncalibDataVars=ROOT.RooArgSet(ws.var('en'), ws.var('eta'), ws.var('phi'))
    for wType in weights: 
        uncalibDataVars.add( ws.var('%sEn_%s'%(wType,calibPostFix) ) )
        for subDet in subDets:
            uncalibDataVars.add( ws.var('%sEn_%s_%s'%(wType,subDet,calibPostFix) ) )
    getattr(ws,'import')( ROOT.RooDataSet('data_%s'%calibPostFix,'data_%s'%calibPostFix,uncalibDataVars) )

    #fill the dataset
    for ientry in xrange(0,ws.data('data').numEntries()):
        entryVars=ws.data('data').get(ientry)
        newEntry=ROOT.RooArgSet()
        for baseVar in ['en','eta','phi']: 
            ws.var(baseVar).setVal( entryVars.find(baseVar).getVal() )
            newEntry.add( ws.var(baseVar) )

            #apply material correction before HGC
            if matParamBeforeHGC :
                genEta = ROOT.TMath.Abs(ws.var('eta').getVal())
                edep0_corr  = ws.var('edep0').getVal()/ROOT.TMath.TanH(genEta) 
                edep0_corr  *= matParamBeforeHGCMap.Eval(genEta)*int(vetoTrackInt) 
                ws.var('edep0').setVal('edep0')

            for wType in weights :
                #per sub-detector
                for subDet in subDets:
                    ienVal=entryVars.find('%sEn_%s_func'%(wType,subDet)).getVal()
                    ws.var('%sEn_%s_%s'%(wType,subDet,calibPostFix)).setVal(ienVal)
                    newEntry.add(ws.var('%sEn_%s_%s'%(wType,subDet,calibPostFix)))

                #summed over all sub-detectors
                ienVal=entryVars.find('%sEn_func'%wType).getVal()
                ws.var('%sEn_%s'%(wType,calibPostFix)).setVal(ienVal)
                newEntry.add(ws.var('%sEn_%s'%(wType,calibPostFix)))

            #all filled, add new row
            ws.data('data_%s'%calibPostFix).add( newEntry )

    newWsUrl=wsUrl.replace('.root','_%s.root'%calibPostFix)
    ws.writeToFile(newWsUrl,true)
    print 'Workspace adapted for pion calibration, stored @ %s'%newWsUrl
    return newWsUrl

"""
Steer the calibration study
"""
def runCalibrationStudy(opt):

    #get the workspace
    wsUrl=adaptWorkspaceForPionCalibration(opt)
    wsOutF=ROOT.TFile.Open(wsUrl)
    ws=wsOutF.Get('w')
    wsOutF.Close()

    #verify sub-detector correlations

    #calibrate the energy estimators (split up in different energies and pseudo-rapidity ranges)
#    nSigmasToFit=2 #2.8
#    calibGr={}
#    resGr={}
#    for wType in weights:
#        calibGr[wType]=ROOT.TMultiGraph()
#        calibGr[wType].SetName('calib_%s'%wType)
#        calibGr[wType].SetTitle( weightTitles[wType] )
#        resGr[wType]=calibGr[wType].Clone('res_%s'%wType)
#
#    for iEtaRange in xrange(0,len(etaRanges)):
#        
#        genEta_min=etaRanges[iEtaRange][0]
#        genEta_max=etaRanges[iEtaRange][1]
#        genEta_mean=0.5*(genEta_max+genEta_min)
#        
#        #keep track of eta slices separately
#        etaSliceCalibGr={}
#        etaSliceResGr={}
#        iwtype=0
#        for wType in weights:
#            iwtype=iwtype+1
#            etaSliceCalibGr[wType]=ROOT.TGraphErrors()
#            etaSliceCalibGr[wType].SetName('calib_%s_%s'%(iEtaRange,wType))
#            etaSliceCalibGr[wType].SetTitle('%3.1f<|#eta|<%3.1f'%(genEta_min,genEta_max))
#            etaSliceCalibGr[wType].SetLineColor(iwtype)
#            etaSliceCalibGr[wType].SetMarkerColor(iwtype)
#            etaSliceCalibGr[wType].SetMarkerStyle(iEtaRange+20)
#            etaSliceResGr[wType]=etaSliceCalibGr[wType].Clone('res_%s_%s'%(iEtaRange,wType))
#
#        for iEnRange in xrange(0,len(enRanges)):
#            genEn_min = enRanges[iEnRange][0]
#            genEn_max = enRanges[iEnRange][1]
#            redData=ws.data('data_%s'%calibPostFix).reduce('en>=%f && en<=%f && eta>=%f && eta<=%f'%(genEn_min,genEn_max,genEta_min,genEta_max))
#            if redData.numEntries()<10 : continue
#            genEn_mean, genEn_sigma = redData.mean(ws.var('en')),  redData.sigma(ws.var('en'))
#
#            for wType in weights :
#
#                #prepare the fit to this slice
#                vName = '%sEn_%s'%(wType,calibPostFix)
#                v_mean, v_sigma    = redData.mean(ws.var(vName)), redData.sigma(ws.var(vName))
#                v_min, v_max       = v_mean-5*v_sigma, v_mean+5*v_sigma
#                v_fitMin, v_fitMax = v_mean-nSigmasToFit*v_sigma, v_mean+nSigmasToFit*v_sigma
#
#                #define PDF
#                fitName          = 'range%d%d_%s'%(iEtaRange,iEnRange,vName)            
#                ws.var(vName).setRange(fitName,v_min,v_max)
#                ws.var(vName).setRange('fit_%s'%fitName,v_fitMin, v_fitMax)
#                ws.factory('RooCBShape::resol_%s(%s,mean_%s[%f,%f,%f],sigma_%s[%f,%f,%f],alpha_%s[0.001.,0.0001,20.0],n_%s[1])'%
#                           (fitName,vName,
#                            fitName,v_mean,v_min,v_max,
#                            fitName,v_sigma,v_sigma*0.001, v_sigma*2,
#                            fitName,
#                            fitName)
#                           )
#                
#                #fit
#                fres = ws.pdf('resol_%s'%fitName).fitTo( redData, ROOT.RooFit.Range('fit_%s'%fitName), ROOT.RooFit.Save(True) )
#                meanFit, meanFit_error   = ws.var('mean_%s'%fitName).getVal(), ws.var('mean_%s'%fitName).getError()
#                sigmaFit, sigmaFit_error = ws.var('sigma_%s'%fitName).getVal(), ws.var('sigma_%s'%fitName).getError()
#
#                #save results
#                np=etaSliceCalibGr[wType].GetN()
#                etaSliceCalibGr[wType].SetPoint(np,genEn_mean,meanFit)
#                etaSliceCalibGr[wType].SetPointError(np,0,meanFit_error)
#                etaSliceResGr[wType].SetPoint(np,genEn_mean,sigmaFit/meanFit)
#                etaSliceResGr[wType].SetPointError(np,0,ROOT.TMath.Sqrt( ROOT.TMath.Power(meanFit*sigmaFit_error,2)+
#                                                                          ROOT.TMath.Power(meanFit_error*sigmaFit,2) ) /
#                                                    ROOT.TMath.Power(meanFit,2) )
#
#                #print genEn_min,genEn_max,genEn_mean,wType
#                #raw_input()
#                
#                #for debug purposes
#                showCalibrationFitResults( theVar=ws.var(vName),
#                                           theData=redData,
#                                           thePDF=ws.pdf('resol_%s'%fitName),
#                                           theLabel='#it{Energy=%d GeV, %3.1f<#eta<%3.1f}\\#mu=%3.2f#pm%3.2f\\#sigma=%3.2f#pm%3.2f'%(genEn_mean,genEta_min,genEta_max,meanFit,meanFit_error,sigmaFit,sigmaFit_error),
#                                           fitName=fitName,
#                                           outDir=outDir)
#
#        for wType in weights: 
#            calibGr[wType].Add(etaSliceCalibGr[wType],'p')
#            resGr[wType].Add(etaSliceResGr[wType],'p')
#
#    #derive calibration
#    calibModel=ROOT.TF1('calibmodel',"[0]*x+[1]",0,800)
#    calibModel.SetLineWidth(1)
#    calibF=ROOT.TFile.Open('%s/calib.root'%outDir,'RECREATE')
#    for wType in weights :
#        calibGr[wType].Fit(calibModel,'MER+')
#        calibGr[wType].GetFunction(calibModel.GetName()).SetLineColor(calibGr[wType].GetListOfGraphs().At(0).GetLineColor())
#        calibGr[wType].GetFunction(calibModel.GetName()).Clone('%s_calib'%wType).Write()
#    calibF.Close()
#
#    #show results
#    showCalibrationCurves(calibGr=calibGr,outDir=outDir,calibPostFix=calibPostFix)
#    showResolutionCurves(resGr=resGr,outDir=outDir,calibPostFix=calibPostFix)
#

"""
steer 
"""
def main():
    
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i',      '--in' ,      dest='input',        help='Input file',                                     default=None)
    parser.add_option('-w',      '--ws' ,      dest='wsUrl',        help='Workspace file',                                 default=None)
    parser.add_option('--emCalib' ,            dest='emCalibUrl',   help='em calibration files (e.g. EE:calib_ee.root,HEF:calib_hef.root', default=None)
    parser.add_option('--vetoTrackInt',        dest='vetoTrackInt', help='flag if tracker interactions should be removed', default=False, action="store_true")
    parser.add_option('--vetoHEBLeaks',        dest='vetoHEBLeaks', help='flag if HEB leaks are allowed',                  default=False, action='store_true')
    parser.add_option('-v',      '--var' ,     dest='treeVarName',  help='Variable to use as energy estimator',            default='edep_sim')
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
