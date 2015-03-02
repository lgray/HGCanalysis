#!/usr/bin/env python

import ROOT
import io,os,sys
import optparse
import commands
from UserCode.HGCanalysis.PlotUtils import *
from UserCode.HGCanalysis.HGCTree2Workspace import *

"""
"""
def adaptWorkspaceForEMCalibration(opt):
    
    wsUrl             = opt.wsUrl

    #prepare workspace (if needed) and output
    if wsUrl is None :

        #init weighting scheme
        matFurl='%s/src/UserCode/HGCanalysis/data/HGCMaterialOverburden.root'%os.environ['CMSSW_BASE']
        matF=ROOT.TFile.Open(matFurl)
        if opt.lambdaWeighting:
            matParamBeforeHGCGr=matF.Get("lambdaOverburden")
            weightingScheme={
                'EE': [([1, 1 ], matParamBeforeHGCGr, 0.010, 1.0),
                       ([2, 11], None,                0.036, 1.0),
                       ([12,21], None,                0.043, 1.0),
                       ([22,30], None,                0.056, 1.0),
                       ([31,31], None,                0.338, 1.0),
                       ([32,42], None,                0.273, 1.0),
                       ([43,54], None,                0.475, 1.0)]
                }
        else:
            matParamBeforeHGCGr=matF.Get("x0Overburden")
            weightingScheme={
                'EE': [([1, 1 ], matParamBeforeHGCGr, 0.080, 1.0),
                       ([2, 11], None,                0.620, 1.0),
                       ([12,21], None,                0.809, 1.0),
                       ([22,30], None,                1.239, 1.0),
                       ([31,31], None,                3.580, 1.0),
                       ([32,42], None,                3.103, 1.0),
                       ([43,54], None,                5.228, 1.0)]
                }
        matF.Close()

        #prepare the workspace and get new url
        wsUrl=prepareWorkspace(url=opt.input,weightingScheme=weightingScheme,vetoTrackInt=opt.vetoTrackInt,vetoHEBLeaks=opt.vetoHEBLeaks,treeVarName=opt.treeVarName)
        
    return wsUrl

"""
draws the correlation shower leakage, and phi resolution
"""
def computeCompensationWeights(enRanges,etaRanges,ws,outDir):

    print '[computeCompensationWeights] will look at correlations between shower leakage, phi and resolutions'

    fitFunc=ROOT.TF1('phicorr','TMath::Abs(x-[0])>0.03 ? [3]*pow([0]+0.03,2)+[2]*([0]+0.03)+[1] : [3]*pow(x-[0],2)+[2]*TMath::Abs(x-[0])+[1]',0,ROOT.TMath.Pi()/9.)
    fitFunc.SetParameter(0,ROOT.TMath.Pi()/18.)
    canvas=ROOT.TCanvas('c','c',1000,500)
    for ien in xrange(0,len(enRanges)):

        genEn_min=enRanges[ien][0]
        genEn_max=enRanges[ien][1]
        genEn_mean=0.5*(genEn_max+genEn_min)

        tailFracH = ROOT.TH2F('tailfrac%d'%ien,';Tail fraction;#DeltaE/E',10,0,0.5,50,0.5,1.5)
        tailFracH.Sumw2()
        tailFracH.SetDirectory(0)
        phiH      = ROOT.TH2F('phi%d'%ien,';#phi [rad];#DeltaE/E',25,0,ROOT.TMath.Pi()/9.,50,0.5,1.5)
        phiH.Sumw2()
        phiH.SetDirectory(0)

        #fill histograms and prepare for quantile computation
        redIncData=ws.data('data').reduce('en>=%f && en<=%f'%(genEn_min,genEn_max))
        nEntries=redIncData.numEntries()
        if nEntries<10: continue
        entryWeight=1./float(nEntries)

        recEn_mean=redIncData.mean(ws.var('en_EE'))
        for ientry in xrange(0,nEntries):
            entryVars=redIncData.get(ientry)
            phi=entryVars.find('phi').getVal()
            nphi=ROOT.TMath.Floor(ROOT.TMath.Abs(9*phi/ROOT.TMath.Pi()))
            phi=ROOT.TMath.Abs(phi)-nphi*ROOT.TMath.Pi()/9.
            tailFrac=entryVars.find('tailfrac_EE').getVal()
            recEn=entryVars.find('en_EE').getVal()
            phiH.Fill(phi,recEn/recEn_mean,entryWeight)
            tailFracH.Fill(tailFrac,recEn/recEn_mean,entryWeight)
        
        #show profile
        canvas.Clear()
        canvas.Divide(2,1)   
        p=canvas.cd(1)
        p.SetRightMargin(0.15)
        tailFracH.Draw('colz')
        MyPaveText('#bf{CMS} #it{simulation} Energy=%d'%genEn_mean)
        p=canvas.cd(2)
        p.SetRightMargin(0.15)
        phiH.Draw('colz')
        phiprofH=phiH.ProfileX('phiprof_%d'%ien)
        phiprofH.SetMarkerStyle(20)
        phiprofH.Draw('e1same')
        phiprofH.Fit(fitFunc,'MRQ+','same')        
        canvas.Modified()
        canvas.Update()
        canvas.SaveAs('%s/em_%d_profile.png'%(outDir,ien))
        canvas.SaveAs('%s/em_%d_profile.C'%(outDir,ien))

    return fitFunc

"""
Steer the calibration study
"""
def runCalibrationStudy(opt):

    #get the workspace
    wsUrl = adaptWorkspaceForEMCalibration(opt)
    print 'Workspace for e.m. calibration, stored @ %s'%wsUrl
    wsOutF=ROOT.TFile.Open(wsUrl)
    ws=wsOutF.Get('w')
    wsOutF.Close()

    #prepare output dir
    outDir=os.path.dirname(wsUrl)
    

    #init phase space regions of interest
    etaRanges = [[1.5,1.6],[1.6,1.75],[1.75,2.0],[2.0,2.25],[2.25,2.5],[2.5,2.75],[2.75,2.95]]
    enRanges  = [[9,11],[19,21],[39,41],[49,51],[74,75],[99,101],[149,151],[249,251],[499,501]] 

    #readout calibration
    weightTitles={'simple':'Simple sum'}
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
            calibMapRes[wType]=calibF.Get('calib_3_%s_res'%wType).Clone() #use the 2.0-2.25 eta range
        calibF.Close()
    except:
        print 'No calibration will be applied'

    #compensation as function of critical variables
    #phiCorrGr=computeCompensationWeights(enRanges,etaRanges,ws,outDir)

    #create dataset for calibration
    uncalibDataVars=ROOT.RooArgSet(ws.var('en'), ws.var('eta'), ws.var('phi'),ws.var('tailfrac_EE'))
    for wType in weightTitles: 
        ws.factory('%sEn[0,0,999999999]'%wType)
        ws.var('%sEn'%wType).SetTitle( '%s'%weightTitles[wType] )
        uncalibDataVars.add( ws.var('%sEn'%wType ) )
    getattr(ws,'import')( ROOT.RooDataSet('data_uncalib_final','data_uncalib_final',uncalibDataVars) )
    for ientry in xrange(0,ws.data('data').numEntries()):

        entryVars=ws.data('data').get(ientry)
        newEntry=ROOT.RooArgSet()

        for baseVar in ['en','eta','phi','tailfrac_EE']: 
            ws.var(baseVar).setVal( entryVars.find(baseVar).getVal() )
            newEntry.add( ws.var(baseVar) )


        #check phi
        phi=entryVars.find('phi').getVal()
        nphi=ROOT.TMath.Floor(ROOT.TMath.Abs(9*phi/ROOT.TMath.Pi()))
        phi=ROOT.TMath.Abs(phi)-nphi*ROOT.TMath.Pi()/9.
        if ROOT.TMath.Abs(phi-ROOT.TMath.Pi()/18.)<0.03 : continue

        enEstimators={'simple':entryVars.find('en_EE').getVal()}
            
        for wType in weightTitles :

            ienVal=enEstimators[wType]

            #calibrated energy estimator
            if wType in calibMap:
                ienVal=calibMap[wType].GetX(ienVal)

            if wType in calibMapRes:
                resCorr=calibMapRes[wType].Eval(ienVal)
                ienVal*=(1-resCorr/100.)

            ws.var('%sEn'%wType).setVal(ienVal)
            newEntry.add(ws.var('%sEn'%wType))

        #all filled, add new row
        ws.data('data_uncalib_final').add(newEntry)


    #calibrate the energy estimators (split up in different energies and pseudo-rapidity ranges)
    nSigmasToFit=2.8 #2.8
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

            for wType in weightTitles:

                #prepare the fit to this slice
                vName = '%sEn'%wType
                v_mean, v_sigma    = redData.mean(ws.var(vName)), redData.sigma(ws.var(vName))
                v_min, v_max       = v_mean-2*v_sigma, v_mean+2*v_sigma
                v_fitMin, v_fitMax = v_mean-nSigmasToFit*v_sigma, v_mean+nSigmasToFit*v_sigma

                if v_mean<=0 or v_sigma<=0: 
                    print 'Something wrong for E in ',enRanges[iEnRange],' eta in ',etaRanges[iEtaRange],' with ',wType,' weights'
                    continue

                #define PDF
                fitName          = 'range%d%d_%s'%(iEtaRange,iEnRange,vName)            
                ws.var(vName).setRange(fitName,v_min,v_max)
                ws.var(vName).setRange('fit_%s'%fitName,v_fitMin, v_fitMax)
                ws.factory('RooCBShape::resol_%s(%s,mean_%s[%f,%f,%f],sigma_%s[%f,%f,%f],alpha_%s[0.001.,0.0001,20.0],n_%s[2,0.5,3])'%
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

                #save result
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

        for wType in weightTitles:
            calibGr[wType].Add(etaSliceCalibGr[wType],'p')
            resGr[wType].Add(etaSliceResGr[wType],'p')

    #derive calibration
    calibModel=ROOT.TF1('calibmodel',"[0]*x",0,1000)
    if opt.useSaturatedCalibModel:
        print 'Saturated calib model will be used'
        calibModel=ROOT.TF1('calibmodel','x<[4] ? [0]*x : [3]*((x-[1])/(x+[2])-([4]-[1])/([4]+[2]))+[0]*[4]',0,1000)
        calibModel.FixParameter(4,80)
    else:
        print 'Using linear approximation'
    calibModel.SetLineWidth(1)
    for wType in weightTitles:
        if opt.useSaturatedCalibModel:
            calibGr[wType].Fit(calibModel,'MER+','',0,1000)
        else:
            calibGr[wType].Fit(calibModel,'MER+','',0,100)
        calibGr[wType].GetFunction(calibModel.GetName()).SetRange(0,1000)
        calibGr[wType].GetFunction(calibModel.GetName()).SetLineColor(calibGr[wType].GetListOfGraphs().At(0).GetLineColor())

    #show results
    resCorrectionGr,resCalibGr=showCalibrationCurves(calibGr=calibGr,calibRanges=etaRanges,outDir=outDir,calibPostFix=calibPostFix)
    showResolutionCurves(resGr=resGr,outDir=outDir,calibPostFix=calibPostFix)

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
    parser.add_option('-i',      '--in' ,      dest='input',        help='Input file',                                     default=None)
    parser.add_option('-w',      '--ws' ,      dest='wsUrl',        help='Workspace file',                                 default=None)
    parser.add_option('-c',      '--calib' ,   dest='calibUrl',     help='Calibration file',                               default=None)
    parser.add_option('--vetoTrackInt',        dest='vetoTrackInt', help='flag if tracker interactions should be removed', default=False, action="store_true")
    parser.add_option('--lambdaWeighting',     dest='lambdaWeighting', help='flag if layers should be weighted by lambda', default=False, action="store_true")
    parser.add_option('--vetoHEBLeaks',        dest='vetoHEBLeaks',  help='flag if HEB leaks are allowed',                                  default=False, action='store_true')
    parser.add_option('--useSaturatedCalibModel',  dest='useSaturatedCalibModel', help='use a calibration model which saturates at high energy', default=False, action='store_true')
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
