#!/usr/bin/env python

import ROOT
import numpy as numpy
import json
import os

from UserCode.HGCanalysis.PlotUtils import *

"""
Converts all to a workspace and computes optimized weights
"""
def prepareWorkspace(url,weightingScheme,treeVarName,vetoTrackInt,vetoHEBLeaks=False):

    #prepare the workspace
    outUrl=os.path.basename(url).replace('.root','')
    os.system('mkdir -p %s'%outUrl)
    ws=ROOT.RooWorkspace("w")
    dsVars=ROOT.RooArgSet( ws.factory('en[0,0,9999999999]'), 
                           ws.factory('eta[1.5,1.45,3.1]'), 
                           ws.factory('phi[0,-3.2,3.2]') )
    for key in weightingScheme:
        dsVars.add( ws.factory('en_%s[0,0,99999999]'%key) )
        dsVars.add( ws.factory('tailfrac_%s[0,0,99999999]'%key) )
        dsVars.add( ws.factory('c_%s[0,0,99999999]'%key) )
        dsVars.add( ws.factory('rho_%s[0,0,99999999]'%key) )
    getattr(ws,'import')( ROOT.RooDataSet('data','data',dsVars) )

    #read all to a RooDataSet
    if url.find('/store/')>=0: url='root://eoscms//eos/cms/%s'%url
    fin=ROOT.TFile.Open(url)
    HGC=fin.Get('analysis/HGC')
    simStep=treeVarName.split('_')[1]
    print 'Reading hits from HGC tree found in %s'%url
    print 'Status will be updated every 1k events read / %dk events available'%(HGC.GetEntriesFast()/1000)
    for entry in xrange(0,HGC.GetEntriesFast()+1):
        HGC.GetEntry(entry)
        if entry%1000==0 : drawProgressBar(float(entry)/float(HGC.GetEntriesFast()))

        #veto interactions in the tracker
        if vetoTrackInt and HGC.hasInteractionBeforeHGC : continue

        #check the amount of energy deposited in the back HEB
        if vetoHEBLeaks:
            sumBackHEB=0
            for ilayer in [51,52,53]:
                sumBackHEB+=(getattr(HGC,treeVarName))[ilayer-1]
            if sumBackHEB>3: continue

        #require generated particle to be in the endcap
        genEn=HGC.genEn
        genEta=ROOT.TMath.Abs(HGC.genEta)
        genPhi=HGC.genPhi
        if genEta<1.4 or genEta>3.0 : continue
        ws.var('en').setVal(genEn)
        ws.var('eta').setVal(genEta)
        ws.var('phi').setVal(genPhi)

        #start a new entry
        newEntry=ROOT.RooArgSet(ws.var('en'), ws.var('eta'),  ws.var('phi'))

        #get the relevant energy deposits and add new row
        for subDet,subDetScheme in weightingScheme.iteritems():
            
            totalEn, totalEnTail, nhits, nhits10mip, nhitsavg = 0, 0, 0, 0, 0
            for subDetRange in subDetScheme:

                #decode corrections to apply to this sub-detector range
                integRange             = subDetRange[0]
                etaDepWeight           = 0.0
                if not (vetoTrackInt or subDetRange[1] is None):
                    etaDepWeight       = subDetRange[1].Eval(genEta)
                geomCorrection         = 1. #/ROOT.TMath.TanH(genEta)
                weight                 = subDetRange[2]
                scaleCorrections       = subDetRange[3]

                #integrate layers applying corrections
                totalEnTail=0
                for ilay in xrange(integRange[0],integRange[1]+1):
                    ien         = (weight*geomCorrection+etaDepWeight)*(getattr(HGC,treeVarName))[ilay-1]
                    if not (scaleCorrections is None):
                        ien=scaleCorrections[0].GetX(ien)
                        ien*=(1-scaleCorrections[1].Eval(ien)/100.)
                    totalEn    += ien
                    if ilay>integRange[1]-3: totalEnTail += ien
                    nhits      += (getattr(HGC,'nhits_%s'%(simStep)))[ilay-1]
                    nhits10mip += (getattr(HGC,'nhits10mip_%s'%(simStep)))[ilay-1]
                    nhitsavg   += (getattr(HGC,'nhitsavg_%s'%(simStep)))[ilay-1]
                        
            #add to the set of variables the values computed for this subdetector
            ws.var('en_%s'%subDet).setVal(totalEn)
            newEntry.add(ws.var('en_%s'%subDet))

            tailfrac=0
            if totalEn>0: tailfrac=totalEnTail/totalEn
            ws.var('tailfrac_%s'%subDet).setVal(tailfrac)
            newEntry.add(ws.var('tailfrac_%s'%subDet))

            cfrac = 0
            if nhits>1:
                cfrac = float(nhits-nhits10mip)/float(nhits-nhitsavg)
            ws.var('c_%s'%subDet).setVal(cfrac)
            newEntry.add(ws.var('c_%s'%subDet))

            rho = 0
            volume=getattr(HGC,'totalVolume%s_%s'%(subDet,simStep))
            if volume>0:
                rho   = totalEn/volume
            ws.var('rho_%s'%subDet).setVal(rho)
            newEntry.add(ws.var('rho_%s'%subDet))

        ws.data('data').add( newEntry )

    fin.Close()

    #all done, write to file
    wsFileUrl='%s/workspace.root'%outUrl
    ws.writeToFile(wsFileUrl,True)
    print '\nCreated the analysis RooDataSet with %d events, stored @ %s'%(ws.data('data').numEntries(),wsFileUrl)
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
def showCalibrationCurves(calibGr,calibRanges,outDir,calibPostFix) :
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
    resCorrectionGr={}
    igr=0
    for wType in calibGr:
        if igr==0:
            calibGr[wType].Draw('a')
            calibGr[wType].GetXaxis().SetTitleSize(0)
            calibGr[wType].GetXaxis().SetLabelSize(0)
            calibGr[wType].GetYaxis().SetTitle('Reconstructed energy [GeV]')
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
        if ffunc.GetNpar()==1:
            MyPaveText('#it{%s} : %3.4f E_{rec} '%(calibGr[wType].GetTitle(),1./ffunc.GetParameter(0)),
                       0.2,0.93-igr*0.05,0.4,0.90-igr*0.05)
        else:
            txt='#it{%s}\\ '%calibGr[wType].GetTitle()
            for ip in xrange(0,ffunc.GetNpar()):
                txt += '[%d]=%3.2f '%(ip,ffunc.GetParameter(ip))
            MyPaveText(txt,0.2,0.93-igr*0.1,0.4,0.85-igr*0.1)

        #compute the residuals
        resCalibGr[wType]=ROOT.TMultiGraph()
        resCorrectionGr[wType]=None
        for gr in calibGr[wType].GetListOfGraphs() :

            newGr=gr.Clone('%s_res'%gr.GetName())
            newGr.Set(0)
            xval, yval = ROOT.Double(0), ROOT.Double(0)
            for ip in xrange(0,gr.GetN()):
                gr.GetPoint(ip,xval,yval)
                xval_error=gr.GetErrorX(ip)
                yval_error=gr.GetErrorY(ip)
                xCalib=ffunc.GetX(yval)
                xCalib_err=ROOT.TMath.Max(ROOT.TMath.Abs(ffunc.GetX(yval+yval_error)-xCalib),ROOT.TMath.Abs(ffunc.GetX(yval-yval_error)-xCalib))
                newGr.SetPoint(ip,xval,100*(xCalib/xval-1))
                newGr.SetPointError(ip,0,100*xCalib_err/xval)
            resCalibGr[wType].Add(newGr,'p')

            #linear approximation to residuals
            newGr.Fit('pol0','QME0+')
            if resCorrectionGr[wType] is None:
                resCorrectionGr[wType]=gr.Clone('%s_calib_residuals'%wType)
                resCorrectionGr[wType].SetTitle(calibGr[wType].GetTitle())
                resCorrectionGr[wType].Set(0)
            ip=resCorrectionGr[wType].GetN()
            calibXmin=calibRanges[ip][0]
            calibXmax=calibRanges[ip][1]
            resCorrectionGr[wType].SetPoint(ip,0.5*(calibXmax+calibXmin),newGr.GetFunction('pol0').GetParameter(0)/100.)
            resCorrectionGr[wType].SetPointError(ip,0.5*(calibXmax-calibXmin),newGr.GetFunction('pol0').GetParError(0)/100.)

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

    #now show residual calibration
    canvas.Clear()
    igr=0
    leg=ROOT.TLegend(0.2,0.7,0.9,0.9)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)
    for wType in resCorrectionGr:
        if igr==0:
            resCorrectionGr[wType].Draw('ap')
            resCorrectionGr[wType].GetXaxis().SetTitle('Pseudo-rapidity')
            resCorrectionGr[wType].GetYaxis().SetTitle('Residual correction')
            resCorrectionGr[wType].GetYaxis().SetTitleOffset(0.9)
            resCorrectionGr[wType].GetYaxis().SetTitleSize(0.05)
            resCorrectionGr[wType].GetYaxis().SetLabelSize(0.04)
            resCorrectionGr[wType].GetXaxis().SetTitleOffset(0.9)
            resCorrectionGr[wType].GetXaxis().SetTitleSize(0.05)
            resCorrectionGr[wType].GetXaxis().SetLabelSize(0.04)
            resCorrectionGr[wType].GetYaxis().SetRangeUser(-1.5*resCorrectionGr[wType].GetYaxis().GetXmax(),resCorrectionGr[wType].GetYaxis().GetXmax()*1.5)
        else:
            resCorrectionGr[wType].Draw('p')
        leg.AddEntry(resCorrectionGr[wType],resCorrectionGr[wType].GetTitle(),"p")
        igr+=1
    leg.Draw()
    MyPaveText('#bf{CMS} #it{simulation}').SetTextSize(0.04)
    canvas.cd()
    canvas.Modified()
    canvas.Update()
    canvas.SaveAs('%s/rescalib%s.png'%(outDir,calibPostFix))
    
    return resCorrectionGr,resCalibGr
    

"""
shows a set of resolution curves
"""
def showResolutionCurves(resGr,outDir,calibPostFix,model=0) :

    resolModel=None
    if model==0  : resolModel=ROOT.TF1('resolmodel',"sqrt([0]*[0]/x+[1]*[1])",0,1000)
    else : 
        resolModel=ROOT.TF1('resolmodel',"sqrt([0]*[0]/x+[1]*[1]+[2]*[2]/(x*x))",0,1000)
        resolModel.SetParameter(2,0)
        resolModel.SetParLimits(2,0.05,0.1)
    resolModel.SetParameter(0,0.2);
    resolModel.SetParLimits(0,0,2);
    resolModel.SetParameter(1,0);
    resolModel.SetParLimits(1,0,1.0);

    canvas=ROOT.TCanvas('c','c',500,500)
    canvas.cd()
    leg=ROOT.TLegend(0.75,0.5,0.9,0.95)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.03)

    igr=0
    pt=[]
    for wType in resGr:
        if igr==0:
            resGr[wType].Draw('a')
            resGr[wType].GetXaxis().SetTitleSize(0.05)
            resGr[wType].GetXaxis().SetLabelSize(0.05)
            resGr[wType].GetXaxis().SetTitle('Generated energy [GeV]')
            resGr[wType].GetYaxis().SetTitle('#sigma_{E} / E')
            resGr[wType].GetYaxis().SetTitleOffset(1.1)
            resGr[wType].GetYaxis().SetTitleSize(0.05)
            resGr[wType].GetYaxis().SetLabelSize(0.04)
            resGr[wType].GetYaxis().SetRangeUser(0.,0.25)
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
        if model==0 :
            pt.append( MyPaveText('#it{%s} :  %3.4f#scale[0.8]{/#sqrt{E}} #oplus %3.4f'%(resGr[wType].GetTitle(),sigmaStoch,sigmaConst),
                                  0.2,0.93-igr*0.05,0.4,0.90-igr*0.05) )
        else :
            sigmaNoise = ffunc.GetParameter(2)
            pt.append( MyPaveText('#it{%s} :  %3.4f#scale[0.8]{/#sqrt{E}} #oplus %3.4f#scale[0.8]{/E} #oplus %3.4f'%(resGr[wType].GetTitle(),sigmaStoch,sigmaNoise,sigmaConst),
                                  0.2,0.93-igr*0.05,0.4,0.90-igr*0.05) )
        pt[igr-1].SetTextColor(lcol)
        pt[igr-1].SetTextSize(0.03)

    leg.Draw()
    MyPaveText('#bf{CMS} #it{simulation}').SetTextSize(0.04)

    canvas.Modified()
    canvas.Update()
    canvas.SaveAs('%s/resol%s.png'%(outDir,calibPostFix))
    canvas.SaveAs('%s/resol%s.C'%(outDir,calibPostFix))
