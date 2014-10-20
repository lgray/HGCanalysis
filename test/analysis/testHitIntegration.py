#!/usr/bin/env python

import ROOT
import sys
import optparse
import commands
import array
import glob
import os
import numpy as np
from ROOT import *

from UserCode.HGCanalysis.PlotUtils import *

"""
helper for hit integration
"""
class HitIntegrator:
  
  """
  book counters
  """
  def __init__(self,nsd,nlay):
    self.nsd=nsd
    self.nlay=[0]*nsd
    for isd in xrange(0,self.nsd): 
      #HEB sub-sectors cut in 1/2
      self.nlay[isd]=nlay[isd]
      if isd==2: self.nlay[isd]=nlay[isd]/2 
    self.totalEnInLayer=[]
    self.totalEn=[]
    self.totalEnIT=[]
    self.totalHits=[]
    self.totalADC=[]
    self.alpha=[]
    for isd in xrange(0,self.nsd) : 
      self.totalEnInLayer.append( [0]*self.nlay[isd] )
      self.totalEn.append( [0]*self.nlay[isd] )
      self.totalEnIT.append( [0]*self.nlay[isd] )
      self.totalADC.append( [0]*self.nlay[isd] )
      self.totalHits.append( [0]*self.nlay[isd] )
      self.alpha.append( [0]*self.nlay[isd] )

  """
  reset counters
  """
  def reset(self):
    for isd in xrange(0,self.nsd) : 
      for ilay in xrange(0,self.nlay[isd]) :
        self.totalEnInLayer[isd][ilay] = 0
        self.totalEn[isd][ilay]        = 0
        self.totalEnIT[isd][ilay]      = 0
        self.totalHits[isd][ilay]      = 0
        self.totalADC[isd][ilay]       = 0
        self.alpha[isd][ilay]          = 0

  """
  increment counters
  """
  def integrate(self,isd,ilay,simHitEn,simHitEnIT,hitADC,deltaR,probeCone):
    self.totalEnInLayer[isd][ilay-1]   += simHitEn
    if deltaR>probeCone: return
    self.totalEn[isd][ilay-1]   += simHitEn
    self.totalEnIT[isd][ilay-1] += simHitEnIT
    self.totalADC[isd][ilay-1]  += hitADC
    if hitADC>1 : self.totalHits[isd][ilay-1] += 1
    if deltaR>0 : self.alpha[isd][ilay-1] += hitADC/deltaR


"""
checks HGC hits for a particle gun, against extrapolated tracks or gen particles
a control region in the opposite pseudo-rapidity is checked as well
"""
def testHitIntegration(urlList,probeCone,useTrackAsRef,outUrl):

  print '[testHitIntegration] with %d files, using dR=%3.2f'%(len(urlList),probeCone)
  if useTrackAsRef : print ' reconstructed track will be used as a reference'

  #general root formats
  customROOTstyle()
  gROOT.SetBatch(False)
  gStyle.SetPalette(51)
  ROOT.TH1.SetDefaultSumw2(True)

  #get the trees
  HGC=ROOT.TChain('analysis/HGC')
  for url in urlList: HGC.Add(url)
  HGC.GetEntry(0)

  #init hit integrators
  signalHitIntegrator=HitIntegrator(HGC.nsd,HGC.nlay)
  ctrlHitIntegrator=[None]*10
  for iphi in xrange(0,10): ctrlHitIntegrator[iphi]=HitIntegrator(HGC.nsd,HGC.nlay)

  #prepare histos
  histos={}
  if useTrackAsRef :
    histos['tk_ptres']  = ROOT.TH1F('tk_ptres',  ';(p_{T}^{track}-p_{T}^{gen})/p_{T}^{gen} [GeV];Events;',50,-0.255,0.245)
    histos['tk_etares'] = ROOT.TH1F('tk_etares',';(#eta^{track}-#eta^{gen});Events;',                    50,-0.0102,0.0098)
    histos['tk_phires'] = histos['tk_etares'].Clone('tk_phires')
    histos['tk_phires'].GetXaxis().SetTitle('(#phi^{track}-#phi^{gen}) [rad]')

  for sd in xrange(0,signalHitIntegrator.nsd):
    for lay in xrange(0,signalHitIntegrator.nlay[sd]):
      for reg in ['','ctrl_']:
        pfix='sd%d_lay%d_%s'%(sd,lay+1,reg)
        #hit alignment
        histos[pfix+'hitwgtdx']      = ROOT.TH1F(pfix+'hitwgtdx',';x^{hit}-x^{track} [mm];ADC weighted hits;',50,-102,98)
        histos[pfix+'hitwgtdy']      = ROOT.TH1F(pfix+'hitwgtdy',';y^{hit}-y^{track} [mm];ADC weighted hits;',50,-102,98)
        histos[pfix+'hitwgtdz']      = ROOT.TH1F(pfix+'hitwgtdz',';z^{hit}-z^{track} [mm];ADC weighted hits;',50,-102,98)

        maxSimHitEn=1000
        if sd==2 : maxSimHitEn=10000
        histos[pfix+'simhiten']      = ROOT.TH2F(pfix+'simhiten',';Energy [keV]; Events',50,0,maxSimHitEn,15,1.5,3.0)
        histos[pfix+'simhitenit']    = histos[pfix+'simhiten'].Clone(pfix+'simhitenit')
        histos[pfix+'hitadc']        = ROOT.TH2F(pfix+'hitadc',';ADC counts; Events',100,-0.5,99.5,15,1.5,3.0)
        histos[pfix+'nhits']         = ROOT.TH2F(pfix+'nhits',';Number of hits; Events',100,0,100,15,1.5,3.0)
        histos[pfix+'hitalpha']      = ROOT.TH2F(pfix+'hitalpha',';#alpha=log #Sigma (ADC_{i}/#Delta R_{i}) x #Theta(#DeltaR_{i}<%3.1f);Events'%(probeCone),50,-5,15,15,1.5,3.0)

  for hname in histos: histos[hname].SetDirectory(0)

  #prepare ntuple
  print '[testHitIntegration] histograms stored @ %s'%outUrl
  fOut=ROOT.TFile(outUrl,'RECREATE')
  ntupleVarNames='genEn:genEta:genPhi'
  for sd in xrange(0,signalHitIntegrator.nsd):
    for lay in xrange(0,signalHitIntegrator.nlay[sd]):
      for ivar in ['toten','en','adc','puadc','hits','puhits']:
        ntupleVarNames += ':%s_s%dl%d'%(ivar,sd,lay+1)
  fOut.cd()
  ntupleVars=[0]*len(ntupleVarNames.split(':'))  
  output_tuple = ROOT.TNtuple("HGC","HGC",ntupleVarNames)
  output_tuple.SetDirectory(fOut)

  #loop over events
  for iev in xrange(0,HGC.GetEntries()):
    HGC.GetEntry(iev)

    sys.stdout.write( '\r[testHitIntegration] status [%d/%d]'%(iev,HGC.GetEntries()))
       
    if HGC.ngen!=1 : continue
    ntupleVars[0]=HGC.gen_en[0]
    ntupleVars[1]=HGC.gen_eta[0]
    ntupleVars[2]=HGC.gen_phi[0]
    
    #track resolution
    if useTrackAsRef:
      if HGC.ntk!=1 : continue
      histos['tk_ptres'].Fill(HGC.tk_pt[0]/HGC.gen_pt[0]-1)
      histos['tk_etares'].Fill(HGC.tk_eta[0]-HGC.gen_eta[0])
      histos['tk_phires'].Fill(HGC.tk_phi[0]-HGC.gen_phi[0])

    #hit analysis
    if HGC.nhits==0 : continue
    signalHitIntegrator.reset()
    for chi in ctrlHitIntegrator : chi.reset()
    for n in xrange(0,HGC.nhits):

      # hit coordinates
      sdType          = HGC.hit_type[n]
      layer           = HGC.hit_layer[n]
      hit_eta         = HGC.hit_eta[n]
      #if abs(abs(hit_eta)-abs(HGC.gen_eta[0]))>probeCone*1.5 : continue # not worth looking into this ones
      hit_phi         = HGC.hit_phi[n]
      etaEnCorrection = ROOT.TMath.Abs(ROOT.TMath.TanH(hit_eta))
      simHitEn        = HGC.hit_edep[n]*etaEnCorrection
      #cf. http://root.cern.ch/phpBB3/viewtopic.php?t=9457
      simHitEnIT      = HGC.hit_edep_sample[n*5+0]*etaEnCorrection
      hitADC          = HGC.digi_adc[n]*etaEnCorrection
      hit_x           = HGC.hit_x[n]
      hit_y           = HGC.hit_y[n]
      hit_z           = HGC.hit_z[n]


      isCtrlRegion,regionType=False,''
      if hit_eta*HGC.gen_eta[0]<0 : regionType,isCtrlRegion='ctrl_',True
        
      pfix='sd%d_lay%d_%s'%(sdType,layer,regionType)

      #by default use the generator level particle as reference
      refEta = HGC.gen_eta[0]
      if isCtrlRegion : refEta = -1*refEta
      refPhi = HGC.gen_phi[0]

      #if track is to be used, find the extrapolation of the track from closest z (only needed for HEB)
      if useTrackAsRef : 
        iForMinDZ=layer-1
        if sdType>0 : iForMinDZ += HGC.nlay[0]
        if sdType>1 : 
          minDZ=999999.
          iForMinDZ=-1
          for itkext in xrange(HGC.nlay[0]+HGC.nlay[1],66):
            #cf. http://root.cern.ch/phpBB3/viewtopic.php?t=9457
            tk_z=HGC.tk_extrapol_z[0*HGC.ntk+itkext]
            if isCtrlRegion : tk_z = -1*tk_z
            dZ=ROOT.TMath.Abs(tk_z-hit_z)
            if dZ > minDZ : continue
            minDZ=dZ
            iForMinDZ=itkext

        #if no extrapolation found, move to the next event
        if iForMinDZ<0 : continue

        #cf. http://root.cern.ch/phpBB3/viewtopic.php?t=9457
        tk_x = HGC.tk_extrapol_x[0*HGC.ntk+iForMinDZ]
        tk_y = HGC.tk_extrapol_y[0*HGC.ntk+iForMinDZ]
        tk_z = HGC.tk_extrapol_z[0*HGC.ntk+iForMinDZ]
        histos[pfix+'hitwgtdx'].Fill(hit_x-tk_x,hitADC)
        histos[pfix+'hitwgtdy'].Fill(hit_y-tk_y,hitADC)
        histos[pfix+'hitwgtdz'].Fill(hit_z-tk_z,hitADC)
        
        #update reference eta and phi
        tk_z   = hit_z
        tk_rho = ROOT.TMath.Sqrt(tk_x*tk_x+tk_y*tk_y+tk_z*tk_z)
        if tk_rho==0 or tk_z==0 : continue
        refPhi = ROOT.TMath.ATan2(tk_y,tk_x)
        refEta = 0
        if tk_rho>tk_z: refEta=0.5*ROOT.TMath.Log( (tk_rho+tk_z)/(tk_rho-tk_z) )


      #distance in phase space: for the control region sample every 36 deg in phi 
      if isCtrlRegion :
        for iphi in xrange(0,10) :
          ctrlPhi = refPhi+iphi*ROOT.TMath.Pi()*36./180.  
          deltaPhi = ROOT.TVector2.Phi_mpi_pi(ctrlPhi-hit_phi)
          deltaEta = refEta-hit_eta
          deltaR   = ROOT.TMath.Sqrt(deltaEta*deltaEta+deltaPhi*deltaPhi)
          ctrlHitIntegrator[iphi].integrate(sdType,layer,simHitEn,simHitEnIT,hitADC,deltaR,probeCone)
      else : 
        deltaPhi = ROOT.TVector2.Phi_mpi_pi(refPhi-hit_phi)
        deltaEta = refEta-hit_eta
        deltaR   = ROOT.TMath.Sqrt(deltaEta*deltaEta+deltaPhi*deltaPhi)
        signalHitIntegrator.integrate(sdType,layer,simHitEn,simHitEnIT,hitADC,deltaR,probeCone)
      

    #fill hit histograms
    ntupleIdx=3
    for sd in xrange(0,signalHitIntegrator.nsd):
      for lay in xrange(0,signalHitIntegrator.nlay[sd]):
        
        pfix='sd%d_lay%d_'%(sd,lay+1)
          
        #signal region
        histos[pfix+'simhiten'].Fill( signalHitIntegrator.totalEn[sd][lay],HGC.gen_eta[0] )
        histos[pfix+'simhitenit'].Fill( signalHitIntegrator.totalEnIT[sd][lay],HGC.gen_eta[0] )
        histos[pfix+'nhits'].Fill( signalHitIntegrator.totalHits[sd][lay],HGC.gen_eta[0] )
        histos[pfix+'hitadc'].Fill( signalHitIntegrator.totalADC[sd][lay],HGC.gen_eta[0] )
        if signalHitIntegrator.alpha[sd][lay]>0 : histos[pfix+'hitalpha'].Fill( ROOT.TMath.Log(signalHitIntegrator.alpha[sd][lay]),HGC.gen_eta[0] )

        #control region
        totalEn_ctrl=[]
        totalEnIT_ctrl=[]
        totalHits_ctrl=[]
        totalADC_ctrl=[]
        alpha_ctrl=[]
        for chi in ctrlHitIntegrator:
          totalEn_ctrl.append( chi.totalEn[sd][lay] )
          totalEnIT_ctrl.append( chi.totalEnIT[sd][lay])
          totalHits_ctrl.append( chi.totalHits[sd][lay] )
          totalADC_ctrl.append( chi.totalADC[sd][lay] )
          alpha_ctrl.append( chi.alpha[sd][lay] )

        medianTotalEn_ctrl   = np.median( totalEn_ctrl )
        medianTotalHits_ctrl = np.median(totalHits_ctrl)
        medianTotalADC_ctrl  = np.median(totalADC_ctrl)
        histos[pfix+'ctrl_simhiten'].Fill( medianTotalEn_ctrl,HGC.gen_eta[0] )
        histos[pfix+'ctrl_simhitenit'].Fill( np.median( totalEnIT_ctrl),HGC.gen_eta[0] )
        histos[pfix+'ctrl_nhits'].Fill( medianTotalHits_ctrl,HGC.gen_eta[0] )
        histos[pfix+'ctrl_hitadc'].Fill( medianTotalADC_ctrl,HGC.gen_eta[0] )
        medianAlpha=np.median(alpha_ctrl)
        if medianAlpha>0 : histos[pfix+'ctrl_hitalpha'].Fill( ROOT.TMath.Log(medianAlpha),HGC.gen_eta[0] )

        #order must match the one declared in ntupleVarNames
        ntupleVars[ntupleIdx]= signalHitIntegrator.totalEnInLayer[sd][lay]
        ntupleIdx+=1
        ntupleVars[ntupleIdx]= signalHitIntegrator.totalEn[sd][lay]
        ntupleIdx+=1
        ntupleVars[ntupleIdx]= signalHitIntegrator.totalADC[sd][lay]
        ntupleIdx+=1
        ntupleVars[ntupleIdx]= medianTotalADC_ctrl
        ntupleIdx+=1
        ntupleVars[ntupleIdx]= signalHitIntegrator.totalHits[sd][lay]
        ntupleIdx+=1
        ntupleVars[ntupleIdx]= medianTotalHits_ctrl
        ntupleIdx+=1

    #fill the ntuple
    output_tuple.Fill(array.array('f',ntupleVars))
  
        
  #write the result to a file
  fOut.cd()
  output_tuple.Write()
  for hname in histos:
    histos[hname].SetDirectory(fOut)
    histos[hname].Write()
  fOut.Close()


"""
checks the input arguments and steers the analysis
"""
def main():

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i',      '--in' ,      dest='input',    help='Input file',                                     default=None)
    parser.add_option('-o',      '--out' ,     dest='output',   help='Output file',                                    default='HitIntegrationAnalysis.root')
    parser.add_option('-r',      '--dR' ,      dest='dR',       help='DeltaR cone used for analysis (0.3 by default)', default=0.3,   type=float)
    parser.add_option('--useTrack',            dest='useTrack', help='If given, use track as reference',               default=False, action="store_true")
    (opt, args) = parser.parse_args()

    #check inputs
    if opt.input is None :
        parser.print_help()
        sys.exit(1)

    #run analysis
    urlList=[]
    if opt.input.find('/store')>=0:
      basename=os.path.basename(opt.input)
      basedir=os.path.dirname(opt.input)
      from UserCode.HGCanalysis.storeTools_cff import fillFromStore
      allFiles=fillFromStore(basedir)
      for f in allFiles:
        if f.find(basename)<0 : continue
        urlList.append(f)
    else:
      urlList=glob.glob('%s*.root'%opt.input)

    testHitIntegration(urlList=urlList,
                       probeCone=opt.dR,
                       useTrackAsRef=opt.useTrack,
                       outUrl=opt.output)

if __name__ == "__main__":
    main()
