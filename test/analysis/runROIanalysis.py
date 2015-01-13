#!/usr/bin/env python

import ROOT
from UserCode.HGCanalysis.PlotUtils import *
from array import array
import os,sys
import optparse
import commands
import numpy

"""
scales rec hit energy (MIP) to e.m. energy scale
"""
def scaleToEM(edep,layer,eta):

    weight=0.01
    if layer>=2 and layer<=11    : weight=0.036
    elif layer>=12 and layer<=21 : weight=0.043
    elif layer>=22 and layer<=30 : weight=0.056
    elif layer==31               : weight=0.338
    elif layer>=32 and layer<=42 : weight=0.273
    elif layer>42:                 weight=0.475

    #weight *= ROOT.TMath.TanH( ROOT.TMath.Abs(eta) )

    if layer<=30               : edep=(edep*0.2339+0.1778)*1.00000*weight
    if layer>30 and layer<=42  : edep=(edep*0.1828+0.9601)*1.29208*weight
    else                       : edep=(edep*0.2471+1.4563)*1.29208*1.0535*weight

    return edep


"""
Loops over the trees and profiles the showers
"""
def runGenLevelAnalysis(opt):

    #create output tuple
    #if opt.output!='./' : 
    #    os.system('mkdir -p %s'%opt.output)
    fOut = ROOT.TFile("%s/ROIanalysis.root"%opt.output, "RECREATE")
    fOut.cd()
    subDets=['EE','FH','BH']
    roiTypes=['dr01','dr02','dr03','dr04','2x2','3x3','5x5']
    keysForHit=['raw','em','em_wgt','em_wgt_t']
    varNames={'genEn':numpy.zeros(1,dtype=float),
              'genEta':numpy.zeros(1,dtype=float),
              'genPhi':numpy.zeros(1,dtype=float),
              'genEmFrac':numpy.zeros(1,dtype=float),
              'genHadFrac':numpy.zeros(1,dtype=float)}
    for roiKey in roiTypes:
        for det in subDets:
            for integKey in keysForHit:
                varNames['%s_%s_%s_en'%(det,roiKey,integKey)]=numpy.zeros(1,dtype=float)
                varNames['%s_%s_%s_eta'%(det,roiKey,integKey)]=numpy.zeros(1,dtype=float)
                varNames['%s_%s_%s_phi'%(det,roiKey,integKey)]=numpy.zeros(1,dtype=float)
    ntuple=ROOT.TTree('ROITuple','ROITuple')
    for var in varNames:
        ntuple.Branch(var, varNames[var], '%s/D'%var)
    ntuple.SetAutoSave(100000)

    #book histograms
    histos={
        'cutflow'    : ROOT.TH1F('cutflow','Selection;Events',5,0,5),
        'mjj'        : ROOT.TH1F('mjj',';Dijet invariant mass [GeV];Events',20,0,2000),
        'detajj'     : ROOT.TH1F('detajj',';#Delta#eta;Events',10,0,7),
        'etaLead'    : ROOT.TH1F('etaLead',';Pseudo-rapidity;Events',20,0,5),
        'etaTrailer' : ROOT.TH1F('etaTrailer',';Pseudo-rapidity;Events',20,0,5),
        'enLead'     : ROOT.TH1F('enLead',';Energy [GeV];Events',20,0,3000),
        'enTrailer'  : ROOT.TH1F('enTrailer',';Energy [GeV];Events',20,0,3000),
        'emFrac'     : ROOT.TH2F('emFrac',';Energy [GeV];e.m. fraction;Jets',20,0,1000,10,0,1),
        'hadFrac'    : ROOT.TH2F('hadFrac',';Energy [GeV];hadronic fraction;Jets',20,0,1000,10,0,1),
        'invFrac'    : ROOT.TH2F('invFrac',';Energy [GeV];invisible fraction;Jets',20,0,1000,10,0,1)      
    }
    histos['cutflow'].GetXaxis().SetBinLabel(1,'Total')
    histos['cutflow'].GetXaxis().SetBinLabel(2,'#geq 2j')
    histos['cutflow'].GetXaxis().SetBinLabel(3,'M_{jj}, #Delta#eta')
    histos['cutflow'].GetXaxis().SetBinLabel(4,'1j in HGC')
    histos['cutflow'].GetXaxis().SetBinLabel(5,'2j in HGC')
    for var in histos:
        histos[var].Sumw2()
        histos[var].SetDirectory(fOut)

    #fill histos
    nEvts=0
    fIn=ROOT.TFile.Open(opt.input)
    HGC=fIn.Get('analysis/HGC')
    for i in xrange(0,HGC.GetEntriesFast()):
        HGC.GetEntry(i)

        if i%10==0 : drawProgressBar(float(i)/float(HGC.GetEntriesFast()))
            
        histos['cutflow'].Fill(0)

        #find the leading-energy jets
        jetP4=[]
        nLepMatches=0
        for nj in xrange(0,HGC.njgen):
            if HGC.genj_pt[nj]<15 : continue
            if ROOT.TMath.Abs(HGC.genj_eta[nj])>4.7 : continue

            #check if not matched to a tau
            closestMatchId=0
            closestDR=1.0
            for ngen in xrange(0,HGC.ngen):
                deta=HGC.genj_eta[nj]-HGC.gen_eta[ngen]
                dphi=ROOT.TVector2.Phi_mpi_pi(HGC.genj_phi[nj]-HGC.gen_phi[ngen])
                dr=ROOT.TMath.Sqrt(deta*deta+dphi*dphi)
                if dr>0.5 : continue
                if dr>closestDR : continue
                closestDR=dr
                closestMatchId=HGC.gen_id[ngen]
            if ROOT.TMath.Abs(closestMatchId)>=11 and ROOT.TMath.Abs(closestMatchId)<=16 : 
                nLepMatches+=1
                continue

            nselJets=len(jetP4)
            jetP4.append([nj,ROOT.TLorentzVector()])
            jetP4[nselJets][1].SetPtEtaPhiE( HGC.genj_pt[nj], HGC.genj_eta[nj], HGC.genj_phi[nj], HGC.genj_en[nj] )
        if len(jetP4)<2 : continue
        histos['cutflow'].Fill(1)
        
        #define the tag jets with highest Mjj
        tagJet1Idx, tagJet2Idx = 0, 1
        for nj in xrange(0,len(jetP4)):
            for mj in xrange(nj+1,len(jetP4)):
                newDijet=jetP4[nj][1]+jetP4[mj][1]
                dijet=jetP4[tagJet1Idx][1]+jetP4[tagJet2Idx][1]                
                if newDijet.M()<dijet.M() : 
                    continue
                tagJet1Idx=nj
                tagJet2Idx=mj

        tagJet1,tagJet2=jetP4[tagJet1Idx][1], jetP4[tagJet2Idx][1]
        dijet=tagJet1+tagJet2
        deta=ROOT.TMath.Abs(tagJet1.Eta()-tagJet2.Eta())
        histos['mjj'].Fill(dijet.M())
        histos['detajj'].Fill(deta)
        if dijet.M()<250 or deta<2.5 : continue
        histos['cutflow'].Fill(2)
        
        #requirement for jets in HGC
        njetsInHGC=0
        for jetIdx in [jetP4[tagJet1Idx][0],jetP4[tagJet2Idx][0]]:

            absEta=ROOT.TMath.Abs(HGC.genj_eta[jetIdx])

            if absEta>1.5 and absEta<3 : 

                njetsInHGC+=1
                histos['emFrac'].Fill(HGC.genj_en[jetIdx],HGC.genj_emfrac[jetIdx])
                histos['hadFrac'].Fill(HGC.genj_en[jetIdx],HGC.genj_hadfrac[jetIdx])
                histos['invFrac'].Fill(HGC.genj_en[jetIdx],HGC.genj_invfrac[jetIdx])

                hitsPerLayer={}
                xyCentres={}
                for nlay in xrange(1,55): 
                    hitsPerLayer[nlay]=[]
                    xyCentres[nlay]=[0,0]
                for ihit in xrange(0,HGC.nhits):
                    deta=HGC.genj_eta[jetIdx]-HGC.hit_eta[ihit]
                    dphi=ROOT.TVector2.Phi_mpi_pi(HGC.genj_phi[jetIdx]-HGC.hit_phi[ihit])
                    dr=ROOT.TMath.Sqrt(deta*deta+dphi*dphi)
                    if dr>0.6 : continue
                    hitsPerLayer[ HGC.hit_layer[ihit] ]. append( [dr,ihit] )
                    refRho=ROOT.TMath.Abs(HGC.hit_z[ihit]/ROOT.TMath.SinH(HGC.genj_eta[jetIdx]))
                    xyCentres[ HGC.hit_layer[ihit] ] = [refRho*ROOT.TMath.Cos( HGC.genj_phi[jetIdx] ), refRho*ROOT.TMath.Sin(HGC.genj_phi[jetIdx] ) ]

                #fill variables of interest
                for var in varNames:
                    if var=='genEn'        : varNames[var][0]=HGC.genj_en[jetIdx]
                    elif var=='genEta'     : varNames[var][0]=HGC.genj_eta[jetIdx]
                    elif var=='genPhi'     : varNames[var][0]=HGC.genj_phi[jetIdx]
                    elif var=='genPt'      : varNames[var][0]=HGC.genj_pt[jetIdx]
                    elif var=='genEmFrac'  : varNames[var][0]=HGC.genj_emfrac[jetIdx]
                    elif var=='genHadFrac' : varNames[var][0]=HGC.genj_hadfrac[jetIdx]
                    elif var=='genInvFrac' : varNames[var][0]=HGC.genj_invfrac[jetIdx]
                    else                   : varNames[var][0]=0 

                for lay in hitsPerLayer:
                    hitsPerLayer[lay]=sorted( hitsPerLayer[lay], key=lambda drhit: drhit[0])

                    #sub-detector
                    det='EE'
                    if lay>30 and lay<=42: det='FH'
                    if lay>42: det='BH'

                    hitCtr=0
                    for hit in hitsPerLayer[lay]:
                        hitCtr+=1

                        dx=ROOT.TMath.Abs( xyCentres[lay][0]-HGC.hit_x[hit[1]] )
                        dy=ROOT.TMath.Abs( xyCentres[lay][1]-HGC.hit_y[hit[1]] )

                        #check distance to centre
                        keysForHit=[]
                        if hit[0]<0.1        : keysForHit.append('dr01')
                        if hit[0]<0.2        : keysForHit.append('dr02')
                        if hit[0]<0.3        : keysForHit.append('dr03')
                        if hit[0]<0.4        : keysForHit.append('dr04')
                        if dx<1.0 and dy<1.0 : keysForHit.append('2x2')
                        if dx<1.5 and dy<1.5 : keysForHit.append('3x3')
                        if dx<4.5 and dy<4.5 : keysForHit.append('5x5')

                        #variables of interest for the hit
                        eta           = HGC.hit_eta[ hit[1] ]
                        phi           = HGC.hit_phi[ hit[1] ]
                        edep          = HGC.hit_edep[ hit[1] ]
                        edep_em       = scaleToEM(edep, lay, eta)
                        edep_em_wgt   = edep_em*HGC.hit_wgt[ hit[1] ]
                        edep_em_wgt_t = edep_em*HGC.hit_wgt_t[ hit[1] ]

                        for integKey in keysForHit:
                            varNames['%s_%s_raw_en'%(det,integKey)][0]       += edep
                            varNames['%s_%s_raw_eta'%(det,integKey)][0]      += edep*eta
                            varNames['%s_%s_raw_phi'%(det,integKey)][0]      += edep*phi
                            varNames['%s_%s_em_en'%(det,integKey)][0]        += edep_em
                            varNames['%s_%s_em_eta'%(det,integKey)][0]       += edep_em*eta
                            varNames['%s_%s_em_phi'%(det,integKey)][0]       += edep_em*phi
                            varNames['%s_%s_em_wgt_en'%(det,integKey)][0]    += edep_em_wgt
                            varNames['%s_%s_em_wgt_eta'%(det,integKey)][0]   += edep_em_wgt*eta
                            varNames['%s_%s_em_wgt_phi'%(det,integKey)][0]   += edep_em_wgt*phi
                            varNames['%s_%s_em_wgt_t_en'%(det,integKey)][0]  += edep_em_wgt_t
                            varNames['%s_%s_em_wgt_t_eta'%(det,integKey)][0] += edep_em_wgt_t*eta
                            varNames['%s_%s_em_wgt_t_phi'%(det,integKey)][0] += edep_em_wgt_t*phi

                ntuple.Fill()

        if njetsInHGC==1 : histos['cutflow'].Fill(3)
        if njetsInHGC==2 : histos['cutflow'].Fill(4)
        histos['etaLead'].Fill( ROOT.TMath.Max( ROOT.TMath.Abs(tagJet1.Eta()), ROOT.TMath.Abs(tagJet2.Eta()) ) )
        histos['etaTrailer'].Fill( ROOT.TMath.Min( ROOT.TMath.Abs(tagJet1.Eta()), ROOT.TMath.Abs(tagJet2.Eta()) ) )
        histos['enLead'].Fill( ROOT.TMath.Max(tagJet1.E(),tagJet2.E()) )
        histos['enTrailer'].Fill( ROOT.TMath.Min(tagJet1.E(),tagJet2.E()) )
        
    #all done here
    fIn.Close()

    fOut.Write()
    fOut.Close()


"""
steer 
"""
def main():
    
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i',      '--in' ,      dest='input',        help='Input file',              default=None)
    parser.add_option('-o',      '--out' ,     dest='output',       help='Output directory',        default='./')
    (opt, args) = parser.parse_args()
                                       
    #check inputs
    if opt.input is None:
        parser.print_help()
        sys.exit(1)

    #basic ROOT customization
    customROOTstyle()
    ROOT.gROOT.SetBatch(False)
    #ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    runGenLevelAnalysis(opt)
    
if __name__ == "__main__":
    sys.exit(main())
