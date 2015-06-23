#! /usr/bin/env python

import ROOT

from ROOT import TFile, TH1F,TH2F,TH2I, TCanvas, TLegend,SetOwnership, TGraph, TF1
from math import tanh

import collections
import functools

class memoized(object):
   '''Decorator. Caches a function's return value each time it is called.
   If called later with the same arguments, the cached value is returned
   (not reevaluated).
   '''
   def __init__(self, func):
      self.func = func
      self.cache = {}
   def __call__(self, *args):
      if not isinstance(args, collections.Hashable):
         # uncacheable. a list, for instance.
         # better to not cache than blow up.
         return self.func(*args)
      if args in self.cache:
         return self.cache[args]
      else:
         value = self.func(*args)
         self.cache[args] = value
         return value
   def __repr__(self):
      '''Return the function's docstring.'''
      return self.func.__doc__
   def __get__(self, obj, objtype):
      '''Support instance methods.'''
      return functools.partial(self.__call__, obj)

@memoized
def fast_tanh(x):
    return tanh(x)

energies = {}
resolutions = [20,50,80,100,150,200]
upper_tolerance = 0.1

#infile = TFile.Open('/tmp/lgray/step3_perfect.root_SimHits_0.root')

infile = TFile.Open('/tmp/lgray/Single22_CMSSW_6_2_0_SLHC25_patch4_RECO-PU0__SimHits_0.root')

tree = infile.Get('analysis').Get('HGC')

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1)

eff_sigma_bound = 0.1585
def gef_eff_rms(histo):
    bin_start = 0
    bin_stop  = 0
    tot_integral = histo.Integral()
    while histo.Integral(0,bin_start)/tot_integral < eff_sigma_bound: bin_start += 1

    while (1 - histo.Integral(0,bin_stop)/tot_integral) > eff_sigma_bound: bin_stop += 1
    
    val_start = histo.GetBinCenter(bin_start)
    val_stop = histo.GetBinCenter(bin_stop)
    print histo.GetName(), (val_stop - val_start)/2.0
    return (val_stop - val_start)/2.0

def vertex_corrected_time(time,vtxZ,genEta):
    return time + vtxZ*0.0333564095*tanh(genEta) - 1

def make_histos(histos,energy):
    harr = []
    harr.append(TH2I('ntdchits_vs_eta_%d'%energy,
                     ';|#eta|;TDC Hit Multiplicity',
                     10,1.5,3.0,
                     500,0,500))
    # per layer times
    harr.append(TH1F('vtx_corr_layer_time_0ps_%d'%energy,
                     ';Per-Layer Vertex-Corrected Arrival Time (ns)',
                     800,-2,6))
    for reso in resolutions:
        harr.append(TH1F('vtx_corr_layer_time_%dps_%d'%(reso,energy),
                     ';Per-Layer Vertex-Corrected Arrival Time (ns)',
                     800,-2,6))
    # cluster times
    harr.append(TH1F('vtx_corr_cluster_time_0ps_%d'%energy,
                     ';Cluster Vertex-Corrected Arrival Time (ns)',
                     40000,-2,6))
    for reso in resolutions:
        harr.append(TH1F('vtx_corr_cluster_time_%dps_%d'%(reso,energy),
                     ';Cluster Vertex-Corrected Arrival Time (ns)',
                     40000,-2,6))

    print harr
    histos[energy] = harr
    

def fill_histos(tree,histos,energy):
    cut_string_cluster = 'emeanTime_rec > 0 && hasInteractionBeforeHGC == 0 &&'\
        ' abs(genEn - %d) < 0.1'%(energy)
    cut_string_cluster_hits = 'emeanTime_rec > 0 && hasInteractionBeforeHGC == 0 &&'\
        ' abs(genEn - %d) < 0.1'%(energy)
    cut_string_layer = 'emeanTimeLayer_rec > 0 && hasInteractionBeforeHGC == 0 &&'\
        ' abs(genEn - %d) < 0.1'%(energy)

    print cut_string_layer

    vertex_corr = '(%s - 1) >> %s'
    
    tree.Draw('Sum$(nhitstdc_rec) : abs(genEta) >> %s'%(histos[0].GetName()),
              cut_string_layer,'goff')

    #nTDCHits = 0
    #for nhits_layer in tree.nhitstdc_rec:
    #    nTDCHits += nhits_layer
    #histos[0].Fill(abs(tree.genEta),nTDCHits)
    #time distributions

    tree.Draw(vertex_corr%('emeanTimeLayer_rec',histos[1].GetName()),
              cut_string_layer,'goff')

    for k,reso in enumerate(resolutions):
        tree.Draw(vertex_corr%('emeanTimeLayer%d_rec'%reso,histos[2+k].GetName()),
                  cut_string_layer,'goff')

    tree.Draw(vertex_corr%('emeanTime_rec',histos[8].GetName()),
              cut_string_cluster,'goff')

    for k,reso in enumerate(resolutions):
       cut_string_cluster_reso = 'emeanTime%d_rec > 0 && hasInteractionBeforeHGC == 0 &&  abs(genEn - %d) < 0.1'%(reso,energy)
       
       tree.Draw(vertex_corr%('emeanTime%d_rec'%reso,histos[9+k].GetName()),
                 cut_string_cluster_reso,'goff')

    #for time in tree.emeanTimeLayer_rec:
    #    if time < -0.0 : continue
    #    histos[1].Fill(vertex_corrected_time(time,tree.genVertexZ,tree.genEta))
    #for k,reso in enumerate(resolutions):
    #    times = getattr(tree,'emeanTimeLayer%d_rec'%reso)
    #    for time in times:
    #        if time < -0.0 : continue
    #        histos[2+k].Fill(vertex_corrected_time(time,tree.genVertexZ,tree.genEta))
               

def make_occupancy_plot(histos):
    occ = TCanvas('occ','',600,600)
    occ.SetLogy()
    occ.cd()
    legd = TLegend(0.508,0.7185,0.9631,0.9685,'','brNDC')
    legd.SetNColumns(3)
    legd.SetFillColor(ROOT.kWhite)
    SetOwnership( legd, 0 )
    for i,energy in enumerate(reversed(sorted(histos.keys()))):
        print i,energy
        tdc_occ_histo = histos[energy][0]
        tdc_occ_histo.SetLineColor(i+1)
        legd.AddEntry(tdc_occ_histo,'%d GeV'%energy)
        if i == 0 :
            prof = tdc_occ_histo.ProfileX()
            prof.GetYaxis().SetRangeUser(5e-1,1100)
            prof.GetYaxis().SetTitle('<TDC Hit Multiplicity>')
            prof.GetYaxis().SetTitleOffset(1.35)
            prof.Draw()
        else:
            tdc_occ_histo.ProfileX().Draw("same")    

    legd.Draw()
    occ.Update()
    occ.Print('tdc_occupancy_vs_eta.pdf')
    del occ

def make_reso_plot(histos):
    reso_layer = TCanvas('reso_layer','',600,600)
    reso_cluster = TCanvas('reso_cluster','',600,600)
    reso_80ps_energy = TCanvas('reso_80ps_energy','',600,600)

    reso_layer.SetLogx()
    reso_cluster.SetLogx()

    reso_layer.SetGridy()
    reso_cluster.SetGridy()

    reso_layer.cd()
    legd = TLegend(0.508,0.7185,0.9631,0.9685,'','brNDC')
    legd.SetNColumns(2)
    legd.SetFillColor(ROOT.kWhite)
    SetOwnership( legd, 0 )
    graphs_layer = []
    graphs_cluster = []
    
    legd2 = TLegend(0.508,0.7185,0.9631,0.9685,'','brNDC')
    legd2.SetNColumns(3)
    legd2.SetFillColor(ROOT.kWhite)
    SetOwnership( legd2, 0 )

    mygaus = TF1('mygaus','gaus')

    for ires in reversed(range(2,8)): #loop over the resolutions
        graphs_layer.append(TGraph(len(histos.keys())))
        graphs_cluster.append(TGraph(len(histos.keys())))
        layer_graph = graphs_layer[-1]
        layer_graph.SetTitle('')
        cluster_graph = graphs_cluster[-1]
        cluster_graph.SetTitle('')
        layer_graph.SetName('layer_reso_vs_energy_%d'%ires)
        cluster_graph.SetName('cluster_reso_vs_energy_%d'%ires)
        
        #cluster_graph.SetMarkerStyle(2)
        if( ires +1 != 5 ):
           layer_graph.SetLineColor(ires+1)
           cluster_graph.SetLineColor(ires+1)
        else:
           layer_graph.SetLineColor(9)
           cluster_graph.SetLineColor(9)

        legd.AddEntry(layer_graph,'Cell #Delta_{t} = %dps'%resolutions[ires-2],'lp')

        for i,energy in enumerate(sorted(histos.keys())):        
            histos[energy][ires].Fit(mygaus,'nodraw')   

            
            
            sigma = mygaus.GetParameter(2)
            
            layer_graph.SetPoint(i,energy,sigma)

            histos[energy][ires+7].Fit(mygaus,'nodraw')
            sigma = gef_eff_rms(histos[energy][ires+7]) #mygaus.GetParameter(2)
            cluster_graph.SetPoint(i,energy,sigma)

        
        if ires == 7: 
            reso_layer.cd()            
            layer_graph.Draw('APL')
            layer_graph.GetYaxis().SetRangeUser(0.0,0.2)
            layer_graph.GetYaxis().SetTitle("Layer #Delta_{t} (ns)")
            layer_graph.GetYaxis().SetTitleOffset(1.5)
            layer_graph.GetXaxis().SetTitle("Pion Energy (GeV)")            
            layer_graph.GetXaxis().SetTitleOffset(1.2)
            reso_cluster.cd()
            cluster_graph.Draw('APL')
            cluster_graph.GetYaxis().SetRangeUser(0.0,0.2)
            cluster_graph.GetYaxis().SetTitle("Cluster #Delta_{t} (ns)")
            cluster_graph.GetYaxis().SetTitleOffset(1.5)
            cluster_graph.GetXaxis().SetTitle("Pion Energy (GeV)")
            cluster_graph.GetXaxis().SetTitleOffset(1.2)
        else:
            reso_layer.cd()
            layer_graph.Draw('same')
            reso_cluster.cd()
            cluster_graph.Draw('same')

        if resolutions[ires-2] == 80:
           reso_80ps_energy.cd()
           for i,energy in enumerate(sorted(histos.keys())):
              
              hist = histos[energy][ires+7]
              rebinned = hist.Rebin(100,'%s_rebin'%hist.GetName())              
              rebinned.SetLineColor(i)
              legd2.AddEntry(rebinned,'E = %d GeV'%energy,'lp')
              if i == 0:
                 #rebinned.GetYaxis().SetRangeUser(1e-3,0.35)
                 rebinned.GetYaxis().SetTitle('a.u.')
                 rebinned.DrawNormalized()
                 
              else:
                 rebinned.DrawNormalized('same')
    

    reso_layer.cd()
    legd.Draw()
    reso_layer.Update()
    reso_layer.Print('resolayer_vs_pt_vs_cellreso.pdf')

    reso_cluster.cd()
    legd.Draw()
    reso_cluster.Update()
    reso_cluster.Print('resocluster_vs_pt_vs_cellreso.pdf')
    #del reso

    reso_80ps_energy.cd()
    legd2.Draw()
    reso_80ps_energy.Update()
    reso_80ps_energy.Print('80ps_resolution_overlay.pdf')

    return reso_layer, reso_cluster, reso_80ps_energy, graphs_layer, graphs_cluster

for event in xrange(tree.GetEntries()):
    tree.GetEntry(event)
    
    genEn = tree.genEn
    energy_bin = int(genEn + upper_tolerance)
    if energy_bin not in energies.keys():
        make_histos(energies,energy_bin)

print energies.keys()


for energy_bin in energies.keys():
    fill_histos(tree,energies[energy_bin],energy_bin)

make_occupancy_plot(energies)

canvas1, canvas2, canvas3, layers, clusters = make_reso_plot(energies)
