import ROOT

#get other things from ROOT
from ROOT import TH1F, TH2F, Math

from math import tanh, floor, sqrt, cos

import numpy as np
from numpy import random as rng

#load library with dictionaries for objects in tree
ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.AutoLibraryLoader.enable()

#open ROOT file
fIn=ROOT.TFile.Open('HGCROIAnalyzer_pion10.root')

#read tree from file
tree=fIn.Get('analysis/HGC')
print 'Preparing to analyse %d events'%tree.GetEntriesFast()

cm_per_ns = 29.9792458

# get the histogram of a single shower

def integral(bins_and_values):     
     result = 0.0
     for bin in bins_and_values:
          result += bin[1]
     return result
          

def mean(bins_and_values):
     result = 0.0     
     for bin in bins_and_values:
          result += bin[0]*bin[1]
     norm = integral(bins_and_values)
     if norm <= 0.0: return -1
     return result/norm

eff_sigma_bound = 0.1585
def get_eff_rms(histogram):
     idx_start = 0
     idx_stop  = 0
     tot_integral = integral(histogram)
     if( tot_integral <= 0.0 ): return -1, -1, -1
     while integral(histogram[:idx_start+1])/tot_integral < eff_sigma_bound: idx_start += 1     
     while (1 - integral(histogram[:idx_stop+1])/tot_integral) > eff_sigma_bound: idx_stop += 1
     val_start = histogram[idx_start][0]
     val_stop = histogram[idx_stop][0]
     (val_stop - val_start)/2.0
     return (val_stop - val_start)/2.0, idx_start, idx_stop


cluster_times = {}
cluster_times['all'] = TH1F('cluster_times',';cluster #Deltat (ns);',600,-0.1,0.5)
cluster_times['3x3'] = TH1F('cluster_times_3x3',';cluster #Deltat (ns);',600,-0.1,0.5)
cluster_times['5x5'] = TH1F('cluster_times_5x5',';cluster #Deltat (ns);',600,-0.1,0.5)
cluster_times['7x7'] = TH1F('cluster_times_7x7',';cluster #Deltat (ns);',600,-0.1,0.5)

cluster_times_20up = {}
cluster_times_20up['all'] = TH1F('cluster_times_20up',';cluster #Deltat (ns);',600,-0.1,0.5)
cluster_times_20up['3x3'] = TH1F('cluster_times_20up_3x3',';cluster #Deltat (ns);',600,-0.1,0.5)
cluster_times_20up['5x5'] = TH1F('cluster_times_20up_5x5',';cluster #Deltat (ns);',600,-0.1,0.5)
cluster_times_20up['7x7'] = TH1F('cluster_times_20up_7x7',';cluster #Deltat (ns);',600,-0.1,0.5)

cluster_times_20down = {}
cluster_times_20down['all'] = TH1F('cluster_times_20down',';cluster #Deltat (ns);',600,-0.1,0.5)
cluster_times_20down['3x3'] = TH1F('cluster_times_20down_3x3',';cluster #Deltat (ns);',600,-0.1,0.5)
cluster_times_20down['5x5'] = TH1F('cluster_times_20down_5x5',';cluster #Deltat (ns);',600,-0.1,0.5)
cluster_times_20down['7x7'] = TH1F('cluster_times_20down_7x7',';cluster #Deltat (ns);',600,-0.1,0.5)

cluster_hitmult = {}
cluster_hitmult['all'] = TH1F('cluster_hitmult',';cluster used-hit multiplicity;',200,0,200)
cluster_hitmult['3x3'] = TH1F('cluster_hitmult_3x3',';cluster used-hit multiplicity;',200,0,200)
cluster_hitmult['5x5'] = TH1F('cluster_hitmult_5x5',';cluster used-hit multiplicity;',200,0,200)
cluster_hitmult['7x7'] = TH1F('cluster_hitmult_7x7',';cluster used-hit multiplicity;',200,0,200)

hit_times = {}
hit_times['all'] = TH1F('hit_times',';hit #Deltat (ns);',1000,-0.5,0.5)
hit_times['3x3'] = TH1F('hit_times_3x3',';hit #Deltat (ns);',1000,-0.5,0.5)
hit_times['5x5'] = TH1F('hit_times_5x5',';hit #Deltat (ns);',1000,-0.5,0.5)
hit_times['7x7'] = TH1F('hit_times_7x7',';hit #Deltat (ns);',1000,-0.5,0.5)

hit_times_inlayers = {}
hit_times_inlayers['all'] = TH2F('hit_times_inlayers',';layer;hit #Deltat (ns)',600,317,617,1000,-0.5,0.5)
hit_times_inlayers['3x3'] = TH2F('hit_times_3x3_inlayers',';layer;hit #Deltat (ns)',600,317,617,1000,-0.5,0.5)
hit_times_inlayers['5x5'] = TH2F('hit_times_5x5_inlayers',';layer;hit #Deltat (ns)',600,317,617,1000,-0.5,0.5)
hit_times_inlayers['7x7'] = TH2F('hit_times_7x7_inlayers',';layer;hit #Deltat (ns)',600,317,617,1000,-0.5,0.5)

hit_times_fromaxis = {}
hit_times_fromaxis['all'] = TH2F('hit_times_fromaxis',';transverse distance from shower axis (cm);hit #Deltat (ns)',600,0,6,1000,-0.5,0.5)
hit_times_fromaxis['3x3'] = TH2F('hit_times_3x3_fromaxis',';transverse distance from shower axis (cm);hit #Deltat (ns)',600,0,6,1000,-0.5,0.5)
hit_times_fromaxis['5x5'] = TH2F('hit_times_5x5_fromaxis',';transverse distance from shower axis (cm);hit #Deltat (ns)',600,0,6,1000,-0.5,0.5)
hit_times_fromaxis['7x7'] = TH2F('hit_times_7x7_fromaxis',';transverse distance from shower axis (cm);hit #Deltat (ns)',600,0,6,1000,-0.5,0.5)

hit_times_inlayers_20up = {}
hit_times_inlayers_20up['all'] = TH2F('hit_times_inlayers_20up',';layer;hit #Deltat (ns)',600,317,617,1000,-0.5,0.5)
hit_times_inlayers_20up['3x3'] = TH2F('hit_times_3x3_inlayers_20up',';layer;hit #Deltat (ns)',600,317,617,1000,-0.5,0.5)
hit_times_inlayers_20up['5x5'] = TH2F('hit_times_5x5_inlayers_20up',';layer;hit #Deltat (ns)',600,317,617,1000,-0.5,0.5)
hit_times_inlayers_20up['7x7'] = TH2F('hit_times_7x7_inlayers_20up',';layer;hit #Deltat (ns)',600,317,617,1000,-0.5,0.5)

hit_times_inlayers_20down = {}
hit_times_inlayers_20down['all'] = TH2F('hit_times_inlayers_20down',';layer;hit #Deltat (ns)',600,317,617,1000,-0.5,0.5)
hit_times_inlayers_20down['3x3'] = TH2F('hit_times_3x3_inlayers_20down',';layer;hit #Deltat (ns)',600,317,617,1000,-0.5,0.5)
hit_times_inlayers_20down['5x5'] = TH2F('hit_times_5x5_inlayers_20down',';layer;hit #Deltat (ns)',600,317,617,1000,-0.5,0.5)
hit_times_inlayers_20down['7x7'] = TH2F('hit_times_7x7_inlayers_20down',';layer;hit #Deltat (ns)',600,317,617,1000,-0.5,0.5)

variations = ['all','3x3','5x5','7x7']

##RooFit stuff for dealing with the admixture of expected resolutions
exp_reso = ROOT.RooRealVar('exp_reso','Expected Resolution (ns)',0.0,1.0)
cluster_time = ROOT.RooRealVar('clus_time','Vertex Corrected Cluster Arrival Time (ns)',-0.2,0.7)
cluster_data = ROOT.RooDataSet('cluster_data','Cluster Time Data',ROOT.RooArgSet(cluster_time,exp_reso))

#PDF definition
zero = ROOT.RooConstVar('zero','zero',0.0)
mean = ROOT.RooRealVar('mean','Mean Bias',0.0,-5,5)
sigma_scal = ROOT.RooRealVar('sigma_scal','Scaling of Expected Sigma',1.0,1e-6,5)
scaled_reso = ROOT.RooLinearVar('scaled_reso','Resolution',exp_reso,sigma_scal,zero)
gaus = ROOT.RooGaussian('gauss','gaussian resolution',cluster_time,mean,scaled_reso)

total = ROOT.RooProdPdf('total','total',
                        ROOT.RooFit.Conditional(ROOT.RooArgSet(gaus),ROOT.RooArgSet(cluster_time)))

reso_smear = 0.001
for i in xrange(0,tree.GetEntriesFast()):
     tree.GetEntry(i)
     
     clusters_in_rois = {}
     for icl,cluster in enumerate(tree.Clusters):
          if cluster.roiidx_ not in clusters_in_rois.keys():
               clusters_in_rois[cluster.roiidx_] = []
          clusters_in_rois[cluster.roiidx_].append(icl)
     
     #sort the clusters in ROIs by energy
     for roi in clusters_in_rois.keys():
          clusters_in_rois[roi].sort(key=lambda icl: tree.Clusters[icl].en_)
          clusters_in_rois[roi].reverse()
     
     for icl,cluster in enumerate(tree.Clusters):
          interaction_z = abs(tree.ROIs[cluster.roiidx_].stablez_)
          if( icl != clusters_in_rois[cluster.roiidx_][0] ): continue
          if( interaction_z < 317.0 ): continue
          tdc_hits = {} 
          tdc_hits_axis = {}
          for var in variations:
               tdc_hits[var] = []
               tdc_hits_axis[var] = []
          for hit in tree.RecHits:
               if hit.clustId_ != icl: continue
               if hit.t_<0 : continue
               
               # smear the time and re-digitize it
               time = (hit.t_+0.0025)/0.005
               time = floor(rng.normal(time,1e-9 + int(reso_smear/0.005)))
               time = 0.005*time
               
               #correct the time to the gen vertex
               hit_pos = Math.XYZVector(hit.x_,hit.y_,hit.z_)
               axis_vec = Math.XYZVector(cluster.axis_x_,cluster.axis_y_,cluster.axis_z_)
               time_corr = tree.GenVertex.z()*tanh(hit_pos.eta())/cm_per_ns - 1.0
               time = time + time_corr
               time_check = hit.t_ + time_corr
               
               # cos = adj / hyp , hyp = adj / cos = adj / tanh(eta)
               for var in variations:
                    if var == 'all':
                         keep = True
                    else:
                         keep = getattr(hit,'isIn%s_'%var)
                    if keep: 
                         hit_times[var].Fill(time)
                         hit_times_inlayers[var].Fill(abs(hit.z_),time)
                         if( abs(hit.z_) > 322 and time < 0.3 ): 
                              tdc_hits[var].append((hit.en_,time,abs(hit.z_)))
                              cl_axis   = np.array([cluster.axis_x_,cluster.axis_y_,cluster.axis_z_])
                              cl_center = np.array([cluster.center_x_,cluster.center_y_,cluster.center_z_])
                              cl_at_hit = (hit.z_ - cluster.center_z_)/cos(axis_vec.theta())*cl_axis + cl_center #
                              #print hit.x_, hit.y_, hit.z_, " ", cl_at_hit[0], cl_at_hit[1], cl_at_hit[2]
                              delta_xy = sqrt( (cl_at_hit[0] - hit.x_)**2 + (cl_at_hit[1] - hit.y_)**2 )
                              tdc_hits_axis[var].append( (hit.en_,time,hit.x_,hit.y_,hit.z_,delta_xy) )

          for var in variations:
               tdc_hits[var].sort(key=lambda hit: hit[1])
          
          #process the time-sorted view of the shower
          for var in variations:
               hits = tdc_hits[var]
               nhits = len(hits)
               sum_times = sum([hit[1] for hit in hits])
               if nhits > 0:
                    cluster_hitmult[var].Fill(nhits)
                    cluster_times[var].Fill(sum_times/nhits)   
                    
                    # fill fitting data
                    if var == '5x5':
                         pass
       
                    if len(hits) >= 20:
                         cluster_times_20up[var].Fill(sum_times/nhits)
                         for hit in hits:
                              hit_times_inlayers_20up[var].Fill(hit[2],hit[1])
                    elif len(hits) >= 10:
                         cluster_times_20down[var].Fill(sum_times/nhits)
                         for hit in hits:
                              hit_times_inlayers_20down[var].Fill(hit[2],hit[1])
     
          #process the axis-distance sorted view of the shower
          for var in variations:               
               tdc_hits_axis[var].sort(key=lambda hit: hit[5])
               #print "got %i tdc hits"%(len(tdc_hits_axis[var]))
               usable_hits = []
               for hit in tdc_hits_axis[var]:
                    hit_times_fromaxis[var].Fill(hit[5],hit[1])
                    if hit[5] < 0.8:
                         usable_hits.append(hit)
                             
               if var == 'all':
                    nhits = len(usable_hits)
                    sum_times = sum([hit[1] for hit in usable_hits])
                    if nhits > 0 :
                         #print 'total tdc hits = %d , total usable hits = %d'%(len(tdc_hits_axis[var]),nhits)
                         exp_reso.setVal(reso_smear/sqrt(nhits))
                         cluster_time.setVal(sum_times/nhits)
                         cluster_data.add(ROOT.RooArgSet(cluster_time,exp_reso))
                    #print hit[0], hit[1], 'pos info: ', hit[2], hit[3], hit[4],
                    #print cluster.center_x_, cluster.center_y_, cluster.center_z_ ,hit[5]
                    



## now to do some fitting
total.fitTo(cluster_data,ROOT.RooFit.ConditionalObservables(ROOT.RooArgSet(exp_reso)))

frame = cluster_time.frame(-0.1,0.5)
print frame.GetName()

cluster_data.plotOn(frame)
total.plotOn(frame,ROOT.RooFit.ProjWData(cluster_data))


#make it pretty
from ROOT import TPaveText, TCanvas
fitcanv = TCanvas('boop','',600,600)
fitcanv.cd()
frame.SetTitle('')
frame.GetYaxis().SetTitleOffset(1.40)
txt = TPaveText(0.6264,0.6525,0.8261,0.8538,'NDC')
txt.AddText('mean = %.3g'%mean.getVal())
txt.AddText('scaling = %.3g'%sigma_scal.getVal())
txt.SetFillColor(ROOT.kWhite)
frame.Draw()
txt.Draw()

errcanv = TCanvas('beep','',600,600)
errcanv.cd()
err_frame = exp_reso.frame(0,0.200)
err_frame.SetTitle('')
err_frame.GetYaxis().SetTitleOffset(1.40)
cluster_data.plotOn(err_frame)
err_frame.Draw()


"""
          unique_times = {}
          for hit in tdc_hits:
               if hit[1] not in unique_times.keys():
                    unique_times[hit[1]] = 0
               unique_times[hit[1]] += 1
          #for time in sorted(unique_times.keys()):
          #     print icl, time, unique_times[time]
          histogram = [(key,unique_times[key]) for key in sorted(unique_times.keys())]
          eff_rms, bin_start, bin_stop = get_eff_rms(histogram)
          
          print 'Input histo: ', eff_rms
          print 'Prune histo: ', get_eff_rms(histogram[bin_start:bin_stop+1])[0]

          nhits_used = integral(histogram[bin_start:bin_stop+1])
          mean_clus_time = mean(histogram[bin_start:bin_stop+1])

          if nhits_used > 1.0:
               cluster_times.Fill(mean_clus_time)
"""

