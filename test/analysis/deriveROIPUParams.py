#!/usr/bin/python

import ROOT
from array import array

"""
reads of the distributions, projects, computes median and sigma
"""
def deriveROIPUParams(url, dist='csi', ndRbins=12, nCsiBins=56):

	print 'Deriving ROI PU params for %s'%dist

	#open file
	fIn=ROOT.TFile.Open(url)
	csidr      = fIn.Get('analysis/%sdr'%dist)
	csiMin     = csidr.GetYaxis().GetBinLowEdge(1)
	csiMax     = csidr.GetYaxis().GetBinUpEdge(nCsiBins)
	drMin      = csidr.GetXaxis().GetBinLowEdge(1)
	drMax      = csidr.GetXaxis().GetBinUpEdge(ndRbins)
	regs       = fIn.Get('analysis/regs')
	nLayerBins = regs.GetXaxis().GetNbins()
	nEtaBins   = regs.GetYaxis().GetNbins()

	#projections per eta range
	csiH=[]
	for ieta in xrange(0,nEtaBins):
		csiH.append( ROOT.TH1F('csi_%d'%ieta,'|#eta|=%3.2f;#xi;PDF',nCsiBins,csiMin,csiMax) )
		csiH[ieta].SetMarkerStyle(ieta+20)
		csiH[ieta].Sumw2()
		csiH[ieta].SetDirectory(0)

	#final output histograms characterizing the median and left width as function of dR in the ROI and eta of the ROI
	medianH=ROOT.TH2F('medianPU_%s'%dist,';Layer;#Delta R',nLayerBins,1,nLayerBins+1,ndRbins*nEtaBins,  drMin, drMax+(drMax-drMin)*(nEtaBins-1))
	medianH.Sumw2()
	medianH.SetDirectory(0)
	widthH=medianH.Clone('widthPU_%s'%dist)
	widthH.SetDirectory(0)
	sigma1H=medianH.Clone('sigma1PU_%s'%dist)
	sigma1H.SetDirectory(0)
	sigma2H=medianH.Clone('sigma2PU_%s'%dist)
	sigma2H.SetDirectory(0)

	#quantiles to determine per slice
	probSum   = array('d',  [0.5,0.683,0.954])
	quantiles = array('d', [0.0]*len(probSum))

	#loop over slices in dR,layer,eta
	for xbin in xrange(1,csidr.GetXaxis().GetNbins()+1):

		#determine dR for this slice
		normXbin=(xbin-1)%ndRbins+1
		dr=csidr.GetXaxis().GetBinCenter(normXbin)

		#determine layer for this slice
		ilay=(xbin-1)/ndRbins

		#project csi PDFs for different etas
		proj=csidr.ProjectionY("p_csi",xbin,xbin)
		for h in csiH : h.Reset('ICE')

		#fill csi PDFs per eta
		for ybin in xrange(1,csidr.GetYaxis().GetNbins()+1):

			#determine eta for this bin
			ieta=(ybin-1)/nCsiBins

			#determine true csi for this bin 
			normYbin=(ybin-1)%nCsiBins+1
			csi=csidr.GetYaxis().GetBinCenter(normYbin)

			#fill eta PDF
			csiH[ieta].SetBinContent(normYbin, proj.GetBinContent(ybin) )
			csiH[ieta].SetBinError(normYbin, proj.GetBinError(ybin) )

		#normalize PDFs and compute quantiles
		for ieta in xrange(0,nEtaBins):
			totalCts=csiH[ieta].Integral()
			if totalCts<=0: continue
			csiH[ieta].Scale(1./totalCts)
			csiH[ieta].GetQuantiles(len(probSum), quantiles, probSum)
			medianH.SetBinContent(ilay+1,normXbin+ieta*ndRbins,quantiles[0])
			sigma1H.SetBinContent(ilay+1,normXbin+ieta*ndRbins,quantiles[1])
			sigma2H.SetBinContent(ilay+1,normXbin+ieta*ndRbins,quantiles[2])
						
			leftRMS=0
			npForLeftRMS=0
			for csiXbin in xrange(1,csiH[ieta].GetXaxis().GetNbins()):
				icsi=csiH[ieta].GetXaxis().GetBinLowEdge(csiXbin)
				ncounts=csiH[ieta].GetBinContent(csiXbin)
				if icsi<quantiles[0]:
					npForLeftRMS+=ncounts
					leftRMS+=ncounts*ROOT.TMath.Power(icsi-quantiles[0],2)
			if npForLeftRMS>0:
				leftRMS=ROOT.TMath.Sqrt(leftRMS/npForLeftRMS)
			widthH.SetBinContent(ilay+1,normXbin+ieta*ndRbins,leftRMS)
	fIn.Close()

	return medianH, widthH, sigma1H, sigma2H

"""
mani function
"""
def main() :

	ROOT.gStyle.SetOptStat(0)
	ROOT.gStyle.SetOptTitle(0)

	medianCsiH, widthCsiH, sigma1CsiH, sigma2CsiH   = deriveROIPUParams( url='muontag.root', dist='csi' )
	medianCsiTH, widthCsiTH,sigma1CsiTH, sigma2CsiTH = deriveROIPUParams( url='muontag.root', dist='csit' )

	#save output
	fOut=ROOT.TFile.Open('ROIPUparams.root','RECREATE')
	medianCsiH.Write()
	widthCsiH.Write()
	medianCsiTH.Write()
	widthCsiTH.Write()
	sigma1CsiH.Write()
	sigma2CsiH.Write()
	sigma1CsiTH.Write()
	sigma2CsiTH.Write()
	fOut.Close()


if __name__ == "__main__":
    main()
