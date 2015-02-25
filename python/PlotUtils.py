import ROOT
import sys

"""
simple progress bar
"""
def drawProgressBar(percent, barLen = 100):
    sys.stdout.write("\r")
    progress = ""
    for i in range(barLen):
        if i < int(barLen * percent):
            progress += "="
        else:
            progress += " "
    sys.stdout.write("[ %s ] %.1f%%" % (progress, percent * 100))
    sys.stdout.flush()

"""

"""
class MyPaveText(ROOT.TPaveText):
    def __init__(self,title,x1=0.12,y1=0.95,x2=0.6,y2=0.99):
        ROOT.TPaveText.__init__(self,x1,y1,x2,y2,'brNDC')
        self.SetFillStyle(0)
        self.SetBorderSize(0)
        self.SetTextAlign(12)
        self.SetTextFont(42)
        self.SetTextSize(0.04)
        for t in title.split('\\') : self.AddText(t)
        self.Draw()


"""
computes a chi2 between two TGraphERrors
"""
def getChi2(gr1,gr2):

    chi2,ndf=0,0
    x,y1,y2 = ROOT.Double(0), ROOT.Double(0),ROOT.Double(0)
    for np in xrange(0,gr1.GetN()) :

        gr1.GetPoint(np,x,y1)
        y1_err=gr1.GetErrorY(np)

        gr2.GetPoint(np,x,y2)
        y2_err=gr2.GetErrorY(np)

        den=(y1_err*y1_err+y2_err*y2_err)
        if den==0: continue

        ndf+=1
        chi2+=ROOT.TMath.Power(y1-y2,2)/den
        
    return chi2,ndf

"""
Computes the ratio of the two graphs
"""
def getRatio(gr1,gr2):
	grratio=ROOT.TGraphErrors()
	grratio.SetName('%s_ratio_%s'%(gr1.GetName(),gr2.GetName()))
        grratio.SetTitle(gr2.GetTitle())
	grratio.SetLineColor(gr2.GetLineColor())
	grratio.SetMarkerColor(gr2.GetMarkerColor())
	grratio.SetMarkerStyle(gr2.GetMarkerStyle())
	grratio.SetLineStyle(gr2.GetLineStyle())
        x,y1,y2 = ROOT.Double(0), ROOT.Double(0),ROOT.Double(0)
	for np in xrange(0,gr1.GetN()):
             
            gr1.GetPoint(np,x,y1)
            y1_err=gr1.GetErrorY(np)
            
            gr2.GetPoint(np,x,y2)
            y2_err=gr2.GetErrorY(np)

            if y1==0 : continue

            nnp=grratio.GetN()
            grratio.SetPoint(nnp,x,y2/y1)
            grratio.SetPointError(nnp,0,ROOT.TMath.Sqrt(ROOT.TMath.Power(y2*y1_err,2)+ROOT.TMath.Power(y1*y2_err,2))/(y1*y1))

        return grratio



"""
Style options mostly from CMS's tdrStyle.C
"""
def customROOTstyle() :
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(False)
    ROOT.gStyle.SetPadTopMargin(0.05);
    ROOT.gStyle.SetPadBottomMargin(0.13);
    ROOT.gStyle.SetPadLeftMargin(0.12);
    ROOT.gStyle.SetPadRightMargin(0.04);
    ROOT.gStyle.SetLabelColor(1, "XYZ");
    ROOT.gStyle.SetLabelFont(42, "XYZ");
    ROOT.gStyle.SetLabelOffset(0.007, "XYZ");
    ROOT.gStyle.SetLabelSize(0.035, "XYZ");
    ROOT.gStyle.SetTitleSize(0.04, "XYZ");
    ROOT.gStyle.SetAxisColor(1, "XYZ");
    ROOT.gStyle.SetStripDecimals(True);
    ROOT.gStyle.SetTickLength(0.03, "XYZ");
    ROOT.gStyle.SetNdivisions(510, "XYZ");
    ROOT.gStyle.SetPadTickX(0);
    ROOT.gStyle.SetPadTickY(0);
    ROOT.gStyle.SetMarkerStyle(20);
    ROOT.gStyle.SetHistLineColor(1);
    ROOT.gStyle.SetHistLineStyle(0);
    ROOT.gStyle.SetHistLineWidth(1);
    ROOT.gStyle.SetFrameBorderMode(0);
    ROOT.gStyle.SetFrameBorderSize(1);
    ROOT.gStyle.SetFrameFillColor(0);
    ROOT.gStyle.SetFrameFillStyle(0);
    ROOT.gStyle.SetFrameLineColor(1);
    ROOT.gStyle.SetFrameLineStyle(1);
    ROOT.gStyle.SetFrameLineWidth(1);

"""
Add overflows to the bins
"""
def fixExtremities(h):
   
    fbin = h.GetBinContent(0) + h.GetBinContent(1);
    fbine = ROOT.TMath.Sqrt(h.GetBinError(0)*h.GetBinError(0)
                            + h.GetBinError(1)*h.GetBinError(1))
    h.SetBinContent(1,fbin)
    h.SetBinError(1,fbine)
    h.SetBinContent(0,0)
    h.SetBinError(0,0)

    nbins = h.GetNbinsX()
    fbin = h.GetBinContent(nbins) + h.GetBinContent(nbins+1)
    fbine = ROOT.TMath.Sqrt(h.GetBinError(nbins)*h.GetBinError(nbins)
                          + h.GetBinError(nbins+1)*h.GetBinError(nbins+1))
    h.SetBinContent(nbins,fbin)
    h.SetBinError(nbins,fbine)
    h.SetBinContent(nbins+1,0)
    h.SetBinError(nbins+1,0)
