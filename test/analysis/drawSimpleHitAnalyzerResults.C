/**
   
 */

#include "TSystem.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TString.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TObjArray.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TF1.h"
#include "TH2.h"

void fixExtremities(TH1* h,bool addOverflow, bool addUnderflow);
void showSummary(TObjArray plots,TString name,TString title,TString outDir);

void fixExtremities(TH1* h,bool addOverflow, bool addUnderflow)
{
  if(h==0) return;

  if(addUnderflow){
      double fbin  = h->GetBinContent(0) + h->GetBinContent(1);
      double fbine = sqrt(h->GetBinError(0)*h->GetBinError(0)
                          + h->GetBinError(1)*h->GetBinError(1));
      h->SetBinContent(1,fbin);
      h->SetBinError(1,fbine);
      h->SetBinContent(0,0);
      h->SetBinError(0,0);
    }
  
  if(addOverflow){  
      int nbins = h->GetNbinsX();
      double fbin  = h->GetBinContent(nbins) + h->GetBinContent(nbins+1);
      double fbine = sqrt(h->GetBinError(nbins)*h->GetBinError(nbins) 
                          + h->GetBinError(nbins+1)*h->GetBinError(nbins+1));
      h->SetBinContent(nbins,fbin);
      h->SetBinError(nbins,fbine);
      h->SetBinContent(nbins+1,0);
      h->SetBinError(nbins+1,0);
    }
}


//
void drawSimpleHitAnalyzerResults(TString inURL="SimpleHitAnalysis.root",TString outDir="~/www/HGCal/SimpleHitAnalysis/",int version=6)
{
  //prepare output
  outDir += "v"; outDir+= version;
  gSystem->Exec("mkdir -p " +outDir);
  gSystem->Exec("cp $CMSSW_BASE/src/UserCode/HGCanalysis/test/analysis/simplehitanalyzer_index.html " + outDir +"/index.html");

  TString dists[]={ "en","time","deta","dphi" };
  TFile *_file0 = TFile::Open(inURL);

  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  Int_t maxLayers(44);
  if(version==5) maxLayers=54;
  std::vector<TGraphErrors *> profiles;
  TF1 *lan=new TF1("lan","[0]*TMath::Landau(x,[1],[2])",0,2000);
  for(Int_t ilay=1; ilay<=maxLayers; ilay++)
    {

      TObjArray plots;       
      for(size_t idist=0; idist<sizeof(dists)/sizeof(TString); idist++)
	{  
	  bool doLandauFit(dists[idist]=="en");
	 
	  TString hname(dists[idist]); hname += "_layer"; hname+=ilay;
	  TH1 *h=(TH1 *)_file0->Get("analysis/"+hname);

	  float profVal(h->GetMean()),profValUnc(h->GetMeanError());
	  if(doLandauFit)
	    {
	      Int_t maxBin=h->GetMaximumBin();
	      Float_t mpvest=h->GetXaxis()->GetBinCenter(maxBin);
	      Float_t fitmin=h->GetXaxis()->GetBinCenter(maxBin-5);
	      Float_t fitmax=h->GetXaxis()->GetBinCenter(maxBin+15);
	      h->Fit(lan,"WQR+","",fitmin,fitmax);
	      lan->SetParameter(1,mpvest);
	      profVal=lan->GetParameter(1);
	      profValUnc=lan->GetParError(1);
	    }
	  
	  if(profiles.size()<=idist)
	    {
	      profiles.push_back(new TGraphErrors);  
	      profiles[idist]->SetName(dists[idist]+"_prof");
	      profiles[idist]->SetTitle( h->GetXaxis()->GetTitle() );
	      profiles[idist]->SetMarkerStyle(20); 
	      profiles[idist]->SetLineWidth(2);
	    }
	  Int_t np=profiles[idist]->GetN();
	  profiles[idist]->SetPoint(np,ilay,profVal);
	  profiles[idist]->SetPointError(np,0,profValUnc);
	  
	  fixExtremities(h,true,true);
	  plots.Add(h);
	}

      TString title("HGC-EE");
      if(ilay>30) title="HGC-FH"; 
      if(ilay>42) title=(version==5 ? "HGC-BH" : "HE rebuild");
      title += ", layer=";
      title += ilay;
      TString name("layersummary_"); name+= ilay;
      showSummary(plots,name,title,outDir);
    }

  
      
  TObjArray plots; 
  for(size_t i=0; i<profiles.size(); i++) plots.Add( profiles[i] );
  TString title("");
  showSummary(plots,"layersummary_prof","HGCal simulation profile",outDir);
}

//
void showSummary(TObjArray plots,TString name,TString title,TString outDir)
{
  Int_t nplots(plots.GetEntriesFast());

  TCanvas *c=new TCanvas("c","c",2000,400);
  c->SetLeftMargin(0);
  c->SetTopMargin(0);
  c->SetRightMargin(0);
  c->SetBottomMargin(0);
  c->Divide(nplots,1);
  for(Int_t iplot=0; iplot<nplots; iplot++)
    {
      TPad *p=(TPad *)c->cd(iplot+1);
      p->SetTopMargin(0.05);
      p->SetRightMargin(0.02);
      p->SetBottomMargin(0.1);
      p->SetLeftMargin(0.12);
      
      TObject *obj=plots.At(iplot);
      TString className=obj->ClassName();
      if(className.Contains("TH1")) 
	{
	  TH1 *h=(TH1 *)obj;
	  h->Draw("hist");
	  h->SetLineWidth(2);
	  h->SetFillStyle(1001);
	  h->SetFillColor(kGray);
	  h->GetYaxis()->SetTitleOffset(1.2);
	}
      else if(className.Contains("TGraph"))
	{
	  TGraphErrors *gr=(TGraphErrors *)obj;
	  gr->Draw("ap");
	  gr->GetXaxis()->SetTitle("Endcap layer number");
	  TString ytitle(gr->GetTitle());
	  if(ytitle.Contains("Energy")) ytitle="MPV from Landau fit";
	  else ytitle="Average " + ytitle;
	  gr->GetYaxis()->SetTitle( ytitle );
	  gr->GetYaxis()->SetTitleOffset(1.2);
	}
      
      if(iplot>0) continue;
      TPaveText *pt=new TPaveText(0.12,0.96,0.6,0.99,"brNDC");
      pt->SetBorderSize(0);
      pt->SetFillStyle(0);
      pt->SetTextFont(42);
      pt->SetTextAlign(12);
      pt->SetTextSize(0.04);
      pt->AddText("#bf{CMS} #it{simulation}");
      pt->Draw();
      
      pt=new TPaveText(0.65,0.955,0.95,0.99,"brNDC");
      pt->SetBorderSize(1);
      pt->SetFillStyle(0);
      pt->SetTextFont(42);
      pt->SetTextAlign(12);
      pt->SetTextSize(0.035);
      pt->AddText(title);
      pt->Draw();
    }
  
  //all done here
  c->SaveAs(outDir+"/"+name+".png");
}
