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
#include "TMultiGraph.h"
#include "TMath.h"
#include "TF1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TSystem.h"

#include <iostream>
#include <vector>

using namespace std;

Int_t comparisonColors[]={1,kAzure+2,kRed+2,kGreen-5,kGray+1};

void drawOccupancyResults(TString outDir="~/public/html/HGCal/Occupancy");
void fixExtremities(TH1* h,bool addOverflow, bool addUnderflow);
void simpleComparison(TString distName, TString distTitle,	std::vector<TString> &urlList, std::vector<TString> &subTitles, TString outDir);
void compareProfilesFor(TString distName,TString distTitle,Float_t yranMin,Float_t yranMax, bool drawLog,
			std::vector<Int_t> &thrScan,
			std::vector<TString> &urlList,
			std::vector<TString> &subTitles,
			TString outDir,
			bool isV4geometry=false);
void drawProfiles(TMultiGraph *gr,TMultiGraph *grEvol,
		  Float_t yranMin,Float_t yranMax,bool drawLog,
		  TString ytitle,TString name,
		  TString title,std::vector<int> &layerBoundaries,TString outDir);
void drawDistributions(std::vector< std::vector<TH1 *> > &plots,
		       std::vector<TString> subTitles,
		       TString name, TString ytitle, 
		       Float_t yranMin,Float_t yranMax, bool drawLog, 
		       TString title, TString outDir);

//
void drawOccupancyResults(TString outDir)
{
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  //prepare output directory
  gSystem->Exec("mkdir -p "+outDir);
  gSystem->Exec("cp $CMSSW_BASE/src/UserCode/HGCanalysis/test/analysis/occ_index.html "+ outDir+"/index.html");

  std::vector<TString> urlList,subTitles;
  urlList.push_back("/tmp/psilva/__step2.root_Occupancy_0.root");   subTitles.push_back("FE model");
  //  urlList.push_back("Single13_CMSSW_6_2_0_SLHC18_v2_tau_0_Occupancy_0.root");   subTitles.push_back("No shaping");
  //urlList.push_back("Single13_CMSSW_6_2_0_SLHC18_v2_tau_10_Occupancy_0.root");  subTitles.push_back("#tau=10 ns");
  //urlList.push_back("Single13_CMSSW_6_2_0_SLHC18_v2_tau_20_Occupancy_0.root");  subTitles.push_back("#tau=20 ns");

  //occupancies
  std::vector<Int_t> thrScan;
  thrScan.push_back(2);
  thrScan.push_back(4);
  thrScan.push_back(20);
  thrScan.push_back(40);
  //compareProfilesFor("occ",     "Occupancy / cell",             1e-3, 1,   true,  thrScan,   urlList, subTitles, outDir, false);
  
  //data volumes
  thrScan.clear();
  //  compareProfilesFor("datavol", "Expected readout cell data volume / event [bit]", 0, 12,  false,  thrScan,  urlList, subTitles, outDir, false);
  //compareProfilesFor("trigvol", "Expected trigger cell data volume / event [bit]", 0, 4,   false,  thrScan,  urlList, subTitles, outDir, false);
  //compareProfilesFor("mip",     "Energy [MIP]",           0,    25, false, thrScan,  urlList, subTitles, outDir, false);

  compareProfilesFor("busy", "Fraction of busy cells", 0, 0.001,   false,  thrScan,  urlList, subTitles, outDir, false);
  compareProfilesFor("tot", "Fraction of cells in ToT mode", 0, 0.001,   false,  thrScan,  urlList, subTitles, outDir, false);
  //simpleComparison("evtsize","Event size [log_{10} byte]",urlList,subTitles, outDir);
 }



//
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
void simpleComparison(TString distName, TString distTitle,	std::vector<TString> &urlList, std::vector<TString> &subTitles, TString outDir)
{
  std::vector< std::vector<TH1 *> > plots;
  for(size_t ifile=0; ifile<urlList.size(); ifile++)
    {
      TFile *fIn=TFile::Open(urlList[ifile]);
      TH2F *distH=(TH2F *)fIn->Get("analysis/"+distName);
      Int_t nybins(distH->GetYaxis()->GetNbins());
      if(ifile==0) plots=std::vector< std::vector<TH1 *> >(nybins);
      for(Int_t ybin=1; ybin<=nybins; ybin++)
	{
	  //project
	  TString projName(distH->GetName()); projName += "f"; projName += ifile; projName += "_iy"; projName += ybin; projName += "_proj";
    	  TH1D *h=distH->ProjectionX(distName+projName,ybin,ybin);
	  h->SetDirectory(0);
	  fixExtremities(h,true,true);
	  h->SetTitle( distH->GetYaxis()->GetBinLabel(ybin) );
	  h->SetLineWidth(2);
	  h->SetLineColor(comparisonColors[ifile]);
	  if(ifile==0)
	    {
	      h->SetFillStyle(1001);
	      h->SetFillColor(kGray);
	    }
	  else h->SetFillStyle(0);
	  plots[ybin-1].push_back(h);
	}
      fIn->Close();
    }
  
  drawDistributions(plots,subTitles,distName,distTitle,0,10,false,"",outDir);
}

//
void compareProfilesFor(TString distName,TString distTitle,Float_t yranMin,Float_t yranMax,bool drawLog,
			std::vector<Int_t> &thrScan,
			std::vector<TString> &urlList,
			std::vector<TString> &subTitles,
			TString outDir,
			bool isV4geometry)
{
  size_t nthr=thrScan.size();
  if(nthr==0) nthr=1;
  std::vector< std::vector< std::vector<TGraphErrors *> > > medianDistPerLayer,q95DistPerLayer; //file #, threshold #,eta bin #
  std::vector< std::vector< std::vector< std::vector<TH1 *> > > > distPerLayer; //file #, threshold #, layer #, eta bin #
  std::vector<int> layerBoundaries;
  if(isV4geometry) { layerBoundaries.push_back(31); layerBoundaries.push_back(12); layerBoundaries.push_back(10); }
  else             { layerBoundaries.push_back(30); layerBoundaries.push_back(12); layerBoundaries.push_back(12); }
  for(size_t ifile=0; ifile<urlList.size(); ifile++)
    {
      //open new file
      TFile *inF = TFile::Open(urlList[ifile]);
      medianDistPerLayer.push_back( std::vector< vector<TGraphErrors *> >(nthr) );
      q95DistPerLayer.push_back( std::vector< vector<TGraphErrors *> >(nthr) );
      distPerLayer.push_back( std::vector< std::vector< std::vector<TH1 *> > >(nthr) );
      
      //iterate over thresholds
      for(size_t ithr=0; ithr<nthr; ithr++)
	{	  
	  //built the profile graphs per layer (need a template from the file)
	  TH2F *distH= thrScan.size()==0?
	    (TH2F *) inF->Get("analysis/sd_0_layer1_"+distName) :
	    (TH2F *) inF->Get("analysis/sd_0_layer1_thr2_"+distName);
	  Int_t netaBins(distH->GetYaxis()->GetNbins());
	  for(Int_t ybin=1; ybin<=netaBins; ybin++)
	    {
	      TString grName(distName+"_ithr"); grName += ithr; grName+="_ieta"; grName+=ybin; grName += "_f"; grName += ifile;
	      Float_t etaMin=distH->GetYaxis()->GetBinLowEdge(ybin);
	      Float_t etaMax=distH->GetYaxis()->GetBinUpEdge(ybin);
	      char grTitle[50]; 
	      sprintf(grTitle,"%3.1f<|#eta|<%3.1f",etaMin,etaMax);
	      medianDistPerLayer[ifile][ithr].push_back(new TGraphErrors);
	      medianDistPerLayer[ifile][ithr][ybin-1]->SetName(grName);
	      medianDistPerLayer[ifile][ithr][ybin-1]->SetTitle(grTitle);
	      medianDistPerLayer[ifile][ithr][ybin-1]->SetMarkerStyle(20+(ybin-1)/2);
	      medianDistPerLayer[ifile][ithr][ybin-1]->SetLineColor(27+ybin*2);
	      medianDistPerLayer[ifile][ithr][ybin-1]->SetMarkerColor(27+ybin*2);
	      medianDistPerLayer[ifile][ithr][ybin-1]->SetMarkerSize(1.0);
	      medianDistPerLayer[ifile][ithr][ybin-1]->SetFillStyle(0);
	      medianDistPerLayer[ifile][ithr][ybin-1]->SetFillColor(27+ybin*2);
	      q95DistPerLayer[ifile][ithr].push_back( (TGraphErrors *)medianDistPerLayer[ifile][ithr][ybin-1]->Clone("width_"+grName) );
	    }

	  //iterate over sub-detectors
	  int layerCtr(0);
	  for(size_t isd=0;isd<=2; isd++)
	    {
	      Int_t nlayers=layerBoundaries[isd];
	      for(Int_t ilay=1; ilay<=nlayers; ilay++)
		{
		  layerCtr++;
		  TString pfix("sd_"); pfix+=isd; pfix += "_layer"; pfix += ilay; 
		  if(thrScan.size()) { pfix += "_thr"; pfix += thrScan[ithr]; }
		  distH= (TH2F *)inF->Get("analysis/"+pfix+"_"+distName);

		  //profile at different etas
		  distPerLayer[ifile][ithr].push_back( std::vector<TH1 *>(netaBins) );
		  for(Int_t ybin=1; ybin<=netaBins; ybin++)
		    {
		      //project
		      TString projName(distH->GetName()); projName += "f"; projName += ifile; projName += "_iy"; projName += ybin; projName += "_proj";
		      TH1D *h=distH->ProjectionX(projName,ybin,ybin);
		      h->SetDirectory(0);
		      fixExtremities(h,true,true);
		      char hTitle[50]; 
		      TString sdName("EE");
		      if(isd==1) sdName="HEF";
		      if(isd==2) sdName="HEB";
		      sprintf(hTitle,"%s, layer #%d, %3.1f<|#eta|<%3.1f",
			      sdName.Data(),
			      ilay,
			      distH->GetYaxis()->GetBinLowEdge(ybin),distH->GetYaxis()->GetBinUpEdge(ybin));
	      	      h->SetTitle(hTitle);
		      h->SetLineWidth(2);
		      h->SetLineColor(comparisonColors[ifile]);
		      if(ifile==0){
			h->SetFillStyle(1001);
			h->SetFillColor(kGray);
		      }
		      else h->SetFillStyle(0);
		      distPerLayer[ifile][ithr][layerCtr-1][ybin-1]=h;
		      
		      //compute quantiles for the profiles
		      Double_t xq[3]={0.05,0.5,0.95};
		      Double_t yq[3];
		      h->GetQuantiles(3,yq,xq);
		      Int_t np=medianDistPerLayer[ifile][ithr][ybin-1]->GetN();
		      medianDistPerLayer[ifile][ithr][ybin-1]->SetPoint(np,layerCtr,yq[1]);
		      medianDistPerLayer[ifile][ithr][ybin-1]->SetPointError(np,0,1.253*h->GetMeanError());
		      if(yq[2]>0){
			np=q95DistPerLayer[ifile][ithr][ybin-1]->GetN();
			q95DistPerLayer[ifile][ithr][ybin-1]->SetPoint(np,layerCtr,yq[2]);
			q95DistPerLayer[ifile][ithr][ybin-1]->SetPointError(np,0,2*1.253*h->GetMeanError());
		      }
		    }
		}
	    }
	}
    }

  //now let's plot all of this...
  //distributions
  for(size_t ithr=0; ithr<nthr; ithr++)
    {
      size_t nlay=distPerLayer[0][ithr].size();
      for(size_t ilay=0; ilay<nlay; ilay++)
	{
	  size_t netabins=distPerLayer[0][ithr][ilay].size();
	  std::vector< std::vector<TH1 *> > plots(netabins); //eta bin #, file #
	  for(size_t ietabin=0; ietabin<netabins; ietabin++)
	    {
	      size_t nfiles(distPerLayer.size());
	      plots[ietabin]=std::vector<TH1 *>(nfiles);
	      for(size_t ifile=0; ifile<distPerLayer.size(); ifile++)
		{
		  plots[ietabin][ifile]=distPerLayer[ifile][ithr][ilay][ietabin];
		}
	    }
	  TString title("Layer #"); title += (ilay+1);
	  if(thrScan.size())
	    {
	      char buftitle[50];
	      sprintf(buftitle,"Threshold: %3.1f MIPs",thrScan[ithr]*0.25);
	      title+=","; title+=TString(buftitle);
	    }
	  TString pfix("_"); pfix+=ilay;
	  if(thrScan.size()) { pfix += "_thr"; pfix += thrScan[ithr]; }
	  drawDistributions(plots,subTitles,distName+pfix, distTitle, yranMin,yranMax,drawLog, title,outDir);
	}
    }

  //distribution profiles
  for(size_t iplot=0; iplot<2; iplot++)
    {
      std::vector< std::vector< vector<TGraphErrors *> > > *allPlots=0;
      TString yTitle("");
      if(iplot==0) { allPlots=&medianDistPerLayer;      yTitle=distTitle+" median";       }
      if(iplot==1) { allPlots=&q95DistPerLayer;         yTitle=distTitle+" 95% quantile"; }
      
      for(size_t ithr=0; ithr<nthr; ithr++)
	{
	  TMultiGraph *oddEtaGraphs     = new TMultiGraph;
	  TMultiGraph *evenEtaGraphs    = new TMultiGraph;
	  TMultiGraph *oddEtaGraphEvol  = urlList.size()>0 ? new TMultiGraph : 0; 
	  TMultiGraph *evenEtaGraphEvol = urlList.size()>0 ? new TMultiGraph : 0;
	  for(size_t ieta=0; ieta<(*allPlots)[0][ithr].size(); ieta++)
	    {
	      if(ieta%2==0) evenEtaGraphs->Add( (*allPlots)[0][ithr][ieta], "p" );
	      else         oddEtaGraphs->Add( (*allPlots)[0][ithr][ieta], "p" );

	      Double_t ix,iy;
	      for(Int_t ip=0; ip<(*allPlots)[0][ithr][ieta]->GetN(); ip+=4)
		{
		  TGraph *gr=new TGraph;
		  gr->SetLineColor( (*allPlots)[0][ithr][ieta]->GetLineColor() );
		  gr->SetMarkerColor( (*allPlots)[0][ithr][ieta]->GetMarkerColor() );
		  gr->SetMarkerStyle( (*allPlots)[0][ithr][ieta]->GetMarkerStyle() );
		  gr->SetMarkerSize(0.8);
		  gr->SetFillStyle(0);
		  for(size_t ifile=0; ifile<urlList.size(); ifile++)
		    {
		      (*allPlots)[ifile][ithr][ieta]->GetPoint(ip,ix,iy);
		      gr->SetPoint(ifile,ix,iy);
		    }
		  if(ieta%2==0) evenEtaGraphEvol->Add(gr,"lp");
		  else          oddEtaGraphEvol->Add(gr,"lp");
		}
	    }
	  
	  TString name("summary"); 
	  if(thrScan.size()) { name +="_"; name+=Int_t(thrScan[ithr]*0.25); }
	  if(iplot==0) name=distName+"median_"+name;
	  if(iplot==1) name=distName+"q95_"+name;
	  TString title("");
	  if(thrScan.size())
	    {
	      char buftitle[50];
	      sprintf(buftitle,"Threshold: %3.1f MIPs",thrScan[ithr]*0.25);
	      title=TString(buftitle);
	    }
	  drawProfiles(oddEtaGraphs,  oddEtaGraphEvol,  yranMin, yranMax, drawLog, yTitle, name+"_a", title, layerBoundaries, outDir);
	  drawProfiles(evenEtaGraphs, evenEtaGraphEvol, yranMin, yranMax, drawLog, yTitle, name+"_b", title, layerBoundaries, outDir);
	}
    }
}


//
void drawProfiles(TMultiGraph *gr,TMultiGraph *grEvol,Float_t yranMin,Float_t yranMax, bool drawLog, TString ytitle,TString name,TString title,std::vector<int> &layerBoundaries,TString outDir)
{
  TCanvas *c=new TCanvas("c","c",800,500);
  c->SetLeftMargin(0.15);
  c->SetTopMargin(0.05);
  c->SetRightMargin(0.05);
  c->SetBottomMargin(0.12);
  if(drawLog) c->SetLogy();
  c->SetGridy();


  //draw
  gr->Draw("a");
  gr->GetYaxis()->SetTitle(ytitle);
  gr->GetXaxis()->SetTitle("Layer number");
  gr->GetYaxis()->SetRangeUser(yranMin,yranMax);
  TLegend *leg=c->BuildLegend();
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetNColumns(2);
  if(grEvol) grEvol->Draw("l");
  
  TPaveText *pt=new TPaveText(0.13,0.96,0.6,0.99,"brNDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  pt->SetTextFont(42);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.04);
  pt->AddText("#bf{CMS} #it{simulation}");
  pt->Draw();

  Int_t nlayers(0);
  TLatex *txt=new TLatex;
  txt->SetTextFont(42);
  txt->SetTextSize(0.04);
  txt->SetTextColor(kBlue);
  for(size_t i=0; i<layerBoundaries.size(); i++)
    {
      TString sdName("EE");
      if(i==1) sdName="HEF";
      if(i==2) sdName="HEB";
      txt->DrawLatex( nlayers+layerBoundaries[i]*0.4, drawLog ? 0.7*yranMax : 0.95*yranMax, sdName);

      nlayers+=layerBoundaries[i];
      TLine *l=new TLine(nlayers,yranMin,nlayers,yranMax);
      l->SetLineColor(kBlue);
      l->SetLineStyle(7);
      l->Draw();
    }

  if(title!="")
    {
      pt=new TPaveText(0.65,0.955,0.95,0.99,"brNDC");
      pt->SetBorderSize(1);
      pt->SetFillStyle(0);
      pt->SetTextFont(42);
      pt->SetTextAlign(12);
      pt->SetTextSize(0.035);
      pt->AddText(title);
      pt->Draw();
    }
    
  c->SaveAs(outDir+"/"+name+".png");
}

//
void drawDistributions(std::vector< std::vector<TH1 *> > &plots,
		       std::vector<TString> subTitles,
		       TString name, TString ytitle, 
		       Float_t yranMin,Float_t yranMax, bool drawLog, 
		       TString title, TString outDir)
{
  TCanvas *c=new TCanvas("c","c",2000,1000);

  size_t nx=plots.size()/2;
  while( 2*nx<plots.size() ) nx++;
  c->Divide(nx,2);

  for(size_t i=0; i<plots.size(); i++)
    {
      TPad *p=(TPad *)c->cd(i+1);
      p->SetLeftMargin(0.15);
      p->SetTopMargin(0.05);
      p->SetRightMargin(0.05);
      p->SetBottomMargin(0.12);
      //if(drawLog) p->SetLogy();
      p->SetGridy();
      TLegend *leg=new TLegend(0.65,0.7,0.95,0.9);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextFont(42);
      leg->SetTextSize(0.04);
      float ymax(-999999999.);
      for(size_t j=0; j<plots[i].size(); j++)
	{
	  plots[i][j]->Draw(j==0 ? "hist": "histsame");
	  float iymax=plots[i][j]->GetBinContent(plots[i][j]->GetMaximumBin());
	  ymax=TMath::Max(ymax,iymax);
	  leg->AddEntry(plots[i][j],subTitles[j],"f");
	}
      leg->Draw();
      plots[i][0]->GetXaxis()->SetTitle(ytitle);
      plots[i][0]->GetYaxis()->SetTitle("Probability distribution");
      plots[i][0]->GetYaxis()->SetRangeUser(0,ymax*1.1);
      plots[i][0]->GetXaxis()->SetRangeUser(yranMin,yranMax);
 

      TPaveText *pt=new TPaveText(0.5,0.955,0.95,0.99,"brNDC");
      pt->SetBorderSize(1);
      pt->SetFillStyle(0);
      pt->SetTextFont(42);
      pt->SetTextAlign(12);
      pt->SetTextSize(0.035);
      pt->AddText(plots[i][0]->GetTitle());
      pt->Draw();

      //if(i) continue;
      pt=new TPaveText(0.13,0.96,0.6,0.99,"brNDC");
      pt->SetBorderSize(0);
      pt->SetFillStyle(0);
      pt->SetTextFont(42);
      pt->SetTextAlign(12);
      pt->SetTextSize(0.04);
      pt->AddText("#bf{CMS} #it{simulation}");
      pt->Draw();
    }
  
  c->SaveAs(outDir+"/"+name+".png");
}
