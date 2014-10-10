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

void showOccupancies(TObjArray plots,TString ytitle,TString name,TString title,std::vector<int> &layerBoundaries,TString outDir);
void compareOccupancyResults(std::vector<TString> &urlList,TString outDir);
void fixExtremities(TH1* h,bool addOverflow, bool addUnderflow);

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



void drawOccupancyResults(TString outDir="~/www/HGCal/Occupancy/v5")
{
  //gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  std::vector<TString> urlList;
  urlList.push_back("/tmp/psilva/Single13_CMSSW_6_2_0_SLHC18_v2_tau_0_Occupancy_0.root");
  compareOccupancyResults(urlList,outDir);
 }

//
void compareOccupancyResults(std::vector<TString> &urlList,TString outDir)
{
  Int_t adc_thr[]={2,4,20,40};
  typedef std::pair<int,int> ThresholdKey_t; //file,threshold
  std::map<ThresholdKey_t, vector<TGraphErrors *> > occPerLayer,occwidthPerLayer; //threshold,eta profiles
  std::vector<int> layerBoundaries;
  //v4
  //layerBoundaries.push_back(31); layerBoundaries.push_back(12); layerBoundaries.push_back(10);
  //v5
  layerBoundaries.push_back(30); layerBoundaries.push_back(12); layerBoundaries.push_back(12);

  for(size_t ifile=0; ifile<urlList.size(); ifile++)
    {
      TFile *inF = TFile::Open(urlList[ifile]);

      for(size_t ithr=0; ithr<sizeof(adc_thr)/sizeof(Int_t); ithr++)
	{
	  ThresholdKey_t thr_key(ifile,adc_thr[ithr]);
	  if(occPerLayer.find(thr_key)==occPerLayer.end())
	    {
	      vector<TGraphErrors *> templateVec;
	      occPerLayer[thr_key]=templateVec;
	      occwidthPerLayer[thr_key]=templateVec;
	    }
	  
	  int layerCtr(0);
	  for(size_t isd=0;isd<=2; isd++)
	    {
	      Int_t nlayers=layerBoundaries[isd];

	      for(Int_t ilay=1; ilay<=nlayers; ilay++)
		{
		  layerCtr++;
		  TString pfix("sd_"); pfix+=isd; pfix += "_layer"; pfix += ilay; pfix += "_thr"; pfix += adc_thr[ithr];
		  TH2F *occH= (TH2F *)inF->Get("analysis/"+pfix+"_occ");

		  //init histograms
		  if(layerCtr==1)
		    {
		      for(Int_t ybin=1; ybin<=occH->GetYaxis()->GetNbins(); ybin++)
			{
			  occPerLayer[thr_key].push_back(new TGraphErrors);
			  Float_t etaMin=occH->GetYaxis()->GetBinLowEdge(ybin);
			  Float_t etaMax=occH->GetYaxis()->GetBinUpEdge(ybin);
			  char buf[50]; 
			  sprintf(buf,"%3.1f<|#eta|<%3.1f",etaMin,etaMax);
			  occPerLayer[thr_key][ybin-1]->SetTitle(buf);
			  occPerLayer[thr_key][ybin-1]->SetMarkerStyle(20+(ybin-1)%4);
			  occPerLayer[thr_key][ybin-1]->SetLineColor((ybin-1)%4==0 ? 1 : 20+10*(ybin-1)%4);
			  occPerLayer[thr_key][ybin-1]->SetMarkerColor((ybin-1)%4==0 ? 1 : 20+10*(ybin-1)%4);
			  occwidthPerLayer[thr_key].push_back( (TGraphErrors *)occPerLayer[thr_key][ybin-1]->Clone() );
			}
		    }

		  //profile at different etas
		  for(Int_t ybin=1; ybin<=occH->GetYaxis()->GetNbins(); ybin++)
		    {
		      TH1D *h=occH->ProjectionX("occproj",ybin,ybin);
		      fixExtremities(h,true,true);
		      Double_t xq[3]={0.05,0.5,0.95};
		      Double_t yq[3];
		      h->GetQuantiles(3,yq,xq);

		      Int_t np=occPerLayer[thr_key][ybin-1]->GetN();
		      occPerLayer[thr_key][ybin-1]->SetPoint(np,layerCtr,yq[1]);
		      occPerLayer[thr_key][ybin-1]->SetPointError(np,0,1.253*h->GetRMS()/TMath::Sqrt(h->Integral()));
		      occwidthPerLayer[thr_key][ybin-1]->SetPoint(np,layerCtr,yq[2]-yq[1]);
		      occwidthPerLayer[thr_key][ybin-1]->SetPointError(np,0,h->GetRMSError());
		      
		    }
		}
	    }
	}
    }

  for(std::map<ThresholdKey_t, vector<TGraphErrors *> >::iterator it=occPerLayer.begin();
      it!=occPerLayer.end();
      it++)
    {
      TObjArray plots,widthPlots; 
      for(size_t i=0; i<it->second.size(); i++) 
	{
	  plots.Add( (it->second)[i] );
	  widthPlots.Add( (occwidthPerLayer[it->first])[i] );
	}

      TString name("occsummary_"); name+= it->first.second;
      char buf[50];
      sprintf(buf,"Threshold: %3.1f MIPs",it->first.second*0.25);
      TString title(buf);
      showOccupancies(plots,"Median occupancy",name,title,layerBoundaries,outDir);
      showOccupancies(widthPlots,"Occupancy width (q_{95}-q_{50})","width_"+name,title,layerBoundaries,outDir);
    }
}


//
void showOccupancies(TObjArray plots,TString ytitle,TString name,TString title,std::vector<int> &layerBoundaries,TString outDir)
{
  TCanvas *c=new TCanvas("c","c",800,500);
  c->SetLeftMargin(0.15);
  c->SetTopMargin(0.05);
  c->SetRightMargin(0.05);
  c->SetBottomMargin(0.12);
  c->SetLogy();
  c->SetGridy();

  TLegend *leg=new TLegend(0.15,0.8,0.7,0.94);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);

  //draw
  for(Int_t i=0; i<plots.GetEntriesFast(); i++)
    {
      TGraphErrors *gr=(TGraphErrors *)plots.At(i);
      gr->Draw(i==0 ? "ap" : "p");
      gr->GetXaxis()->SetTitle("HGC layer number");
      gr->GetYaxis()->SetRangeUser(1e-2,1);
      gr->GetXaxis()->SetTitleSize(0.05);
      gr->GetXaxis()->SetLabelSize(0.04);
      gr->GetYaxis()->SetTitleSize(0.05);
      gr->GetYaxis()->SetLabelSize(0.04);
      gr->GetYaxis()->SetTitle(ytitle);
      leg->AddEntry(gr,gr->GetTitle(),"p");
    }
  leg->Draw();
  
  TPaveText *pt=new TPaveText(0.12,0.96,0.6,0.99,"brNDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  pt->SetTextFont(42);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.04);
  pt->AddText("CMS simulation");
  pt->Draw();

  Int_t nlayers(0);
  for(size_t i=0; i<layerBoundaries.size(); i++)
    {
      nlayers+=layerBoundaries[i];
      TLine *l=new TLine(nlayers,1e-2,nlayers,1);
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
