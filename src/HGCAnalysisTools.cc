#include "UserCode/HGCanalysis/interface/HGCAnalysisTools.h"
#include <vector>

#include "TVirtualFitter.h"

using namespace std;

//
float getLambdaForHGCLayer(int hit_layer)
{
  if (hit_layer==1)                       return 0.01; 
  else if(hit_layer>=2 && hit_layer<=11)  return 0.036;
  else if(hit_layer>=12 && hit_layer<=21) return 0.043;
  else if(hit_layer>=22 && hit_layer<=30) return 0.056;
  else if(hit_layer==31)                  return 0.338;
  else if(hit_layer>=32 && hit_layer<=42) return 0.273;
  else if(hit_layer>42)                   return 0.475;
  return 0;
}


//
G4InteractionPositionInfo getInteractionPosition(const std::vector<SimTrack> *SimTk, 
						 const std::vector<SimVertex> *SimVtx, 
						 int barcode)
{
  G4InteractionPositionInfo toRet;
  toRet.pos=math::XYZVectorD(0,0,0);
  toRet.info=0;

  //loop over vertices
  for (const SimVertex &simVtx : *SimVtx) 
    {
      //require the parent to be the given barcode
      bool noParent( simVtx.noParent() );
      if(noParent) continue;
      int pIdx( simVtx.parentIndex() );
      if( pIdx!=barcode) continue;

      int vtxIdx(simVtx.vertexId());
      int rawTkMult(0),eTkMult(0),gTkMult(0),nTkMult(0),nucleiTkMult(0),pTkMult(0);
      for (const SimTrack &vtxTk : *SimTk)
	{
	  int tkVtxIdx( vtxTk.vertIndex() ); 
	  if(tkVtxIdx!=vtxIdx) continue;
	  
	  int tkType=vtxTk.type();
	  rawTkMult++;
	  eTkMult      += (abs(tkType)==11);
	  gTkMult      += (abs(tkType)==22);
	  nTkMult      += (abs(tkType)==2112 || abs(tkType)==2212);
	  nucleiTkMult += (abs(tkType)>1000000000);
	  pTkMult      += (abs(tkType)==211 || abs(tkType)==111);
	}
      
      if(rawTkMult<2) continue;

      toRet.pos=math::XYZVectorD(simVtx.position());
      if(rawTkMult==3 && nTkMult==1 && nucleiTkMult==1 && pTkMult==1) toRet.info=1;
      return toRet;
    }

  return toRet;
}


//
std::pair<float,float> getEffSigma(RooRealVar *var, RooAbsPdf *pdf, float wmin,float wmax, float step, float epsilon)
{
  //get cdf points
  RooAbsReal *cdf = pdf->createCdf(RooArgList(*var));
  float point=wmin;
  vector<pair<float,float> > points;
  while (point <= wmax){
    var->setVal(point);
    if (pdf->getVal() > epsilon){
      points.push_back(pair<float,float>(point,cdf->getVal()));
    }
    point+=step;
  }

  float low = wmin;
  float high = wmax;
  float width = wmax-wmin;
  for (unsigned int i=0; i<points.size(); i++){
    for (unsigned int j=i; j<points.size(); j++){
      float wy = points[j].second - points[i].second;
      if (TMath::Abs(wy-0.683) < epsilon){
	float wx = points[j].first - points[i].first;
	if (wx < width){
	  low = points[i].first;
	  high = points[j].first;
	  width=wx;
	}
      }
    }
  }
  pair<float,float> result(low,high);
  return result;
}


//
CircleParams_t fitCircleTo(TGraphErrors *gr)
{
  CircleParams_t result;
  result.isvalid=false;
  if(gr==0 || gr->GetN()<3) return result;
  result.ndf=gr->GetN()-3;
  
  //Fit a circle to the graph points
  TVirtualFitter::SetDefaultFitter("Minuit");  //default is Minuit
  TVirtualFitter *fitter = TVirtualFitter::Fitter(0, 3+2*gr->GetN());
  fitter->SetFCN(circle_fcn);
  fitter->SetParameter(0, "x0",   0, 0.1, 0,0);
  fitter->SetParameter(1, "y0",   0, 0.1, 0,0);
  fitter->SetParameter(2, "R",    1, 0.1, 0,0);
  for(int ip=0; ip<gr->GetN(); ip++)
    {
      Double_t x,y;
      gr->GetPoint(ip,x,y);
      Double_t xErr=gr->GetErrorX(ip);
      Double_t yErr=gr->GetErrorY(ip);
      TString pfix("_"); pfix+=ip;
      fitter->SetParameter(3+ip*4+0,"x"+pfix,x,0,x,x);
      fitter->FixParameter(3+ip*4+0);
      fitter->SetParameter(3+ip*4+1,"xerr"+pfix,xErr,0,xErr,xErr);
      fitter->FixParameter(3+ip*4+1);
      fitter->SetParameter(3+ip*4+2,"y"+pfix,y,0,y,y);
      fitter->FixParameter(3+ip*4+2);
      fitter->SetParameter(3+ip*4+3,"yerr"+pfix,yErr,0,yErr,yErr);
      fitter->FixParameter(3+ip*4+3);
    }

  Double_t arglist[1] = {0};
  fitter->ExecuteCommand("MIGRAD", arglist, 0);

  result.x0      = fitter->GetParameter(0);
  result.x0_err  = fitter->GetParError(0);
  result.y0      = fitter->GetParameter(1);
  result.y0_err  = fitter->GetParError(1);
  result.r       = fitter->GetParameter(2);
  result.r_err   = fitter->GetParError(2);
  result.chi2    = 0; //fitter->Chisquare(3);
  result.isvalid = true;
  return result;
}


//
void circle_fcn(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) 
{
  //minimisation function computing the sum of squares of residuals
  Int_t np = sizeof(par)/sizeof(Double_t);
  np-=3;
  f = 0;
  for (Int_t ip=0;ip<np;ip++)
    {
      Double_t x    = par[3+4*ip+0];
      Double_t xErr = par[3+4*ip+1];
      Double_t y    = par[3+4*ip+2];
      Double_t yErr = par[3+4*ip+3];

      Double_t u     = x - par[0];
      Double_t v     = y - par[1];
      Double_t r     = TMath::Sqrt(u*u+v*v);;
      Double_t dr    = par[2] - r;
      Double_t drErr = (xErr*u+yErr*v)/r;

      if(drErr>0) f += pow(dr/drErr,2);
      else        f += pow(dr,2);
    }
}

