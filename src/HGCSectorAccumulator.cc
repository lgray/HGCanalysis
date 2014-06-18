#include "UserCode/HGCanalysis/interface/HGCSectorAccumulator.h"

//
HGCSectorAccumulator::HGCSectorAccumulator(int subDet,int layer, int copy) : gxH_(0), gyH_(0), gzH_(0), edepH_(0)
{
  //identifier for the histos
  sprintf(id_,"s_%d_%d_%d",subDet,layer,copy);
  sprintf(title_,"Layer %d Sector %d",layer,copy);
}


//
HGCSectorAccumulator::~HGCSectorAccumulator()
{
}

//
void HGCSectorAccumulator::configure(edm::Service<TFileService> &fs)
{
  //build the local -> global transformation

  //cache global quantities (deprecated)
  //float rho=sqrt(pow(gx_,2)+pow(gy_,2));

  //init sim histos
  int ndivx=TMath::Floor(bl_/cell_);
  ndivx=(ndivx+TMath::Floor((tl_-ndivx*cell_)/cell_));
  int ndivy=TMath::Floor(h_/cell_);
  edepH_ = fs->make<TH2F>(TString("E_")    +id_,title_+TString(";x [mm];y [mm]"),2*ndivx,-cell_*ndivx,cell_*ndivx,2*ndivy,-cell_*ndivy,cell_*ndivy);
  tH_    = fs->make<TH2F>(TString("AvgT_") +id_,title_+TString(";x [mm];y [mm]"),2*ndivx,-cell_*ndivx,cell_*ndivx,2*ndivy,-cell_*ndivy,cell_*ndivy);
  gxH_   = fs->make<TH2F>(TString("gx_")   +id_,title_+TString(";x [mm];y [mm]"),2*ndivx,-cell_*ndivx,cell_*ndivx,2*ndivy,-cell_*ndivy,cell_*ndivy);
  gyH_   = fs->make<TH2F>(TString("gy_")   +id_,title_+TString(";x [mm];y [mm]"),2*ndivx,-cell_*ndivx,cell_*ndivx,2*ndivy,-cell_*ndivy,cell_*ndivy);
  gzH_   = fs->make<TH2F>(TString("gz_")   +id_,title_+TString(";x [mm];y [mm]"),2*ndivx,-cell_*ndivx,cell_*ndivx,2*ndivy,-cell_*ndivy,cell_*ndivy);

  float slope=2*h_/(tl_-bl_);
  float offset=-h_*(tl_+bl_)/(tl_-bl_);
  for(int xbin=1; xbin<gxH_->GetXaxis()->GetNbins(); xbin++)
    {
      for(int ybin=1; ybin<gxH_->GetYaxis()->GetNbins(); ybin++)
	{
	  //the local coordinates
	  float localX=gxH_->GetXaxis()->GetBinCenter(xbin);
	  float localY=gxH_->GetYaxis()->GetBinCenter(ybin);

	  //exclude points outside the range
	  if(localY<fabs(localX)*slope+offset) continue;

	  //local->global (mm to cm) 
	  const HepGeom::Point3D<float> lcoord(localX/10,localY/10,0);
	  const HepGeom::Point3D<float> gcoord( local2globalTr_*lcoord );
	  gxH_->SetBinContent(xbin,ybin,gcoord.x()*10);
	  gyH_->SetBinContent(xbin,ybin,gcoord.y()*10);
	  gzH_->SetBinContent(xbin,ybin,gcoord.z()*10);
	}
    }

  //init reco histos
  if(recoCell_>0)
    {
      ndivx      = TMath::Floor(bl_/recoCell_);
      ndivx      = (ndivx+TMath::Floor((tl_-ndivx*recoCell_)/recoCell_));
      ndivy      = TMath::Floor(h_/recoCell_);
      adcH_      = fs->make<TH2F>(TString("ADC_") +id_,title_+TString(";x [mm];y [mm]"),  2*ndivx,-recoCell_*ndivx,recoCell_*ndivx,2*ndivy,-recoCell_*ndivy,recoCell_*ndivy);
      gxRecoH_   = fs->make<TH2F>(TString("recogx_")+id_,title_+TString(";x [mm];y [mm]"),2*ndivx,-recoCell_*ndivx,recoCell_*ndivx,2*ndivy,-recoCell_*ndivy,recoCell_*ndivy);
      gyRecoH_   = fs->make<TH2F>(TString("recogy_")+id_,title_+TString(";x [mm];y [mm]"),2*ndivx,-recoCell_*ndivx,recoCell_*ndivx,2*ndivy,-recoCell_*ndivy,recoCell_*ndivy);
      gzRecoH_   = fs->make<TH2F>(TString("recogz_")+id_,title_+TString(";x [mm];y [mm]"),2*ndivx,-recoCell_*ndivx,recoCell_*ndivx,2*ndivy,-recoCell_*ndivy,recoCell_*ndivy);
      for(int xbin=1; xbin<gxRecoH_->GetXaxis()->GetNbins(); xbin++)
	{
	  for(int ybin=1; ybin<gxRecoH_->GetYaxis()->GetNbins(); ybin++)
	    {
	      //local coordinates
	      float localX=gxRecoH_->GetXaxis()->GetBinCenter(xbin);
	      float localY=gxRecoH_->GetYaxis()->GetBinCenter(ybin);
	      
	      //local->global
	      const HepGeom::Point3D<float> lcoord(localX,localY,0);
	      const HepGeom::Point3D<float> gcoord( local2globalTr_*lcoord );
	      gxRecoH_->SetBinContent(xbin,ybin,gcoord.x());
	      gyRecoH_->SetBinContent(xbin,ybin,gcoord.y());
	      gzRecoH_->SetBinContent(xbin,ybin,gcoord.z());
	    }
	}
    }
  else
    {
      adcH_=0;
      gxRecoH_=0;
      gyRecoH_=0;
      gzRecoH_=0;
    }
}

//
int HGCSectorAccumulator::acquire(float edep, float t, float x, float y)
{
  if(edepH_==0) return -1;
  if(x>tl_ || x < -tl_ || y>h_ || y<-h_)
    {
      dumpGeometry();
      std::cout << " Can't accumulate @ (" << x << " " << y << ")" <<  std::endl;
      return 0;
    }
  tH_->Fill(x,y,t*edep);
  return edepH_->Fill(x,y,edep);
}

//
int HGCSectorAccumulator::digitize(float adc, float x, float y)
{
  if(adcH_==0) return -1;
  return adcH_->Fill(x,y,adc);
}


//
void HGCSectorAccumulator::reset()
{
  if(edepH_)  edepH_->Reset("ICE");
  if(tH_)     tH_->Reset("ICE");
  if(adcH_)   adcH_->Reset("ICE");
}

//
TVector3 HGCSectorAccumulator::getGlobalPointAt(int bin)
{
  TVector3 xyz(gxH_->GetBinContent(bin),
	       gyH_->GetBinContent(bin),
	       gzH_->GetBinContent(bin) 
	       );
  return xyz;
}

//
TVector3 HGCSectorAccumulator::getRecoGlobalPointAt(int bin)
{
  TVector3 xyz(gxRecoH_ ? gxRecoH_->GetBinContent(bin) : 0,
	       gyRecoH_ ? gyRecoH_->GetBinContent(bin) : 0,
	       gzRecoH_ ? gzRecoH_->GetBinContent(bin) : 0
	       );
  return xyz;
}

//
TVector2 HGCSectorAccumulator::getLocalPointAt(int bin)
{
  Int_t binx,biny,binz;
  gxH_->GetBinXYZ(bin,binx,biny,binz);
  TVector2 xy(gxH_->GetXaxis()->GetBinCenter(binx), gxH_->GetYaxis()->GetBinCenter(biny));
  return xy;
}

//
TVector2 HGCSectorAccumulator::getRecoLocalPointAt(int bin)
{
  TVector2 xy(0,0);
  if(gxRecoH_)
    {
      Int_t binx,biny,binz;
      gxRecoH_->GetBinXYZ(bin,binx,biny,binz);
      xy=TVector2(gxRecoH_->GetXaxis()->GetBinCenter(binx), gxRecoH_->GetYaxis()->GetBinCenter(biny));
    }
  return xy;
}

//
float HGCSectorAccumulator::getEnergyDepAt(int bin)
{
  return edepH_ ? edepH_->GetBinContent(bin) : 0.;
}

//
float HGCSectorAccumulator::getAverageTimeAt(int bin)
{
  return tH_ ? tH_->GetBinContent(bin) : 0.;
}

//
float HGCSectorAccumulator::getADCsAt(int bin)
{
  return adcH_ ? adcH_->GetBinContent(bin) : 0.;
}
