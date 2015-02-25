#include "UserCode/HGCanalysis/plugins/HGCGeometryAnalyzer.h"

#include "DetectorDescription/OfflineDBLoader/interface/GeometryInfoDump.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

#include "DetectorDescription/Core/interface/DDFilter.h"
#include "DetectorDescription/Core/interface/DDFilteredView.h"
#include "DetectorDescription/Core/interface/DDSolid.h"

#include "CLHEP/Units/GlobalSystemOfUnits.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Plane3D.h"
#include "CLHEP/Geometry/Vector3D.h"

#include "TF1.h"

#include <iostream>

#define SEP "\t\t\t\t"

using namespace std;

//
HGCGeometryAnalyzer::HGCGeometryAnalyzer( const edm::ParameterSet &iConfig )
  : testDone_(false)
{
  geometrySource_   = iConfig.getUntrackedParameter< std::vector<std::string> >("geometrySource");
  edm::Service<TFileService> fs;
  fs_=&fs;
}

//
HGCGeometryAnalyzer::~HGCGeometryAnalyzer()
{
}
  
//
void HGCGeometryAnalyzer::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  if(testDone_) return;

  //read geometry from event setup
  for(size_t i=0; i<geometrySource_.size(); i++)
    {
      edm::ESHandle<HGCalGeometry> hgcGeo;
      iSetup.get<IdealGeometryRecord>().get(geometrySource_[i],hgcGeo);

      const HGCalTopology &topo=hgcGeo->topology();
      const HGCalDDDConstants &dddConst=topo.dddConstants();

      int simLayers( dddConst.layers(false) );
      int recLayers( dddConst.layers(true) );

      //REPORT

      cout << geometrySource_[i] << " " << SEP   << "SIM "  <<  SEP  << "REC "  <<  SEP  << endl
	   << "-----------------------------------------------------------------------------" << endl
	   << "#layers " << SEP  << simLayers << SEP  << recLayers << SEP  << endl
	   << "-----------------------------------------------------------------------------" << endl
	   << "The following lines show the correspondence between sim and rec layers (#cells, #rows, cell size, b, t, h, alpha)" << endl
	   << "b/t=trapezoid half distance at the bottom/top  h=trapezoid height  alpha=half sector degree (0 if full sector)" << endl
	   << "Notice SIM/REC units are in mm/cm" << endl
	   << "-----------------------------------------------------------------------------" << endl;
      
      for(int ilay=1; ilay<=simLayers; ilay++)
	{
	  int simcells=dddConst.maxCells(ilay,false);
	  int simrows=dddConst.maxRows(ilay,false);
	  std::vector<HGCalDDDConstants::hgtrap>::const_iterator simModIt( dddConst.getFirstModule(false) );
	  std::vector<HGCalDDDConstants::hgtrform>::const_iterator trformIt( dddConst.getFirstTrForm() );
	  for(int klay=1; klay<ilay; klay++) { simModIt++; trformIt++; }
	  std::cout.precision(4);
	  cout << " " << SEP  
	       << ilay << " : " 
	       << simcells << "/" 
	       << simrows << "/" 
	       << simModIt->cellSize << "/" 
	       << simModIt->bl << "/"
	       << simModIt->tl << "/"
	       << simModIt->h << "/"
	       << simModIt->alpha 
	       << SEP ;

	  std::pair<int,int>  simToReco=dddConst.simToReco(1,ilay,false);
	  std::vector<HGCalDDDConstants::hgtrap>::const_iterator recModIt( dddConst.getFirstModule(true) );
	  if(simToReco.second>0) 
	    {
	      int reccells=dddConst.maxCells(simToReco.second,true);
	      int recrows=dddConst.maxRows(simToReco.second,true);
	      for(int klay=1; klay<simToReco.second; klay++) recModIt++;
	      cout << simToReco.second << " : " 
		   << reccells << "/" 
		   << recrows << "/" 
		   << recModIt->cellSize << "/" 
		   << recModIt->bl << "/"
		   << recModIt->tl << "/"
		   << recModIt->h << "/"
		   << recModIt->alpha 
		   << "/" << trformIt->h3v.perp() << "/" << trformIt->h3v.perp()-recModIt->h << "/" << trformIt->h3v.perp()+recModIt->h << "/" << trformIt->h3v.z();
	      
	      if(fabs(simModIt->bl-10*recModIt->bl)>1e-3 ||
		 fabs(simModIt->tl-10*recModIt->tl)>1e-3 ||
		 fabs(simModIt->h-10*recModIt->h)>1e-3 )
		{
		  cout << " SIM->RECO NOT OK" << endl;
		}
	      cout << SEP;
	    }
	  else cout << "x" << SEP;
	  
	  cout << endl;

	  //test SIM->RECO coordinate assignment
	  int recLay(simToReco.second);
	  if(recLay<=0) continue;

	  TString name("layer"); name+= ilay; name +="_sd"; name+=i;

	  //sim
	  int ncellsy          = 2*simModIt->h/simModIt->cellSize;
	  int ncellsx          = 2*simModIt->tl/simModIt->cellSize;
	  TString funcFormula  = "[0]*TMath::Abs(x)+[1]";
	  float func_slope     = 2*simModIt->h/(simModIt->tl-simModIt->bl);
	  float func_offset    = -2*simModIt->h*simModIt->bl/(simModIt->tl-simModIt->bl);
	  float xmin           = -simModIt->tl;
	  float xmax           = simModIt->tl;
	  float ymin           = 0;
	  float ymax           = 2*simModIt->h;
	  float boundaryXmin   = -simModIt->tl;
	  float boundaryXmax   = simModIt->tl;
	  float translationX   = 0;
	  float translationY   = simModIt->h;
	  if(simModIt->alpha>0)
	    {
	      ncellsx      = 2*simModIt->tl/simModIt->cellSize;
	      funcFormula  = "[0]*x+[1]";
	      func_slope   = simModIt->h/(simModIt->tl-simModIt->bl);
	      func_offset  = -2*simModIt->h*simModIt->bl/(simModIt->tl-simModIt->bl);
	      xmin         = 0;
	      xmax         = 2*simModIt->tl;
	      ymin         = 0;
	      ymax         = 2*simModIt->h;
	      boundaryXmin = 2*simModIt->bl;
	      boundaryXmax = xmax;
	      translationX = 0.5*(simModIt->tl+simModIt->bl);
	      translationY = simModIt->h;
	    }
	  if(simModIt->alpha<0)
	    {
	      ncellsx      = 2*simModIt->tl/simModIt->cellSize;
	      funcFormula  = "[0]*x+[1]";
	      func_slope   = -simModIt->h/(simModIt->tl-simModIt->bl);
	      func_offset  = -2*simModIt->h*simModIt->bl/(simModIt->tl-simModIt->bl);
	      xmin         = -2*simModIt->tl;
	      xmax         = 0;
	      ymin         = 0;
	      ymax         = 2*simModIt->h;
	      boundaryXmin = xmin;
	      boundaryXmax = -2*simModIt->bl;
	      translationX = -0.5*(simModIt->tl+simModIt->bl);
	      translationY = simModIt->h;
	    }
	  
	  TF1 *func=(*fs_)->make<TF1>("boundary_"+name,funcFormula,boundaryXmin,boundaryXmax);
	  func->SetParameter(0,func_slope);
	  func->SetParameter(1,func_offset);
	  func->SetLineWidth(2);
	  func->SetLineColor(1);
	  func->SetLineStyle(7);

	  TH2F *simCellH = (*fs_)->make<TH2F>("simcell_"+ name, ";Local x [mm];Local y [mm];Sim cell number",
					      ncellsx,  xmin, xmax,
					      ncellsy,  ymin, ymax);
	  TH2F *simixH = (TH2F *)simCellH->Clone("simix_"+name ); simixH = (*fs_)->make<TH2F>( *simixH ); simixH->GetZaxis()->SetTitle("i-x");
	  TH2F *simiyH = (TH2F *)simCellH->Clone("simiy_"+name ); simiyH = (*fs_)->make<TH2F>( *simiyH ); simiyH->GetZaxis()->SetTitle("i-y");
	  
	  //reco 
	  TH2F *recCellH = (TH2F *)simCellH->Clone("reccell_"+name ); recCellH = (*fs_)->make<TH2F>( *recCellH );
	  TH2F *recixH = (TH2F *)simCellH->Clone("recix_"+name ); recixH = (*fs_)->make<TH2F>( *recixH ); recixH->GetZaxis()->SetTitle("i-x");
	  TH2F *reciyH = (TH2F *)simCellH->Clone("reciy_"+name ); reciyH = (*fs_)->make<TH2F>( *reciyH ); reciyH->GetZaxis()->SetTitle("i-y");

	  //sim-reco diff
	  TH2F *dxH = (TH2F *)simCellH->Clone("dx_"+name ); dxH = (*fs_)->make<TH2F>( *dxH ); dxH->GetZaxis()->SetTitle("#Delta x(SIM, RECO) [mm]");
	  TH2F *dyH = (TH2F *)simCellH->Clone("dy_"+name ); dyH = (*fs_)->make<TH2F>( *dyH ); dyH->GetZaxis()->SetTitle("#Delta y(SIM, RECO) [mm]");
	  
	  for(int xbin=1; xbin<=simCellH->GetXaxis()->GetNbins(); xbin++)
	    for(int ybin=1; ybin<=simCellH->GetYaxis()->GetNbins(); ybin++)
	      {
		float x    = simCellH->GetXaxis()->GetBinCenter(xbin)-translationX;
		float y    = simCellH->GetYaxis()->GetBinCenter(ybin)-translationY;

		int subsec(x<0 ? 0 : 1);
		if(simModIt->alpha<0) subsec=0;
		if(simModIt->alpha>0) subsec=1;

		int simcell = dddConst.assignCell(x,y,ilay,subsec,false).second;
		simCellH->SetBinContent(xbin,ybin,simcell);
		int reccell = dddConst.assignCell(x/10.,y/10.,recLay,subsec,true).second;
		recCellH->SetBinContent(xbin,ybin,reccell);
		
		std::pair<int,int> ixy=dddConst.findCell(simcell,ilay,subsec,false);
		simixH->SetBinContent(xbin,ybin,ixy.first);
		simiyH->SetBinContent(xbin,ybin,ixy.second);
		
		std::pair<int,int> irecxy=dddConst.findCell(reccell,recLay,subsec,true);
		recixH->SetBinContent(xbin,ybin,irecxy.first);
		reciyH->SetBinContent(xbin,ybin,irecxy.second);
		
		//local coordinates (sim and reco)
		if(reccell>=0 && simcell>=0){
		  std::pair<float,float> simxy=dddConst.locateCell(simcell,ilay,subsec,false);
		  std::pair<float,float> recxy=dddConst.locateCell(reccell,recLay,subsec,true);
		  float dx=(simxy.first-recxy.first*10);
		  if(fabs(dx)<1e-3) dx=0;
		  dxH->SetBinContent(xbin,ybin,dx);
		  float dy=(simxy.second-recxy.second*10);
		  if(fabs(dy)<1e-3) dy=0;
		  dyH->SetBinContent(xbin,ybin,dy);
		}
	      }
	}
    }
  
  testDone_=true;
}




//define this as a plug-in
DEFINE_FWK_MODULE(HGCGeometryAnalyzer);
