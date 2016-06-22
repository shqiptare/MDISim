/* ~/c/MDISim/main_LHC.C        ~/mad/lhc/lhc_V6_5_survey.madx

 rm -rf "/tmp/$LOGNAME/"*

 Unix make
 proj=MDISim ; targ=LHC ; cmakeUnix $proj ;   cd /tmp/$LOGNAME/$proj ; make -j5 ; make install
 proj=MDISim ; targ=LHC ; cmakeUnix $proj ;   cd /tmp/$LOGNAME/$proj ; make $targ -j5 ; mv $targ /tmp/$LOGNAME/Release/. ; cd /tmp/$LOGNAME ; Release/$targ -h

 Xcode
 proj=MDISim ; targ=LHC ; cmakeXcode $proj


 --- 3 alternatives to make  /tmp/$LOGNAME/Release/LHC    for running in line mode ---
 rm -rf "/tmp/$LOGNAME/Release"* ; mkdir /tmp/$LOGNAME/Release ; cd /tmp/$LOGNAME/Release ; cmake ~/c/MDISim      ; cd /tmp/$LOGNAME/Release ; make LHC -j5 ; /tmp/$LOGNAME/Release/LHC
 rm -rf "/tmp/$LOGNAME/Release"* ; mkdir /tmp/$LOGNAME/Release ; cd ~/c/MDISim/ ; make LHC -j5 ; mv -f LHC /tmp/$LOGNAME/Release/. ; make clean ; /tmp/$LOGNAME/Release/LHC
 cd /tmp/$LOGNAME ; xcodebuild -project ~/c/MDISim/MDISim.xcodeproj -target LHC -configuration Release ; Release/LHC -h

 save geometry files in ~/www/Geom  for import from web    TGeoManager::Import("http://hbu.web.cern.ch/hbu/Geom/LHC_IR_5.root");




 --- run in line mode by
 cd /tmp/$LOGNAME  ; rm -f * ; time Release/LHC ; find . -name "*.eps" -exec pstopdf {} \;  ; open *.pdf


 read b1 and b2 optics and twiss files
 combined and plot both beams around ir

 OTHER_CFLAGS = -Wall -fvisibility=default $(MyClib_CFLAG) $(ROOT_CFLAGS) -I$(HOME)/c/MyRoot -I$(HOME)/myroot -DCheckSyntax
 OTHER_LDFLAGS = $(MyClib_LIB) $(MyRoot_LIB) $(ROOT_LIBS)

 cd /tmp/$LOGNAME ; xcodebuild -project ~/c/MDISim/MDISim.xcodeproj -alltargets clean
 cd /tmp/$LOGNAME ; xcodebuild -project ~/c/MDISim/MDISim.xcodeproj -alltargets -configuration Release; echo good test to make all
 cd /tmp/$LOGNAME ; xcodebuild -project ~/c/MDISim/MDISim.xcodeproj -list
 cd /tmp/$LOGNAME ; xcodebuild -project ~/c/MDISim/MDISim.xcodeproj -target LHC -configuration Debug
 cd /tmp/$LOGNAME ; xcodebuild -project ~/c/MDISim/MDISim.xcodeproj -target LHC -configuration Release; /tmp/$LOGNAME/Release/LHC
 cd /tmp/$LOGNAME ; rm -f * ; time Release/LHC ; find . -name "*.eps" -exec pstopdf {} \;  ; open *.pdf


 */

#ifndef __CLING__
#include "MyCppStdHeaders.H"
#include <TRint.h>
#include <TCanvas.h>
#include <TGeoManager.h>
#endif

#include "myEve.h" // StartEveWithGeomDisplay
#include "Apertures.h"
#include "IPpos.h"
//#include "Lorentz.h"
//#include "Twiss_util.h" // Beam
#include "myGraph.h"    // MyApertureGraph
#include "myGraph3D.h"  // PlotTracksInFiles
#include "MyNtuple2Geom.h" // GeomFromMyNtuple

void FindIPs(const Ntuple &nt_b1,const Ntuple &nt_b2,unsigned int verbose) // for the moment only for info
{
  if(verbose) cout << "FindIPs start verbose=" << verbose <<'\n';
  IPpos *ips_b1=new IPpos(&nt_b1,verbose); ips_b1->Print(); // in Apertures.C   CINT  has problems with variable vector size
  IPpos *ips_b2=new IPpos(&nt_b2,verbose); ips_b2->Print();
  CalcAveBeam(nt_b1,nt_b2,ips_b1,ips_b2,verbose); // construct the average beam from b1, b2, in  IPpos.C
}

void Plot_LHC(Ntuple& nt_b1,Ntuple& nt_b2,bool CleanUpApertures,unsigned int verbose)
{
  if(verbose) cout << "Plot_LHC start CleanUpApertures=" << CleanUpApertures << " verbose=" << verbose << endl;

  Apertures* LHC_aper_b1 = new Apertures(&nt_b1,verbose);
  Apertures* LHC_aper_b2 = new Apertures(&nt_b2,verbose);
  
  const double smax=350; // 320 m just start of arc,

  if(CleanUpApertures)
  {
    // work on Apertures of b1
    if(verbose>1) LHC_aper_b1->Print();
    LHC_aper_b1->CheckConsistence();
    LHC_aper_b1->CleanUpApertures();
    if(verbose>1) LHC_aper_b1->Print();
    LHC_aper_b1->CombineLeftRight(smax);
    LHC_aper_b1->CheckConsistence();
    if(verbose) LHC_aper_b1->Print(0,1000,1000);
    Ntuple* nt_cleaned_aper_b1=LHC_aper_b1->MakeNtuple();
    nt_cleaned_aper_b1->PrintSummary();
    nt_cleaned_aper_b1->WriteAsciiNtFile("","LHC_aper_cleaned_b1.out");

    // work on Apertures of b2
    if(verbose>1) LHC_aper_b2->Print();
    LHC_aper_b2->CheckConsistence();
    LHC_aper_b2->CleanUpApertures();
    if(verbose>1) LHC_aper_b2->Print();
    LHC_aper_b2->CombineLeftRight(smax);
    LHC_aper_b2->CheckConsistence();
    if(verbose) LHC_aper_b2->Print(0,1000,1000);
    Ntuple* nt_cleaned_aper_b2=LHC_aper_b2->MakeNtuple();
    nt_cleaned_aper_b2->PrintSummary();
    nt_cleaned_aper_b2->WriteAsciiNtFile("","LHC_aper_cleaned_b2.out");
  }
  if(verbose>2) cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << __LINE__ << endl;

  // double smax= 50, smin=-smax;  double xmin=-0.20,xmax=0.20;
  double smin=-smax;  double xmin=-0.20,xmax=0.20;
  // double smax=6000; double xmin=-8500,xmax=100;
  LHC_aper_b1->PlotFrame(smin,smax,xmin,xmax); // plot frame and keep limits in Apertures class
  // LHC_aper_b2->PlotFrame(smin,smax,xmin,xmax);
  // LHC_aper_b1->PlotFrame(    0,smax,-xymax,xymax); // only right of IP
  cout << "Plot_LHC LHC_aper_b1->zvec.size()=" << LHC_aper_b1->zvec.size() << " call MyApertureGraph" << '\n';
  MyApertureGraph(LHC_aper_b1->zvec.size(),kBlue,LHC_aper_b1->zvec,LHC_aper_b1->x_2_vec,LHC_aper_b1->aper_1_vec);
  MyApertureGraph(LHC_aper_b2->zvec.size(),kRed ,LHC_aper_b2->zvec,LHC_aper_b2->x_2_vec,LHC_aper_b2->aper_1_vec);

  TCanvas* c2=LHC_aper_b1->GetCanvas();
  c2->SaveAs("MyApertureGraph_2d.pdf");

#ifndef __CLING__
  // c2->Close();
#endif
  if(verbose) cout << "Plot_LHC done" <<'\n';
}


#ifndef __CLING__
int main(int argc,char * const argv[]) // normally this
{
  unsigned int verbose=1;
#if __cplusplus >= 201103L // C++11 on
  if(verbose) cout << "main compiled with C++11 on" << '\n';
#else
  if(verbose) cout << "main compiled with C++11 off" << '\n';
#endif

#ifdef __CLING__
  if(verbose) cout << "__CLING__=" << __CLING__ << '\n';
#else
  if(verbose) cout << "__CLING__ not defined" << '\n';
#endif

#ifdef __CINT__
  if(verbose) cout << "__CINT__=" << __CINT__ << '\n';
#else
  if(verbose) cout << "__CINT__ not defined" << '\n';
  TRint* theApp = new TRint("",0,NULL);
#endif


  bool AddTrackToGeometry=true; // CSPE -- better extra program   AddTracksToGeom
  if(AddTrackToGeometry)
  {
    unsigned int verbose=1;
    // StartEveWithGeomDisplay("LHC_IR_5.root");  // local
    StartEveWithGeomDisplay("http://hbu.web.cern.ch/hbu/Geom/LHC_IR_5.root");  // or local  local LHC_IR_5.root  or  saved /Users/hbu/www/Geom/LHC_IR_5.root  or  http://hbu.web.cern.ch/hbu/Geom/LHC_IR_5.root
    const double zmin=0,zmax=330,Scale_xy=100;
    verbose=0;PlotTracksInFiles("~/mad/lhc/twiss_track_ref/0*_px*b1.tfs","~/mad/lhc/output/lhcb1_survey_from5.tfs",zmin,zmax,kGreen,Scale_xy,verbose); // myGraph3D.C
    verbose=0;PlotTracksInFiles("~/mad/lhc/twiss_track_ref/0*_px*b2.tfs","~/mad/lhc/output/lhcb2_survey_from5.tfs",zmin,zmax,kGreen,Scale_xy,verbose);

    theApp->Run();
    exit(0);
  }

  bool CleanUpApertures=true;
  verbose=0;
  Ntuple *ntptr_b1=NULL,*ntptr_b2=NULL;
  if(CleanUpApertures)
  {
    ntptr_b1=new Ntuple(ReadAndMerge("~/mad/lhc/output/lhcb1_twiss_from5.tfs","~/mad/lhc/output/lhcb1_survey_from5.tfs",verbose));
    ntptr_b1->RemoveVarList({"GLOBALTILT","SLOT_ID","ASSEMBLY_ID"});
    if(verbose) ntptr_b1->PrintSummary();
    ntptr_b2=new Ntuple(ReadAndMerge("~/mad/lhc/output/lhcb2_twiss_from5.tfs","~/mad/lhc/output/lhcb2_survey_from5.tfs",verbose));
    ntptr_b2->RemoveVarList({"GLOBALTILT","SLOT_ID","ASSEMBLY_ID"});
    if(verbose) ntptr_b2->PrintSummary();
  }
  else // use existing cleaned aperture files
  {
    Ntuple* ntptr_b1=new Ntuple("~/mad/lhc/output/LHC_aper_cleaned_b1.out"); if(verbose) ntptr_b1->PrintSummary();
    Ntuple* ntptr_b2=new Ntuple("~/mad/lhc/output/LHC_aper_cleaned_b2.out"); if(verbose) ntptr_b2->PrintSummary();
    if(verbose) FindIPs(*ntptr_b1,*ntptr_b2,verbose);
  }

  bool Plot3D=true;
  if(Plot3D)
  {
    const double zmin=0,zmax=320,Scale_xy=100;
    // const double zmin=-320,zmax=320;
    // const double zmin=-15,zmax=15;
    //const double zmin=-180,zmax=180;
    verbose=0; GeomFromMyNtuple("LHC_IR_5","IP",*ntptr_b1,zmin,zmax,kBlue,Scale_xy,verbose);
    verbose=0; GeomFromMyNtuple("LHC_IR_5","IP",*ntptr_b2,zmin,zmax,kRed ,Scale_xy,verbose); // use same name  LHC_IR_5   to draw together
    gGeoManager->Export("LHC_IR_5.root"); // export the geometry

    if(gGeoManager==NULL) TGeoManager::Import("LHC_IR_5.root");
    gGeoManager->GetTopVolume()->Draw("ogl"); // "ogl" drawing does not allow for track and interferes with drawing in c1  -- so better keep separate

    // Int_t PdgCode=2212; // Proton
    // mygraph3d->verbose=2; mygraph3d->MyTrack3D("LHC_IR_5",ntptr_b1,PdgCode);

    // StartEveWithGeomDisplay();

  }
  else // 2dim plot
  {
    if(verbose) FindIPs(*ntptr_b1,*ntptr_b2,verbose);
    Plot_LHC(*ntptr_b1,*ntptr_b2,CleanUpApertures,verbose); // in main_LHC.C
  }
#ifdef __CINT__
  cout << " running with CINT"<< '\n';
#else
  cout << " running in line mode. When done quit root from the menu bar / File " << '\n';
  theApp->Run();
#endif
};
#endif
