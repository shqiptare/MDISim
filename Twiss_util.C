/*   ~/c/MDISim/Twiss_util.C    H. Burkhardt 5/2013     more general 8/2013, to be also useful for LEP, TLEP with synchrotron radiation
 
 further generalized and combined with machine, using Ntuple rather than NtupleRoot, 11/11/2014   -   could go to MyClib at some point
 
 also called by   ~/mad/LEP/lep.madx    was ~/mad/fcc-ee/toolkit/Twiss_util.C
 
 
 
 test syntax by compilation
 cd ~/c/MDISim ; clang -v -Wall $MyClib_CFLAG -I$ROOT_CFLAGS -I$HOME/myroot  -DCheckSyntax -c Twiss_util.C -o /dev/null
 
 root examples see also ~/c/MyRoot/main_root.C   ~/c/MadRoot/main_TwissTrack_tfs.C  ~/mad/synrad/SingleBend.madx
 
 
 next did not work in root
 .include $(HOME)/c/MyClib
 .include $HOME/c/MyClib
 gSystem->AddIncludePath(" -I$HOME/c/MyClib ")
 
 ok in root
 .include /Users/hbu/c/MyClib   for root5, not needed for root6
 
 
 
 
 .include   shows   -I/opt/local/include/root -I$(HOME)/c/MyClib       so   MyClib should be there already
 
 
 
 cd /tmp/$LOGNAME ; rm -f "/tmp/$LOGNAME/"* ; time madx < ~/mad/LEP/lep.madx >& lep.out ; open lep.out ; ps2pdf madxplot.ps ; open madxplot.pdf
 
 cd /tmp/$LOGNAME ; ln -sf ~/myroot/rootlogon.C; ln -sf ~/c/MyRoot/MyRootConst.h ; time root
 .L ~/c/MDISim/Twiss_util.C
 Beam* LEP = new Beam("~/mad/LEP/lep_twiss.tfs",2);
 // Beam* LEP = new Beam("/tmp/$LOGNAME/lep_twiss.tfs",2);
 LEP->Print();
 .q
 
 
 
 cd /tmp/$LOGNAME ; ln -sf ~/myroot/rootlogon.C; ln -sf ~/c/MyRoot/MyRootConst.h ; time root
 .L ~/c/MDISim/Twiss_util.C
 Beam* FCC_ee= new Beam("/tmp/$LOGNAME/tlep_twiss.tfs",2);
 FCC_ee->Print();
 .q
 
 
 cd /tmp/$LOGNAME ; ln -sf ~/myroot/rootlogon.C; ln -sf ~/c/MyRoot/MyRootConst.h ; time root
 .L ~/c/MDISim/Twiss_util.C
 Beam* FCC_ee_060= new Beam("/tmp/$LOGNAME/tlep_7hb_060.tfs",0);
 Beam* FCC_ee_120= new Beam("/tmp/$LOGNAME/tlep_7hb_120.tfs",0);
 FCC_ee_060->Print();
 FCC_ee_120->Print();
 .q
 
 
 cd /tmp/$LOGNAME ; ln -sf ~/myroot/rootlogon.C; ln -sf ~/c/MyRoot/MyRootConst.h ; time root
 .L ~/c/MDISim/Twiss_util.C
 Beam* FCC_ee= new Beam("~/mad/fcc-ee/BernhardHolzer_2013_10_11/tlep_175_twiss.tfs",2);
 FCC_ee->Get_Part().PrintLorentz(cout);
 FCC_ee->Print();
 .q
 
 cd /tmp/$LOGNAME ; ln -sf ~/myroot/rootlogon.C; ln -sf ~/c/MyRoot/MyRootConst.h ; time root
 .L ~/c/MDISim/Twiss_util.C
 Beam* fcc_hh_race= new Beam("~/mad/fcc_hh/output/fcc_hh_race_v3_twiss.tfs",2); // test for hadron machine
 fcc_hh_race->Get_Part().PrintLorentz(cout);
 fcc_hh_race->Print();
 .q
 
 just needs alternative way to define rho  (or I_2)
 
 
 */


#ifndef __CLING__
#include "MyCppStdHeaders.H"
#endif

// #include "MyRootConst.h"
// #include "Ntuple.h"
#include "Lorentz.h"
#include "Mat6x6.h"
#include "Twiss_util.h"
#include "SynRad.h"

Beam::Beam(const string FnameStr,const unsigned int verbose) : RFHV(0), Harm(0), Cs(Echarg/(3*eps0)) // constructor, reads Energy and synchrotron radiation integrals from the tfs header
{
  const char* Fname=FnameStr.c_str();
  this->verbose=verbose;
  if(verbose) cout << "Beam constructor start verbose=" << verbose << '\n';
  const char* FullFileName=ResolveHome(Fname);
  if(verbose) cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << __LINE__ << endl;
  const char* PartName=GetParamString_tfs("PARTICLE",FullFileName).c_str();
  if(verbose) cout << "from twiss PartName=" << PartName << '\n';
  Eb=GetParamValue_tfs("ENERGY",Fname);
  double mass  =GetParamValue_tfs("MASS",Fname);
  Part=NULL;
  phis=0;Qs=0;sige=0;sigz=0; // init with 0
  if(fabs(mass-MyConstPhys::Mele)<1.e-5) // electron
  {
    Part=new Lorentz(Eb,MyConstPhys::Mele,"e");
  }
  else if(fabs(mass-MyConstPhys::Mproton)<1.e-5) // proton
  {
    Part=new Lorentz(Eb,MyConstPhys::Mproton,"p");
  }
  else
  {
    Part=new Lorentz(Eb,mass,PartName); // by default take name from twiss
    if(verbose) cout << "Part->Name=" << Part->Name << '\n';
  }
  if(verbose) Part->PrintLorentz(cout);
  cqm=55.*Hbarc/(32.*sqrt3*Part->m); // calculate cq with the actual particle mass
  gamma=GetParamValue_tfs("GAMMA",Fname);
  I_1=GetParamValue_tfs("SYNCH_1",Fname);
  I_2=GetParamValue_tfs("SYNCH_2",Fname);
  I_3=GetParamValue_tfs("SYNCH_3",Fname);
  I_4=GetParamValue_tfs("SYNCH_4",Fname);
  I_5=GetParamValue_tfs("SYNCH_5",Fname);
  if(I_2==0)
  {
    cout << " warning synchrotron radiation integrals not defined " << '\n'; //-- calculate myself, write a utility in Twiss_util.C
  }
  else
  {
    rho=2*M_PI/I_2;
  }
  ex    =GetParamValue_tfs("EX",  Fname);    // geometrical emittance
  ey    =GetParamValue_tfs("EY",  Fname);    // geometrical emittance
  exN=Part->EmitN(ex);
  eyN=Part->EmitN(ey);
  alfac =GetParamValue_tfs("ALFA",  Fname);  // momentum compaction
  kbunch=GetParamValue_tfs("KBUNCH",Fname);  // number of bunches / beam
  length=GetParamValue_tfs("LENGTH",Fname);  // total length or circumference of the machine
  npart =GetParamValue_tfs("NPART", Fname);  // number of particles / bunch
  eta =  alfac-1./(gamma*gamma); // slip factor
  J[0]=J[1]=1;J[2]=2;
  CalcDerived();
  if(verbose) cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << __LINE__ << endl;
}

void Beam::CalcDerived()
{
  // derived quantities
  // fRF=Harm*frev;   length * frev =  v = beta * CLight;
  frev = Part->beta*CLight/length;
  ibun = npart * frev * fabs(Part->Z)*Echarg;
  U0=1.e-9* Cs *pow(gamma,4.)/rho;
  if(verbose) cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << __LINE__ << endl;
}

void Beam::Print()
{
  if(verbose) cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << __LINE__ << endl;
  cout << "Eb=" << Eb << " GeV" << '\n';
  cout << "gamma=" << gamma << '\n';
  cout << "alfac=" << alfac << '\n'; // momentum compaction
  cout << "length=" << length << " m   alfac*length=" << alfac*length << " m" << '\n';
  
  cout << "kbunch=" << kbunch << '\n';
  cout << "npart=" << npart << '\n';
  
  cout << "I_1=" << I_1 << " m" << '\n';
  cout << "I_2=" << I_2 << " 1/m rho=" << rho << '\n';
  cout << "I_3=" << I_3 << " 1/m2" << '\n';
  cout << "I_4=" << I_4 << '\n';
  cout << "I_5=" << I_5 << " 1/m" << '\n';
  
  cout << "Echarg=" << Echarg << '\n';
  cout << "Cs=" << Cs << '\n';
  cout << "cqm=" << cqm << '\n';
  cout << "U0=" << U0 << " GeV" << '\n';
  
  double sige=0;
  if(2*I_2+I_4 >0) sige=gamma*sqrt(cqm*I_3/(2*I_2+I_4));
  cout << "sige=" << sige << '\n';
  cout << " damping partition numbers Jx, Jy, Jz =" << J[0] << " " << J[1] << " " << J[2] << '\n';
  cout << "ex=" << 1.e9*ex << " ey=" << 1.e9*ey << " nm  geometrical emittances, normalized exN=" << exN << " eyN=" << eyN << " m" << '\n';
  
  cout << "RFHV=" << GVoltsWithUnits(RFHV) << " Harm=" << Harm << '\n';
  cout << " Qs=" << Qs  << " sige=" << sige  << " sigz=" << LengthWithUnits(sigz) << " sigb=" << sigb << '\n';
  
  cout << "frev=" << FreqWithUnits(frev) << " fRF=" << FreqWithUnits(Harm*frev) << '\n';
  ibeam=kbunch*ibun;
  cout << "ibun=" << CurrWithUnits(ibun) << " ibeam=" << CurrWithUnits(ibeam) << '\n';
  cout << "SR Power / beam =" << PowerWithUnits(ibeam*1.e9*U0) << '\n';
};

void Beam::Set_Eb(double Eb)
{
  this->Eb=Eb;
  Part->E=Eb;
  Part->FromE();
}

void Beam::EmitFromSynrad(double kappa) // calculate emittances and bunch length from synchrotron radiation and ey=kappa*ex
{
  double Jx=J[0];
  if(I_2*Jx>0)
  {
    ex=this->cqm*gamma*gamma*I_5/(I_2*Jx); // calculate the equilibrium emittance from synchrotron radiation
    ey=kappa*ex;
    exN=Part->EmitN(ex);
    eyN=Part->EmitN(ey);
    // calculate also the energy spread and bunch length
    sige=gamma*sqrt(cqm*I_3/(2*I_2+I_4));
    const double U0thres=1.e-10; // below use stationary bucket
    if(U0<U0thres) // energy loss negligible
    {
      if(eta<0) phis=0; else phis=M_PI; // stationary bucket
    }
    else // non negligible energy loss
    {
      double sinphis=U0/(Part->Z*RFHV);
      phis=M_PI-asin(sinphis);
    }
    Qs=Part->QsFromRFPhaseRFHV(phis,alfac,RFHV*Harm);
    sigz=alfac*sige*length/(Qs*2*M_PI); // as beam  Machine::CalcDampPartSigZTauq
    
    if(verbose>1) cout << "RFHV=" << GVoltsWithUnits(RFHV) << " alfac=" << alfac << " sige=" << sige << " length=" << LengthWithUnits(length) << " phis=" << phis << " Qs=" << Qs << " sigz=" << sigz << '\n';
    
    double fv=-cos(phis)+(M_PI_2-phis)*sin(phis);  // M_PI_2 is " pi / 2 " defind in /usr/include/math.h
    double sigbMax2=2*pow(Part->beta,2)*RFHV/(M_PI*abs(eta)*Harm*Part->E);
    sigb=sqrt(sigbMax2*fv); //bucket half-height
    if(verbose>1) cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << __LINE__
      << " ex=" << LengthWithUnits(ex)
      << " ey=" << LengthWithUnits(ey)
      << " Qs=" << Qs
      << " phis=" << phis
      << " sige=" << sige
      << " sigz=" << LengthWithUnits(sigz)
      << " sigb=" << sigb
      << '\n';
  };
}

void Beam::SetEmitN(double exN,double eyN) // normally from tfs header or synrad. Here set explicitely, maybe useful for linacs and protons
{
  this->exN=exN; // set normalized emittances given on input
  this->eyN=eyN;
  ex=Part->Emit(exN); // calculate geometrical emittance from normalized
  ey=Part->Emit(eyN);
}

void Beam::RFHV_from_twiss(const Ntuple& nt) // also sets Harm
{
  if(verbose) cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << __LINE__ << endl;
  // RF volts
  valarray<double> RF_volts=nt.GetVar("VOLT");
  if(verbose) cout << " RF_volts.size()=" << RF_volts.size() << '\n';
  
  RFHV=0;
  for(unsigned int i=0;i<RF_volts.size();++i) RFHV+=RF_volts[i];
  
  RFHV*=0.001; // from MV to GV
  
  if(nt.VarExists("HARMON"))
  {
    Harm=nt.Max("HARMON");
  }
  
  if(Harm==0)
  {
    if(nt.VarExists("FREQ"))
    {
      double fRF=nt.Max("FREQ");
      if(fRF>0) Harm= (int) ((fRF*1.e6/frev)+0.5);
      if(verbose) cout << "HARMON not given, derive from rf fequency fRF=" << fRF << " MHz and frev=" << frev << " Hz Harm=" << Harm << '\n';
    }
  }
}

void Beam::Touschek(Ntuple& nt,unsigned int verbose) //hbu started 21/01/2015
// see also BeamParam::Touschek in ~/c/beam/BeamParm.C
// here rough estimate, more accurate use Piwinski formulas, directly Mad-X..
{
  if(verbose) cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << __LINE__ << '\n';
  valarray<double> betx=nt.GetVar("BETX");
  valarray<double> bety=nt.GetVar("BETY");
  valarray<double> dx=nt.GetVar("DX");
  
  unsigned int n=nt.Noent();
  
  double sigxd=0;
  for(unsigned int i=0;i<n;++i) sigxd+=sqrt( betx[i]*ex + sige*sige*dx[i]*dx[i]); // beam size with dispersion
  
  // calculate the average betas
  double AveBetx=0,AveBety=0;
  for(unsigned int i=0;i<n;++i) AveBetx+=betx[i];
  for(unsigned int i=0;i<n;++i) AveBety+=bety[i];
  AveBetx=AveBetx/n;
  AveBety=AveBety/n;
  cout << "Touschek AveBetx,y=" << AveBetx << " " << AveBety << " m" << '\n';
  // calculate the average betas from inverse betas;
  double AveInvBetx=0,AveInvBety=0;
  for(unsigned int i=0;i<n;++i) AveInvBetx+=1./betx[i];
  for(unsigned int i=0;i<n;++i) AveInvBety+=1./bety[i];
  AveInvBetx=n/AveInvBetx;
  AveInvBety=n/AveInvBety;
  cout << "Touschek 1/AveInvBetx,y=" << AveInvBetx << " " << AveInvBety << " m" << '\n';
  
  double sigx=sqrt(AveInvBetx*ex);
  double sigy=sqrt(AveInvBety*ey);
  //double sigxp\sqrt(); // with disperion
  
  
  cout << "Touschek ex, ey=" << LengthWithUnits(ex) << " " << LengthWithUnits(ey) << '\n';
  cout << "Touschek average sigx,y=" << LengthWithUnits(sigx) << " " << LengthWithUnits(sigy) << '\n';
  
  
  // flat beam, as used by Boscolo Raimondi
  valarray<double> sigxp=nt.GetVar("SIGXP");
  double Ave_sigxp=0;
  for(unsigned int i=0;i<n;++i) Ave_sigxp+=sigxp[i];
  Ave_sigxp=Ave_sigxp/n;
  
  double eps=sigb;
  double umin=eps/(Part->gamma*Ave_sigxp);umin=umin*umin;
  double Volume=sigx*sigy*sigz;
  
  double TauInv= sqrt(M_PI)*Re*Re*CLight*npart/( pow(4*M_PI,1.4)*pow(Part->gamma,3) * Volume*Ave_sigxp*eps*eps)  *F_Touschek_appr(umin);
  cout << "Touschek Volume=" << Volume << " m^3 npart=" << npart << " umin=" << umin << " eps=sigb=" << eps
  << " tau=" << TimeWithUnits(1./TauInv) << " = " << 1./(3600*TauInv) << "h" << '\n';
  
  /*
   // from  bs, coh
   double pmax=sigb*Part->P; // longitudinal momentum acceptance - take the bucket height
   double Volume_bs=pow((double)(4.*M_PI),1.5)*sigx*sigy*sigz; // eff. volume of gaussian, in m**3, see  bs, coh
   double TauInv_bs=8.*M_PI*Re*Re*CLight*npart/Volume_bs * pow(Part->m/pmax,2) *F_Touschek_appr(umin);
   cout << "Touschek Volume=" << TauInv_bs << " m^3 npart=" << npart
   << " tau_bs=" << 1./(3600*TauInv_bs) << "h" << '\n';
   */
}

void ShiftToElement(Ntuple &nt,string ElementName,unsigned int verbose) // shift IP to 0 in s and z
{
  vector<unsigned int> i_IP_List=FindElements(nt,"IP",verbose,"NAME"); // all elements of type S or RBEND -  here all SBEND, there are 1783, of which 1640 are standard bends
  if(i_IP_List.size()>0)
  {
    unsigned int is=nt.ivar("S");
    unsigned int iz=nt.ivar("Z");
    unsigned int i=0; // take the first IP
    unsigned int ient=i_IP_List[i];
    double s_IP =nt.Entry(ient,is);
    double z_IP =nt.Entry(ient,iz);
    if(verbose)
    {
      cout << "i=" << i << " ient=" << ient << " s_IP=" << s_IP << " z_IP=" << z_IP << '\n';
    }
    nt.AddToVar("S",-s_IP);
    nt.AddToVar("Z",-z_IP);
  }
}

void CalcSynrad(Ntuple &nt,const Beam &Mach,unsigned int verbose,const double nsigx,const double nsigy) // actual calculation still mostly in PrintBends, PrintQuads
{
  if(verbose>2) cout << CPUTime() << '\n';
  if(verbose>1) cout << __FILE__ << " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " verbose=" << verbose << " Mach.Part->E=" << Mach.Part->E << " GeV ex=" << 1.e9*Mach.ex << " nm  nsigx,y=" << nsigx << ", " << nsigx << '\n';
  
  //-- calculate for all elements and add to ntuple
  unsigned int n=nt.Noent();
  vector<string> KeywordCol=nt.GetStrCol("KEYWORD");
  valarray<double> EcritBend(n);
  valarray<double> ngamBend(n);
  valarray<double> L  =nt.GetVar("L");
  valarray<double> Angle  =nt.GetVar("ANGLE");
  
  if(nt.VarExists("LRAD"))
  {
    if(verbose) cout << "LRAD exists. Will be used as length if L<=0 and LRAD>0" << '\n';
    valarray<double> LRAD =nt.GetVar("LRAD");
    for(unsigned int i=0;i<LRAD.size();++i) if(L[i]<=0 && LRAD[i]>0) L[i]=LRAD[i];
  }
  nt.SetVar("L",L); // keep the modified
  
  // use orbit correctors like bends
  valarray<double> hkick,vkick;
  if(nt.VarExists("HKICK")) hkick=nt.GetVar("HKICK");
  if(nt.VarExists("VKICK")) vkick=nt.GetVar("VKICK");
  unsigned int iname=nt.ivar("NAME");
  for(unsigned int ient=0;ient<n;++ient)
  {
    if (KeywordCol[ient]==string("HKICKER") && ient<hkick.size())
    {
      if(hkick[ient]!=0)
      {
        if(verbose) cout << "iele " << setw(5) << ient << " " << setw(15) << nt.StrEntry(ient,iname) << " is horizontal corrector with hkick=" << setw(10) << hkick[ient] << " will be used as bending angle L,LRAD=" << L[ient] << " m" << '\n';
        Angle[ient]=hkick[ient];
      }
    }
    else if (KeywordCol[ient]==string("VKICKER") && ient<hkick.size())
    {
      if(vkick[ient]!=0)
      {
        if(verbose) cout << "iele " << setw(5) << ient << " " << setw(15) << nt.StrEntry(ient,iname) << " is vertical   corrector with vkick=" << setw(10)<< vkick[ient] << " will be used as bending angle L,LRAD=" << L[ient] << " m" << '\n';
        Angle[ient]=vkick[ient];
      }
    }
  }
  nt.SetVar("ANGLE",Angle); // keep the modified
  
  for(unsigned int i=0;i<n;++i)
  {
    if(Angle[i]!=0 && L[i]!=0) // rbend, sbend, orbit corrector
    {
      double k0=Angle[i]/L[i]; // 1/rho
      EcritBend[i]=1.5*Hbarc*pow(Mach.Part->gamma,3)*fabs(k0);
      ngamBend[i]=NgamBend(fabs(Angle[i]),Mach.Part->gamma); // number of photons radiated in a bend with Angle
    }
    //--   could also add here quadrupole calculation - still done in fly - or do both on the fly without adding columns ?
  }
  nt.AddVar("EcritBend",EcritBend); if(verbose>2) cout << "after AddVar EcritBend" << endl;
  nt.AddVar("ngamBend", ngamBend);  if(verbose>2) cout << "after AddVar ngamBend" << endl;
  
  vector<unsigned int> i_Rbend_List=FindElements(nt,"BEND",verbose,"KEYWORD"); // all elements of type S or RBEND -  here all SBEND, there are 1783, of which 1640 are standard bends
  if(verbose>1) cout << "i_Rbend_List.size=" << i_Rbend_List.size() << '\n';
  
  // include also orbit correctors
  vector<unsigned int> i_hvcorr_List=FindElements(nt,"KICKER",verbose,"KEYWORD"); // orbit correctors called somewhat misleading HKICKER, VKICKER in MAD-X
  if(verbose>1) cout << "i_hvcorr_List.size=" << i_hvcorr_List.size() << '\n';
  
  // combine bends and correctors into new combined list
  vector<unsigned int> i_bend_corr_list;
  vector<unsigned int>::iterator ic=i_hvcorr_List.begin();
  for(auto ib=i_Rbend_List.begin();ib<i_Rbend_List.end();++ib)
  {
    if(verbose>2) cout << "  *ib=" << setw(5) << *ib << " *ic=" << setw(5) << *ic << '\n';
    while(ic<i_hvcorr_List.end() && *ic < *ib)
    {
      i_bend_corr_list.push_back(*ic);
      ++ic;
    }
    i_bend_corr_list.push_back(*ib);
  }
  if(verbose>1) cout << "i_bend_corr_list.size=" << i_bend_corr_list.size() << '\n';

  vector<unsigned int> i_quad_List=FindElements(nt,"QUADRUPOLE",0,"KEYWORD");
  if(verbose>1) cout << "i_quad_List.size()=" << i_quad_List.size() << '\n';
  // PrintBends(nt, i_Rbend_List,Mach,verbose); // rbend, sbend
  // PrintBends(nt,i_hvcorr_List,Mach,verbose); // correctors "H,VKICKER"
  PrintBends(nt,i_bend_corr_list,Mach,verbose); // both bends and correctors
  PrintQuads(nt,  i_quad_List,Mach,verbose,nsigx,nsigy);
}

vector<unsigned int> FindElements(Ntuple &nt,const string& MatchStr,unsigned int verbose,const string& ColumnName,size_t MatchPos)
{
  if(verbose>1) cout << "start FindElements &nt=" << &nt << " MatchStr=" << MatchStr << " verbose=" << verbose << " ColumnName=" << ColumnName << '\n';
  vector<string> nt_StrCol=nt.GetStrCol(ColumnName);
  if(verbose>1)
  {
    cout << "FindElement matching MatchStr=" << MatchStr
    << " nt_StrCol[0]=" << nt_StrCol[0]
    << " nt_StrCol.size()=" << nt_StrCol.size()
    << '\n';
  }
  unsigned int is=nt.ivar("S");
  unsigned int ibetx=nt.ivar("BETX");
  unsigned int ibety=nt.ivar("BETY");
  unsigned int nfound=0;
  vector<unsigned int> iElemList;
  unsigned int PriMax=20;
  if(verbose>1) PriMax=100;
  if(verbose>2) PriMax=999999; //all
  for(unsigned int ient=0;ient<nt.Noent();++ient)
  {
    size_t MatchStrPos=nt_StrCol[ient].find(MatchStr); // look for MatchStr in nt_StrCol
    if(MatchStrPos!=string::npos)
    {
      if(MatchPos!=string::npos && MatchStrPos!=MatchPos ) // require that the MatchStr was found in the position given as MatchPos
      {
        if(verbose>1) cout << __FILE__ << " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " MatchStr=" << MatchStr << " found but position=" << MatchStrPos << " not the requested MatchStrPos=" << MatchStrPos << '\n';
      }
      else
      {
        nfound++;
        iElemList.push_back(ient);
        double s   =nt.Entry(ient,is);
        double betx=nt.Entry(ient,ibetx);
        double bety=nt.Entry(ient,ibety);
        if(verbose && ient<PriMax) cout << "FindElement Element " << setw(5) << ient << " is " << setw(18) << nt_StrCol[ient] << " at s=" << setw(8) << s << " m"
          << " BETX=" << setw(8) << betx
          << " BETY=" << setw(8) << bety
          << " m"
          << '\n';
      }
    }
  }
  if(verbose>1)
  {
    cout << "FindElement end for MatchStr=" << MatchStr << " nfound=" << nfound <<  " matches" << '\n';
  }
  return iElemList;
}

void PrintQuads(Ntuple& nt, const vector<unsigned int>& elemList, const Beam& Mach, unsigned int verbose,const double nsigx,const double nsigy)
// nsig=0 to get sawtooth only, use 1 to get an estimate of average loss from beam size, use 5 to get extreme case
// calculate and print beam size from beta given in elemList for the emit
// rather generic  for normal (not tilted) quadrupole
// here written for quads, with synchrotron radiation
{
  const vector<string> VarList={"NAME","S","SIGX","SIGY","X","Y","L","K1L"}; if(!nt.VarsExist(VarList)) return;
  
  bool Usekw=true;
  unsigned int isigx =nt.ivar("SIGX");
  unsigned int isigy =nt.ivar("SIGY");
  unsigned int ix=nt.ivar("X"); // x-orbit, from offset,sawtooth
  unsigned int iy=nt.ivar("Y"); // y-orbit, in case of LHC, FCC-hh can be separation for crossing angle bumps
  unsigned int iL=nt.ivar("L");
  unsigned int iK1L=nt.ivar("K1L");
  // for info
  unsigned int idivx =nt.ivar("DIVX");
  unsigned int idivy =nt.ivar("DIVY");
  
  
  if(verbose)
  {
    cout << "Start of PrintQuads nsigx=" << nsigx << " nsigy=" << nsigy << " ";
    Mach.Part->PrintLorentz(cout);
    if(nsigx==0 && nsigy==0) cout << " calculate energy loss from offsets only";
    else if(nsigx==1 && nsigy==1) cout << " default sizes, calculate energy loss from offset and average beam loss from beam size for gaussian beams";
    cout << '\n';
  }
  unsigned int is=nt.ivar("S");
  
  unsigned int PriMax=20;
  if(verbose>1) PriMax=100;
  if(PriMax>elemList.size()) PriMax=elemList.size();
  
  size_t MaxNameLen=0;
  if(verbose>2) cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << __LINE__ << " PriMax=" << PriMax << " elemList.size()=" << elemList.size() << endl;
  
  vector<string> NameCol    =nt.GetStrCol("NAME");
  
  for(unsigned int i=0;i<PriMax;++i)
  {
    unsigned int ient=elemList[i];
    if(verbose>2) cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << __LINE__ << " ient=" << ient << " i=" << i << endl;
    if(NameCol[ient].size() >MaxNameLen) MaxNameLen=NameCol[ient].size();
  }
  
  if(verbose)
  {
    cout << "iele "
    << setw(MaxNameLen) << "Element"
    << "           s        L          x          y        K1L         k0      Ecrit       ngam         Bx         By    Power     SIGX     SIGY    Angle     divx     divy" << '\n'
    << setw(MaxNameLen+16) << " "
    << "m        m         mm         mm        m-2        m-1        keV                     T          T     ";
    if(Usekw) cout << "  kW"; else cout << "Watt";
    cout << "       mm       mm     mrad     mrad     mrad";
    cout << '\n';
  }
  double Angle,k0,K1l,Ecrit,EmeanBend,ngamBend;
  cout << std::defaultfloat;
  double U0Sum=0,U0=0,PowBend=0,PowSum=0;
  double Bx=0,By=0,B=0;
  double brho=Mach.Part->P*1.e9/CLight;
  const double eps=1.e-12;
  for(unsigned int i=0;i<elemList.size();++i)
  {
    unsigned int ient=elemList[i];
    if(ient<nt.Noent())
    {
      K1l=nt.Entry(ient,iK1L);
      if(fabs(K1l)>eps) // skip QTs set to 0      calculate equivalent bend radiation add position and nsigma - adding length in quadrature
      {
        double L   =nt.Entry(ient,iL);
        double s   =nt.Entry(ient,is);
        double x=hypot(nsigx*nt.Entry(ient,isigx),nt.Entry(ient,ix)); // average x   for equivalent bend
        double y=hypot(nsigy*nt.Entry(ient,isigy),nt.Entry(ient,iy)); // average y   for equivalent bend
        double K1=K1l/L;
        Bx=y*brho*K1;
        By=x*brho*K1;
        B=sqrt(Bx*Bx+By*By);
        double r=hypot(x,y);
        k0=r*fabs(K1);    // equivalent bend  k0=1/rho
        Angle=k0*L;
        Ecrit=1.5*Hbarc*pow(Mach.Part->gamma,3)*k0;
        EmeanBend=MeanSyn*Ecrit; // equivalent
        ngamBend=NgamBend(Angle,Mach.Part->gamma); // number of photons radiated in a bend with Angle
        U0=ngamBend*EmeanBend;
        U0Sum+=U0;
        PowBend=Mach.ibeam*1.e9*U0;
        PowSum+=PowBend;
        if(verbose && i<PriMax)
        {
          cout << setw(4) << ient << "  "
          << setw(MaxNameLen) << NameCol[ient] << " "
          << setprecision(4) << std::defaultfloat
          << setw(10) << s << " "
          << setw( 8) << L << " "
          << setw(10) << 1.e3*x << " "
          << setw(10) << 1.e3*y << " "
          << setw(10) << K1l << " "
          << setw(10) << k0 << " "
          << setw(10) << 1.e6*Ecrit << " "
          << setw(10) << ngamBend << " " << setprecision(3)
          << setw(10) << Bx << " "
          << setw(10) << By << " ";
          if(Usekw) cout << setw(8) << 1.e-3*PowBend; else cout << setw(8) << PowBend;
          cout
          << setw(9) << 1.e3*nt.Entry(ient,isigx)
          << setw(9) << 1.e3*nt.Entry(ient,isigy)
          << setw(9) << 1.e3*Angle
          << setw(9) << 1.e3*nt.Entry(ient,idivx)
          << setw(9) << 1.e3*nt.Entry(ient,idivy)
          ;
          cout << '\n';
        }
      }
    }
    else
    {
      cout << "PrintQuads *** error *** ient=" << ient << " > Noent=" << nt.Noent() << '\n';
    }
  }
  cout << "from summing up energy loss in the " << elemList.size() << " quads U0Sum=" << U0Sum << " GeV" << '\n';
  cout << "PowSum=" << PowerWithUnits(PowSum) << '\n';
}

void PrintBends(Ntuple& nt, const vector<unsigned int>& elemList,const Beam& Mach,unsigned int verbose) // with synchrotron radiation calculation,  started 8/2013 with ~/mad/LEP/lep.madx
{
  const vector<string> VarList={"S","L","ANGLE"}; if(!nt.VarsExist(VarList)) return;
  
  bool Usekw=true;
  if(verbose>1)
  {
    cout << "PrintBends start elemList.size()=" << elemList.size() << " ";
    Mach.Part->PrintLorentz(cout);
    cout << "Emean=MeanSyn*Ecrit, where MeanSyn=" << MeanSyn << '\n';
  }
  valarray<double> S        =nt.GetVar("S");
  valarray<double> L        =nt.GetVar("L");
  valarray<double> Angle    =nt.GetVar("ANGLE");
  valarray<double> betx     =nt.GetVar("BETX");
  valarray<double> alfx     =nt.GetVar("ALFX");
  valarray<double> sigx     =nt.GetVar("SIGX");
  valarray<double> divx     =nt.GetVar("DIVX");
  valarray<double> EcritBend=nt.GetVar("EcritBend");
  valarray<double> ngamBend =nt.GetVar("ngamBend");
  vector<string> NameCol    =nt.GetStrCol("NAME");
  vector<string> KeywordCol =nt.GetStrCol("KEYWORD");
  cout << std::defaultfloat;
  double U0sum=0,PowSum=0,PowSum250=0;
  unsigned int PriMax=20;
  if(verbose>1) PriMax=100;
  if(PriMax>elemList.size()) PriMax=elemList.size();
  double Ten_MeV=10.e-3;       // 10 MeV
  double Ten_MeV_over_Ecrit=0;
  
  size_t MaxNameLen=0,MaxKeywordLen=string("ELEMENT").length();
  if(verbose>2) cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << __LINE__ << " PriMax=" << PriMax << " elemList.size()=" << elemList.size() << endl;
  
  for(unsigned int i=0;i<PriMax;++i)
  {
    unsigned int ient=elemList[i];
    if(verbose>2) cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << __LINE__ << " ient=" << ient << " i=" << i << endl;
    
    if(Angle[ient]!=0)
    {
      if(NameCol[ient].size()   >MaxNameLen)       MaxNameLen=NameCol[ient].size();
      if(KeywordCol[ient].size()>MaxKeywordLen) MaxKeywordLen=KeywordCol[ient].size();
    }
    
  }
  if(verbose>2) cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << __LINE__ << endl;
  
  if(verbose)
  {
    cout << "iele  "
    << setw(MaxNameLen) << "NAME" << " "
    << setw(MaxKeywordLen) << "KEYWORD" << " "
    << "      S      L      Angle    Ecrit ngamBend      rho       B     BETX      SIGX     divx      Power frac>10MeV ngam*npart Egamtot Emean " << '\n'
    << setw(10+MaxNameLen+MaxKeywordLen) << " "
    << "    m      m       mrad      keV                 m       T      m        mm       mrad     ";
    if(Usekw) cout << "  kW"; else cout << "Watt";
    cout << "                          GeV      keV";
    cout << '\n';
  }
  
  for(unsigned int i=0;i<elemList.size();++i)
  {
    unsigned int ient=elemList[i];
    if(ient<nt.Noent())
    {
      double k0=Angle[ient]/L[ient];
      double rho=0,Bfield=0;
      if(k0!=0)
      {
        rho=1./k0;
        Bfield=BfieldFromRadius(rho,Mach.Part->P); //  ~/c/MyClib/ConstPhys.h    pmom*1.e9/(CLight*Radius)    "Brho = pmom*1.e9/(CLight"
        rho=fabs(rho);
      }
      double EmeanBend=MeanSyn*EcritBend[ient]; // GeV
      Ten_MeV_over_Ecrit=Ten_MeV/EcritBend[ient];
      double U0=ngamBend[ient]*EmeanBend;
      double PowBend=Mach.ibeam*1.e9*U0;
      PowSum+=PowBend;
      if(S[ient]<250) PowSum250+=PowBend;
      if(verbose && i<PriMax && Angle[ient]!=0 )
      {
        cout << setw(4) << ient << "  "
        << setw(MaxNameLen)    << NameCol[ient] << " "
        << setw(MaxKeywordLen) << KeywordCol[ient] << " "
        << setprecision(4) << std::defaultfloat
        << setw(7) << S[ient] << " "
        << setw(6) << L[ient] << " "
        << setprecision(4)
        << setw(10) << 1.e3*Angle[ient] << " "
        << setprecision(3)
        << setw(8) << 1.e6*EcritBend[ient]  // keV
        << setw(8) << ngamBend[ient]
        << setprecision(1) << std::fixed
        << setw(10) << rho
        << setprecision(3) << std::defaultfloat
        << setw(10) << Bfield
        << setw(9) << betx[ient]
        << setw(9) << 1.e3*sigx[ient]
        << setw(9) << 1.e3*divx[ient] << " " << std::defaultfloat;
        if(Usekw) cout << setw(8) << 1.e-3*PowBend; else cout << setw(8) << PowBend;
        cout << " "
        << setw(9) << SyNorm*SynRadInt(Ten_MeV_over_Ecrit)
        << setw(10) << setprecision(3) << ngamBend[ient]*Mach.npart
        << setw(9) << setprecision(3) << EmeanBend*ngamBend[ient]*Mach.npart
        << setw(8) << 1.e6*EmeanBend // keV
        // << setw(9) << Ten_MeV_over_Ecrit
        ;
        cout << '\n';
      }
    }
    else
    {
      cout << "PrintStrength *** error *** ient=" << ient << " > Noent=" << nt.Noent() << '\n';
    }
  }
  
  // calculate energy gain by rf - should better go to extra function
  valarray<double> deltaE(nt.Noent());
  double detot=0;
  valarray<double> VOLT=nt.GetVar("VOLT");
  if( VOLT.size()>0 ) // H E R E improve -  check to avoid crashes when LAG ,... not existing
  {
    valarray<double>  LAG=nt.GetVar("LAG");
    valarray<double>  FREQ=nt.GetVar("FREQ");
    valarray<double>  HARMON=nt.GetVar("HARMON");
    valarray<double>  T=nt.GetVar("T");  // http://mad.web.cern.ch/mad/madx.old/Introduction/tables.html#normal   in meter units
    double dERF=0;
    double lambdaRF=CLight/704.2013e6;
    for(unsigned int ient=0;ient<nt.Noent();++ient)
    {
      double U0=ngamBend[ient]*MeanSyn*EcritBend[ient];
      U0sum+=U0;
      detot-=U0; //  RF gain - energy loss        here negative, as loss to be balanced by RF
      if(VOLT[ient]>0) // rf
      {
        // double dE=VOLT[ient]*sin(2*M_PI*(LAG[ient]-HARMON[ient]*FREQ[ient]*T[ient])); // formula according to http://mad.web.cern.ch/mad/madx.old/Introduction/cavity.html  does not fit, check with code
        double dE=1.e-3* VOLT[ient]*sin(2*M_PI*(LAG[ient]-T[ient]/lambdaRF)); // formula according to http://mad.web.cern.ch/mad/madx.old/Introduction/cavity.html  does not fit, check with code
        detot+=dE;
        dERF+=dE;
        if(verbose>1) cout << " VOLT[" << ient << "]=" << VOLT[ient]
          << " LAG=" << setw(9) << LAG[ient]
          << " FREQ=" << setw(9) << FREQ[ient]
          << " HARMON=" << setw(9) << HARMON[ient]
          << " T=" << setw(9) << T[ient] << " m"
          << " dE=" << setw(9) << dE
          << " lambdaRF=" << setw(9) << lambdaRF
          << " detot=" << setw(9) << detot
          << " dERF=" << setw(9) << dERF
          << '\n';
      }
    }
    cout << "from summing up energy loss the " << elemList.size() << " bends U0sum=" << U0sum << " GeV   compensated by RF to detot=" << detot << " GeV    RF sum with LAG - formula not quite correct dERF=" << dERF << " GeV  elemList.size()=" << elemList.size() << '\n';
  }
  else // no rf, still sum up bend loss
  {
    for(unsigned int ient=0;ient<nt.Noent();++ient) U0sum+=ngamBend[ient]*MeanSyn*EcritBend[ient];
    cout << "from summing up energy loss the " << elemList.size() << " bends U0sum=" << U0sum << " GeV" << '\n';
  }
  cout << "PowSum=" << PowerWithUnits(PowSum) << " first 250m PowSum250=" << PowerWithUnits(PowSum250) << '\n';
}

void CalcBeamSizeDivergence(Ntuple &nt,const double emitx,const double emity,const double dp,const unsigned int verbose)  // calculate beam sizes and divergences and add to ntuple
{
  const vector<string> VarList={"BETX","BETY","ALFX","ALFY","APER_1"}; if(!nt.VarsExist(VarList)) return;
  
  if(verbose) cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << __LINE__ << " emitx=" << emitx << " emity=" << emity << '\n';
  valarray<double> betx  =nt.GetVar("BETX");
  valarray<double> bety  =nt.GetVar("BETY");
  valarray<double> alfx  =nt.GetVar("ALFX");
  valarray<double> alfy  =nt.GetVar("ALFY");
  valarray<double> aper_1=nt.GetVar("APER_1");
  
  unsigned int n=nt.Noent();
  valarray<double> sigx(n);
  valarray<double> sigy(n);
  for(unsigned int i=0;i<n;++i) sigx[i]=sqrt(emitx*betx[i]);
  for(unsigned int i=0;i<n;++i) sigy[i]=sqrt(emity*bety[i]);
  // aperture -- without dispersion in sigma   and just based on aper1
  valarray<double> nsigx(n);
  valarray<double> nsigy(n);
  for(unsigned int i=0;i<n;++i) if(aper_1[i]>0 && aper_1[i]<5) nsigx[i]=aper_1[i]/sigx[i]; else nsigx[i]=99999; // check as AperOK
  for(unsigned int i=0;i<n;++i) if(aper_1[i]>0 && aper_1[i]<5) nsigy[i]=aper_1[i]/sigy[i]; else nsigy[i]=99999;
  
  valarray<double> divx(n);
  valarray<double> divy(n);
  valarray<double> gamx(n);
  valarray<double> gamy(n);
  for(unsigned int i=0;i<n;++i) { gamx[i]=(1+alfx[i]*alfx[i])/betx[i]; divx[i]=sqrt(emitx*gamx[i]); } // with courant snyder gamma
  for(unsigned int i=0;i<n;++i) { gamy[i]=(1+alfy[i]*alfy[i])/bety[i]; divy[i]=sqrt(emity*gamy[i]); } // with courant snyder gamma
  
  nt.AddVar("SIGX",sigx);
  nt.AddVar("SIGY",sigy);
  nt.AddVar("DIVX",divx);
  nt.AddVar("DIVY",divy);
  nt.AddVar("NSIGX",nsigx);
  nt.AddVar("NSIGY",nsigy);
  
  if(dp>0) // calculate also beam size and divergence with deltap and dispersion in x
    // consider to do same in y if there is significant vertical dispersion
    // more accurate calculation for Touschek see Piwinski Chao Tigner Hndbook p. 125
  {
    valarray<double> dx =nt.GetVar("DX");
    valarray<double> dpx=nt.GetVar("DPX");
    
    valarray<double> hs(n); // h(s)  http://ados.web.psi.ch/slsnotes/sls1897a.pdf  chromatic invariant or dispersion’s emittance
    for(unsigned int i=0;i<n;++i) { hs[i]=  gamx[i]*divx[i]*divx[i] + 2*alfx[i]*divx[i]*dpx[i] + betx[i]*dpx[i]*dpx[i]; }
    valarray<double> sigxd(n); // beam size with dispersion
    valarray<double> sigxp(n);
    valarray<double> divxp(n); // divergence with dispersion
    for(unsigned int i=0;i<n;++i) sigxd[i]=sqrt( emitx*betx[i] + dp*dp*dx[i]*dx[i]);
    for(unsigned int i=0;i<n;++i) { sigxp[i]=emitx* sqrt( 1+hs[i]*dp*dp/emitx )/sigx[i]; }
    for(unsigned int i=0;i<n;++i) { double t=dpx[i]+dx[i]*alfx[i]/betx[i]; divxp[i]= sqrt( emitx/betx[i] + dp*dp*t*t ); } // smaller than divx ??
    nt.AddVar("GAMX",gamx);
    nt.AddVar("GAMY",gamy);
    nt.AddVar("SIGXD",sigxd);
    nt.AddVar("SIGXP",sigxp);
    nt.AddVar("DIVXP",divxp);
    nt.AddVar("HS",hs);
  }
  else
  {
    cout << "*** dp=" << dp << " no GAMX, SIGXD, HS calculation done" << '\n';
  }
}

/*
 
 to look at a track :
 cd /tmp/$LOGNAME ; rm -f "/tmp/$LOGNAME/"*.out ; ntuple2ascii -P /tmp/hbu/trk_fname_x_1_BinFile
 
 
 
 see also  ~/c/lhc/ReadRmatrix.C
 
 MAD-X  canonical variables http://madx.web.cern.ch/madx/madX/doc/usrguide/Introduction/tables.html
 
 
 */

LinTrk::LinTrk(const string twiss_fname,const string survey_fname,unsigned int verbose) : Beam(twiss_fname,verbose) // constructor
{
  if(verbose>1) cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << __LINE__ << " verbose=" << verbose << '\n';
  this->verbose=verbose;
  the_nt=new Ntuple(ReadAndMerge(twiss_fname,survey_fname,verbose));
  if(verbose) the_nt->PrintSummary();
}

/*
 better give correct emittances in twiss
 possible to set emittances in call
 
 lintrk_b1.exN=660.e-9;  lintrk_b1.ex=lintrk_b1.Get_Part().Emit(lintrk_b1.exN); // give emittance directly
 lintrk_b1.eyN=20.e-9;   lintrk_b1.ey=lintrk_b1.Get_Part().Emit(lintrk_b1.eyN); // give emittance directly
 
 */

void LinTrk::SetSteps(unsigned int nx, double dx, unsigned int ny, double dy, unsigned int np, double dp)
{
  this->nx=nx; this->dx=dx;
  this->ny=ny; this->dy=dy;
  this->np=np; this->dp=dp;
}

void LinTrk::trk_rmatrix(const string& trk_fname)
{
  if(verbose>1) cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << __LINE__ << " verbose=" << verbose << '\n';
  valarray<double> BETX=the_nt->GetVar("BETX"); double betx0=BETX[0];
  valarray<double> BETY=the_nt->GetVar("BETY"); double bety0=BETY[0];
  if(verbose) cout << "initial BETX=" << betx0 << " BETY=" << bety0 << " m" << '\n';
  cout << "verbose=" << verbose << '\n';
  double sigx0=sqrt(betx0*ex);
  double sigy0=sqrt(bety0*ey);
  if(verbose) cout << "initial beam size sigx0=" << sigx0 << " sigy0=" << sigy0 << '\n';
  vector<Mat6x6> rmatrices=the_nt->GetRmatrices();  // get all rmatrices
  if(rmatrices.size()<1) exit(1);
  vector<string> tags;
  tags.push_back("X");
  tags.push_back("Y");
  tags.push_back("S");
  const unsigned int nvar=tags.size();
  valarray <double> data(nvar);
  unsigned int is=the_nt->ivar("S");
  
  // for(double x0=0;x0<1.
  if(verbose>1) cout << __FILE__ << " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << endl;
  Ntuple nt_trk(tags,&cout);
  if(verbose>1) cout << __FILE__ << " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << endl;
  for(unsigned int ix=0;ix<nx;++ix)
  {
    if(verbose>1) cout << __FILE__ << " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << endl;
    Vec6 vec0(ix*dx*sigx0,0, 0, 0, 0,0);
    if(verbose) cout << " i=" << 0 << " " << vec0.Print() << '\n';
    for(unsigned int i=0;i<rmatrices.size();++i) { Vec6 xi=rmatrices[i]*vec0; if(verbose && i<10) cout << " i=" << i << " " << xi.Print() << '\n'; data[0]=xi.r[0]; data[1]=xi.r[2]; data[2]=xi.r[4]+the_nt->Entry(i,is); nt_trk.fill(data); }
    if(verbose) nt_trk.PrintSummary();
    string outfname=trk_fname+"_x_"+to_string(ix)+"_BinFile";
    cout << __FILE__ << " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " outfname=" << outfname << '\n';
    nt_trk.WriteBinaryNtFile(outfname);
    nt_trk.ResetData();
  }
  for(unsigned int iy=0;iy<ny;++iy)
  {
    if(verbose>1) cout << __FILE__ << " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << endl;
    Vec6 vec0(0, 0, iy*dy*sigy0, 0, 0, 0);
    if(verbose) cout << " i=" << 0 << " " << vec0.Print() << '\n';
    for(unsigned int i=0;i<rmatrices.size();++i) { Vec6 xi=rmatrices[i]*vec0; if(verbose && i<10) cout << " i=" << i << " " << xi.Print() << '\n'; data[0]=xi.r[0]; data[1]=xi.r[2]; data[2]=xi.r[4]+the_nt->Entry(i,is); nt_trk.fill(data); }
    if(verbose) nt_trk.PrintSummary();
    string outfname=trk_fname+"_y_"+to_string(iy)+"_BinFile";
    cout << __FILE__ << " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " outfname=" << outfname << '\n';
    nt_trk.WriteBinaryNtFile(outfname);
    nt_trk.ResetData();
  }
  for(unsigned int ip=0;ip<np;++ip) // momentum offset
  {
    if(verbose>1) cout << __FILE__ << " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << endl;
    Vec6 vec0(0, 0, 0, 0, 0, ip*dp );
    if(verbose) cout << " i=" << 0 << " " << vec0.Print() << '\n';
    for(unsigned int i=0;i<rmatrices.size();++i) { Vec6 xi=rmatrices[i]*vec0; if(verbose && i<10) cout << " i=" << i << " " << xi.Print() << '\n'; data[0]=xi.r[0]; data[1]=xi.r[2]; data[2]=xi.r[4]+the_nt->Entry(i,is); nt_trk.fill(data); }
    if(verbose) nt_trk.PrintSummary();
    string outfname=trk_fname+"_p_"+to_string(ip)+"_BinFile";
    cout << __FILE__ << " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " outfname=" << outfname << '\n';
    nt_trk.WriteBinaryNtFile(outfname);
    nt_trk.ResetData();
  }
}

//------------  test with mutiple inheritence --- works, but generally not a good idea ---

LatticeAndBeam::LatticeAndBeam(Ntuple &nt,Beam &Mach) : Beam(Mach),Ntuple(nt)
{
  if(verbose>1) cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << __LINE__ << " verbose=" << verbose << '\n';
}

void LatticeAndBeam::Print_LatticeAndBeam()
{
  if(verbose>1) cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << __LINE__ << " verbose=" << verbose << '\n';
  PrintSummary(); // Ntuple PrintSummary
  Print(); // Beam Print
}
