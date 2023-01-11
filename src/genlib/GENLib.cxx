//*****************************************************************************
//* =====================
//*  GENLib 
//* =====================
//*  
//* (Description)
//*    Class library for event generators
//*
//* (Update Record)
//*    2007/01/27  K.Fujii	Original version.
//*
//*****************************************************************************

#include "physsim/GENLib.h"
#include "physsim/GENNumCon.h"

#include <sstream>
#include <iomanip>
//#define __DEBUG__


ClassImp(GENBranch)
ClassImp(GENFrame)
ClassImp(GENPhase2)
ClassImp(GENDecayMode)
ClassImp(GENModePicker)
ClassImp(GENPDTEntry)
ClassImp(GENPDTWBoson)
ClassImp(GENPDTZBoson)
ClassImp(GENPDTPhoton)
ClassImp(GENPDTGluon)

//-----------------------------------------------------------------------------
// ==============================
//  class GENDecayMode
// ==============================
//_____________________________________________________________________________
//_____________________________________________________________________________
// --------------------------
//  DebugPrint
// --------------------------
void GENDecayMode::DebugPrint(const Option_t *)
{
   using namespace std;
   cerr << " Gamma = " << setw(11) << setprecision(9) 
                       << setiosflags(ios::fixed)
                       << setiosflags(ios::showpoint)
        << GetGamma() 
        << " BR = "    << setw(11) << setprecision(9)
                       << setiosflags(ios::fixed)
                       << setiosflags(ios::showpoint)
        << GetBR()     
        << " : " << (IsLocked() ? "off" : "on ")
        << " : --> ";
   TIter next(this);
   GENPDTEntry *ep;
   while ((ep = static_cast<GENPDTEntry *>(next()))) {
      cerr <<  ep->GetName() << (ep->GetPID() > 0 ? " " : "bar ");
   }
   cerr << endl;
}

//-----------------------------------------------------------------------------
// ==============================
//  class GENModePicker
// ==============================
//_____________________________________________________________________________
// --------------------------
//  Add
// --------------------------
void GENModePicker::Add(GENDecayMode *mp)
{
   TObjArray::Add(mp);
   fDone = kFALSE;
}

//_____________________________________________________________________________
// --------------------------
//  PickMode
// --------------------------
GENDecayMode *GENModePicker::PickMode(Double_t x,
                                      Double_t &weight,
                                      Int_t    &mode)
{
   if (!fDone) Update();
   mode = 0;
   TIter next(this);
   GENDecayMode *mp;
   while ((mp = static_cast<GENDecayMode *>(next()))) {
      if ((!mp->IsLocked()) && x*fBRsum <= mp->fCumBR) {
         weight = fBRsum/mp->fBR; 
         break;
      }
      mode++;
   }
   return static_cast<GENDecayMode *>(mp);
}

//_____________________________________________________________________________
// --------------------------
//  GetMode
// --------------------------
const GENDecayMode *GENModePicker::GetMode(Int_t m) const
{
   if (m<1 || m>GetEntries()) return 0; 
   return static_cast<const GENDecayMode *>(At(m-1));
}

GENDecayMode *GENModePicker::GetMode(Int_t m)
{
   fDone = kFALSE;
   if (m<1 || m>GetEntries()) return 0; 
   return static_cast<GENDecayMode *>(At(m-1));
}

//_____________________________________________________________________________
// --------------------------
//  Update
// --------------------------
void GENModePicker::Update()
{
   if (fDone) return;
   fDone = kTRUE;

   fGamma = 0.;
   TIter next(this);
   GENDecayMode *mp;
   while ((mp = static_cast<GENDecayMode *>(next()))) {
      fGamma += mp->fGamma;
   }

   Double_t cum = 0.;
   next.Reset();
   while ((mp = static_cast<GENDecayMode *>(next()))) {
      mp->fBR = mp->fGamma/fGamma;
      if (!mp->IsLocked()) cum += mp->fBR;
      mp->fCumBR = cum;
   }
   fBRsum = cum;
}

//-----------------------------------------------------------------------------
// ==============================
//  class GENPDTEntry
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
GENPDTEntry::GENPDTEntry(const Char_t     *name,
                               Int_t       pid,
                               Double_t    charge,
                               Double_t    spin,
                               Double_t    mass,
                               Int_t       gen,
                               Double_t    ispin,
                               Double_t    color)
           : fName(name),
             fPID(pid),
             fCharge(charge),
             fSpin(spin),
             fMass(mass),
             fGen(gen),
             fIsoSpin(ispin),
             fColor(color)
{
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
GENPDTEntry::~GENPDTEntry()
{
   TIter next(this);
   GENDecayMode *dmp;
   while ((dmp = static_cast<GENDecayMode *>(next()))) {
      dmp->Delete();
   }
   Delete();
}

//_____________________________________________________________________________
// --------------------------
//  GetQ2BW
// --------------------------
Double_t GENPDTEntry::GetQ2BW(Double_t    qmin,   // Q_min
                              Double_t    qmax,   // Q_max
                              Double_t       x,   // integration variable
                              Double_t &weight)   // Jacobian weight
{
   Double_t m    = GetMass ();
   Double_t gm   = GetWidth();
   Double_t mgm  = m*gm;
   Double_t m2   = m*fMass;
   Double_t mgm2 = mgm*mgm;

   Double_t thmin = TMath::ATan((qmin-m)*(qmin+m)/mgm);
   Double_t thmax = TMath::ATan((qmax-m)*(qmax+m)/mgm);
   Double_t theta = thmin + (thmax - thmin)*x;
   Double_t q2    = mgm * TMath::Tan(theta) + m2;

   weight = (thmax - thmin)*(TMath::Power(q2-m2,2) + mgm2)/mgm;

   return q2;
}

//_____________________________________________________________________________
// --------------------------
//  SetQ2BW
// --------------------------
void GENPDTEntry::SetQ2BW(Double_t    qmin,   // Q_min
                              Double_t    qmax,   // Q_max
                              Double_t      q2,   // Q2
                              Double_t &weight)   // Jacobian weight
{
   Double_t m    = GetMass ();
   Double_t gm   = GetWidth();
   Double_t mgm  = m*gm;
   Double_t m2   = m*fMass;
   Double_t mgm2 = mgm*mgm;

   Double_t thmin = TMath::ATan((qmin-m)*(qmin+m)/mgm);
   Double_t thmax = TMath::ATan((qmax-m)*(qmax+m)/mgm);
   //   Double_t theta = thmin + (thmax - thmin)*x;
   //   Double_t q2    = mgm * TMath::Tan(theta) + m2;

   weight = (thmax - thmin)*(TMath::Power(q2-m2,2) + mgm2)/mgm;

   //   return q2;
}

//_____________________________________________________________________________
// --------------------------
//  DebugPrint
// --------------------------
void GENPDTEntry::DebugPrint(const Option_t *opt)
{
   using namespace std;
   cerr << "=======================" << endl
        << " " << fName              << endl
        << "=======================" << endl;
   cerr << " --------------------------------------------------------- " << endl;
   Update();
   TIter next(this);
   GENDecayMode *mp;
   while ((mp = static_cast<GENDecayMode *>(next()))) {
      mp->DebugPrint(opt);
   }
   cerr << " --------------------------------------------------------- " << endl
        << " Gamma_tot = " << GetWidth() << " [GeV]"                     << endl
        << " Mass      = " << GetMass()  << " [GeV]"                     << endl
	<< endl;
}

//-----------------------------------------------------------------------------
// ==============================
//  class GENPDTWBoson
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
GENPDTWBoson::GENPDTWBoson()
{
   Initialize();
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
GENPDTWBoson::~GENPDTWBoson()
{
}

//_____________________________________________________________________________
// --------------------------
//  Initialize
// --------------------------
void GENPDTWBoson::Initialize()
{
   fName   = TString("W");
   fPID    = 24;
   fCharge = +1.;
   fSpin   =  1.;
   fMass   = kM_w;
   //--
   //  Set decay modes
   //--
   for (Int_t ig=0; ig<3; ig++) {
      const Char_t   *namu = kName[0][0][ig];
      const Char_t   *namd = kName[0][1][ig];
            Int_t     pidu = kPID [0][0][ig];  
            Int_t     pidd = kPID [0][1][ig];  
            Double_t  qfu  = 0.;
            Double_t  qfd  = -1.;
            Double_t  spin = 0.5;
            Double_t  t3u  = +0.5;
            Double_t  t3d  = -0.5;
            Double_t  mu   = kMass[0][0][ig];
            Double_t  md   = kMass[0][1][ig];
            Double_t  cf   = 1.;
            Double_t  vff  = 1.;
      Double_t gam  = GamToFF(mu, md, vff, cf);
      if (gam == 0.) continue;
      GENPDTEntry  *d1p = new GENPDTEntry(namu, pidu, qfu, spin, mu, ig+1,  t3u, cf);
      GENPDTEntry  *d2p = new GENPDTEntry(namd,-pidd,-qfd, spin, md, ig+1, -t3d, cf);
      GENDecayMode *dmp = new GENDecayMode(gam);
      dmp->Add(d1p);
      dmp->Add(d2p);
      Add(dmp); 
   }
   for (Int_t igu=0; igu<3; igu++) {
      for (Int_t igd=0; igd<3; igd++) {
         const Char_t   *namu = kName[1][0][igu];
         const Char_t   *namd = kName[1][1][igd];
               Int_t     pidu = kPID [1][0][igu];  
               Int_t     pidd = kPID [1][1][igd];  
               Double_t  t3u  = +0.5;
               Double_t  t3d  = -0.5;
               Double_t  qfu  = kChrg[1][0][igu];
               Double_t  qfd  = kChrg[1][1][igd];
               Double_t  spin = 0.5;
               Double_t  mu   = kMass[1][0][igu];
               Double_t  md   = kMass[1][1][igd];
               Double_t  cf   = 3*(1. + kAlphaS/kPi);
               Double_t  vff  = kVkm [igu][igd];
         Double_t gam  = GamToFF(mu, md, vff, cf);
         if (gam == 0.) continue;
         GENPDTEntry  *d1p  = new GENPDTEntry(namu, pidu, qfu, spin, mu, igu+1,  t3u, cf);
         GENPDTEntry  *d2p  = new GENPDTEntry(namd,-pidd,-qfd, spin, md, igd+1, -t3d, cf);
         GENDecayMode *dmp = new GENDecayMode(gam);
         dmp->Add(d1p);
         dmp->Add(d2p);
         Add(dmp); 
      }
   }
}

//_____________________________________________________________________________
// --------------------------
//  GamToFF
// --------------------------
Double_t GENPDTWBoson::GamToFF(Double_t mu,  // up type mass
                               Double_t md,  // down type mass
                               Double_t vff, // Vkm
                               Double_t cf)  // color factor
{
   if (mu+md >= fMass) return 0.;
   Double_t mw2  = fMass*fMass;
   Double_t mu2  = mu*mu;
   Double_t md2  = md*md;
   Double_t p1p2 = (mw2 - mu2 - md2)/2.;
   Double_t x1   = mu2/mw2;
   Double_t x2   = md2/mw2;
   Double_t p1   = (fMass/2.)*TMath::Sqrt(1. - 2*(x1+x2) + (x1-x2)*(x1-x2));

   Double_t tta  = (p1p2 + 2*(mu2+p1p2)*(md2+p1p2)/mw2)/3.;
   Double_t fac  = (kAlpha/kSin2W)/2;
   Double_t gam  = fac*tta*cf*vff*vff*p1/mw2;

   return gam;
}

//-----------------------------------------------------------------------------
// ==============================
//  class GENPDTZBoson
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
GENPDTZBoson::GENPDTZBoson()
{
   Initialize();
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
GENPDTZBoson::~GENPDTZBoson()
{
}

//_____________________________________________________________________________
// --------------------------
//  Initialize
// --------------------------
void GENPDTZBoson::Initialize()
{
   fName   = TString("Z");
   fPID    = 23;
   fCharge = 0.;
   fSpin   = 1.;
   fMass   = kM_z;
   //--
   //  Set decay modes
   //--
   for (Int_t ic=0; ic<2; ic++) {
      for (Int_t it=0; it<2; it++) { 
         for (Int_t ig=0; ig<3; ig++) {
            const Char_t   *name = kName[ic][it][ig];
                  Int_t     pid  = kPID [ic][it][ig];  
                  Double_t  t3   = (1 - 2*it)/2.;
                  Double_t  qf   = kChrg[ic][it][ig];
                  Double_t  spin = 0.5;
                  Double_t  cf   = 2*ic + 1;
                  Double_t  mass = kMass[ic][it][ig];
            GENDecayMode *dmp;
            GENPDTEntry  *d1p, *d2p;
            if (ic) cf  *= 1 + kAlphaS/kPi;
            d1p  = new GENPDTEntry(name, pid, qf, spin, mass, ig+1,  t3, cf);
            d2p  = new GENPDTEntry(name,-pid,-qf, spin, mass, ig+1, -t3, cf);
            Double_t gam  = GamToFF(t3, qf, cf, mass);
            if (gam == 0.) continue;
            dmp = new GENDecayMode(gam);
            dmp->Add(d1p);
            dmp->Add(d2p);
            Add(dmp); 
         }
      }
   }
}

//_____________________________________________________________________________
// --------------------------
//  GamToFF
// --------------------------
Double_t GENPDTZBoson::GamToFF(Double_t t3, // weak isospin
                               Double_t qf, // charge
                               Double_t cf, // color factor
                               Double_t m)  // mass
{
   Double_t mz2 = fMass*fMass;
   Double_t p1  = mz2/4 - m*m;

   if (p1 <= 0.) return 0.;

   p1 = TMath::Sqrt(p1);
   Double_t gv  =  t3/2 - kSin2W*qf;
   Double_t ga  = -t3/2;
   Double_t tta = 2*((gv*gv + ga*ga)*(mz2 - 4*p1*p1/3) - 4*ga*ga*m*m); 
   Double_t fac = (kAlpha/kSin2W/kCos2W)/2;
   Double_t gam = fac*tta*cf*p1/mz2;

   return gam;
}

//-----------------------------------------------------------------------------
// ==============================
//  class GENPDTPhoton
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
GENPDTPhoton::GENPDTPhoton()
{
   Initialize();
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
GENPDTPhoton::~GENPDTPhoton()
{
}

//_____________________________________________________________________________
// --------------------------
//  Initialize
// --------------------------
void GENPDTPhoton::Initialize()
{
   fName   = TString("Gamma");
   fPID    = 22;
   fCharge = 0.;
   fSpin   = 1.;
   fMass   = 0.;
}

//-----------------------------------------------------------------------------
// ==============================
//  class GENPDTGluon
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
GENPDTGluon::GENPDTGluon()
{
   Initialize();
}

//_____________________________________________________________________________
// --------------------------
//  D-tor
// --------------------------
GENPDTGluon::~GENPDTGluon()
{
}

//_____________________________________________________________________________
// --------------------------
//  Initialize
// --------------------------
void GENPDTGluon::Initialize()
{
   fName   = TString("Gluon");
   fPID    = 21;
   fCharge = 0.;
   fSpin   = 1.;
   fMass   = 0.;
   fColor  = 8.;
}

//-----------------------------------------------------------------------------
// ==============================
//  class GENBranch
// ==============================

//-----------------------------------------------------------------------------
// ==============================
//  class GENBranch
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
GENBranch::GENBranch(Double_t q2,
                     Double_t costh,
                     Double_t phi,
                     Double_t m12,
                     Double_t m22)
         : fQ2(q2),
           fCosTheta(costh),
           fPhi(phi), 
           fM12(m12),
           fM22(m22),
           fBR1Ptr(0),
           fBR2Ptr(0)
{
   Double_t x1 = fM12/fQ2;
   Double_t x2 = fM22/fQ2;
   fBetaBar    = TMath::Sqrt(1. - 2*(x1+x2) + TMath::Power(x1-x2,2));
}

GENBranch::GENBranch(Double_t   q2,
                     Double_t   costh,
                     Double_t   phi,
                     GENBranch *br1p,
                     GENBranch *br2p)
         : fQ2(q2),
           fCosTheta(costh),
           fPhi(phi), 
           fM12(br1p->GetQ2()),
           fM22(br2p->GetQ2()),
           fBR1Ptr(br1p),
           fBR2Ptr(br2p)
{
   Double_t x1 = fM12/fQ2;
   Double_t x2 = fM22/fQ2;
   fBetaBar    = TMath::Sqrt(1. - 2*(x1+x2) + TMath::Power(x1-x2,2));
}

GENBranch::GENBranch(Double_t   q2,
                     Double_t   costh,
                     Double_t   phi,
                     GENBranch *br1p,
                     Double_t   m22)
         : fQ2(q2),
           fCosTheta(costh),
           fPhi(phi), 
           fM12(br1p->GetQ2()),
           fM22(m22),
           fBR1Ptr(br1p),
           fBR2Ptr(0)
{
   Double_t x1 = fM12/fQ2;
   Double_t x2 = fM22/fQ2;
   fBetaBar    = TMath::Sqrt(1. - 2*(x1+x2) + TMath::Power(x1-x2,2));
}

GENBranch::GENBranch(Double_t   q2,
                     Double_t   costh,
                     Double_t   phi,
                     Double_t   m12,
                     GENBranch *br2p)
         : fQ2(q2),
           fCosTheta(costh),
           fPhi(phi), 
           fM12(m12),
           fM22(br2p->GetQ2()),
           fBR1Ptr(0),
           fBR2Ptr(br2p)
{
   Double_t x1 = fM12/fQ2;
   Double_t x2 = fM22/fQ2;
   fBetaBar    = TMath::Sqrt(1. - 2*(x1+x2) + TMath::Power(x1-x2,2));
}

// ==============================
//  class GENFrame
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
GENFrame::GENFrame()
{
   fEV[0].SetXYZ(1., 0., 0.);
   fEV[1].SetXYZ(0., 1., 0.);
   fEV[2].SetXYZ(0., 0., 1.);
}

GENFrame::GENFrame(const ANL4DVector &q, const GENFrame &eb)
{
   fEV[2] = q.Vect().Unit();
   fEV[1] = eb.fEV[2].Cross(fEV[2]);
   Double_t ae2  = fEV[1].Mag();
   static const Double_t kXmin = 1.e-12;
   if (ae2 < kXmin) {
      fEV[0] = eb.fEV[0]; fEV[1] = eb.fEV[1]; fEV[2] = eb.fEV[2];
      Double_t csth = fEV[2] * eb.fEV[2];
      if (csth <= 0.) {
         fEV[2] = -eb.fEV[2];
      }
      return;
   } else {
      fEV[1] = fEV[1].Unit();
   }
   fEV[0] = fEV[1].Cross(fEV[2]).Unit();
}

//_____________________________________________________________________________
// --------------------------
//  Transform
// --------------------------
ANL4DVector GENFrame::Transform(const ANL4DVector &pb)
{
   TVector3 pb3v = pb.X()*fEV[0] + pb.Y()*fEV[1] + pb.Z()*fEV[2];
   return ANL4DVector(pb.E(),pb3v.Px(), pb3v.Py(), pb3v.Pz());
}

//-----------------------------------------------------------------------------
// ==============================
//  class GENPhase2
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
GENPhase2::GENPhase2(const ANL4DVector &q,
                           Double_t     m12,
                           Double_t     m22,
                     const GENFrame    &eb,
                           Double_t     costh,
                           Double_t     phi,
                           Int_t        mode)
         : fQ(q), fM12(m12), fM22(m22), 
           fEb(eb), fEa(eb), 
           fCosTheta(costh), fPhi(phi),
           fBetaBar(0.), 
           fMode(mode),
           fDone(kFALSE)
{
}

//_____________________________________________________________________________
// --------------------------
//  GetFourMomentum
// --------------------------
ANL4DVector GENPhase2::GetFourMomentum(Int_t i)
{
   if (!fDone) Update(); return i ? fP2 : fP1;
}

//_____________________________________________________________________________
// --------------------------
//  GetFrame
// --------------------------
GENFrame GENPhase2::GetFrame(Int_t i)
{
   if (!fDone) Update(); return i ? fEa : fEb;
}

//_____________________________________________________________________________
// --------------------------
//  GetBetaBar
// --------------------------
Double_t GENPhase2::GetBetaBar()
{
   if (!fDone) Update(); return fBetaBar;
}

//_____________________________________________________________________________
// --------------------------
//  Beta2
// --------------------------
Double_t GENPhase2::Beta2(Double_t x1, Double_t x2)
{
   return 1. - 2*(x1+x2) + (x1-x2)*(x1-x2);
}

//_____________________________________________________________________________
// --------------------------
//  Update
// --------------------------
void GENPhase2::Update()
{
   // -------------------------
   //  Check if update needed
   // -------------------------
   if (fDone) return;
   fDone = kTRUE;

   // -------------------------
   //  Calculate Beta_bar
   // -------------------------
   Double_t amq2 = fQ.Mag2();
   if (amq2 <= 0.) {
      cerr << " >>>>> Error in GENPhase2::GENPhase2() >>>>>>>> " << endl
           << " q = ("  << fQ.E() << ", " 
                        << fQ.X() << ", " 
                        << fQ.Y() << ", "
                        << fQ.Z() << ")" << endl
           << " q2  = " << amq2   
           << " m12 = " << fM12
           << " m22 = " << fM22 << endl;
      fBetaBar = 0.;
      return;
   }

   Double_t amq = TMath::Sqrt(amq2);
   fBetaBar = Beta2(fM12/amq2, fM22/amq2);   
   if (fBetaBar < 0.) {
      fBetaBar = 0.;
      return;
   }
   fBetaBar = TMath::Sqrt(fBetaBar);

   // -------------------------
   //  Daughter momenta
   // -------------------------
   Double_t ap1  = (amq/2) * fBetaBar;
   Double_t snth = TMath::Sqrt((1.-fCosTheta)*(1.+fCosTheta));
   fP1.SetXYZT(ap1*snth*TMath::Cos(fPhi),
               ap1*snth*TMath::Sin(fPhi),
               ap1*fCosTheta,
               TMath::Sqrt(ap1*ap1+fM12));
   fP2.SetXYZT(-fP1.X(), -fP1.Y(), -fP1.Z(), amq - fP1.E());

   // -------------------------
   //  Boost them to lab. frame
   // -------------------------
   if (!fMode) {
      fEa = fEb;
   } else {
      fEa = GENFrame(fQ, fEb);
      fP1 = fEa.Transform(fP1);
      fP2 = fEa.Transform(fP2);
      TVector3 boostv = fQ.BoostVector();
      fP1.Boost(boostv);
      fP2.Boost(boostv);
   }

   // -------------------------
   //  Fix round-off errors
   // -------------------------
   if (fP1.E() <= 0. || fP2.E() <= 0.) {
      fBetaBar = 0;
   } else {
      Double_t ap = fP1.Vect().Mag();
      if (fP1.E() < ap) fP1.SetE(TMath::Sqrt(ap*ap+fM12));
               ap = fP2.Vect().Mag();
      if (fP2.E() < ap) fP2.SetE(TMath::Sqrt(ap*ap+fM22));
   }
}
