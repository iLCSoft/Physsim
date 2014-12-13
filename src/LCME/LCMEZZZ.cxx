//*****************************************************************************
//* =====================
//*  LCMEZZZ
//* =====================
//*  
//* (Description)
//*    e+e- --> ZZZ Matrix Element
//*
//* (Update Record)
//*    2007/03/29  K.Fujii	Original version in Physsim as generator.
//*    2014/12/08  J.Tian	Modified to calculate Matrix Element only.
//*
//*****************************************************************************

#include "LCMEZZZ.h"

#include <sstream>
#include <iomanip>
#ifdef __PHASESPACE__
#define __NODECAY__
#endif

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "GENNumCon.h"

ClassImp(lcme::LCMEZZZ)

//-----------------------------------------------------------------------------
// ==============================
//  class LCMEZZZ
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
namespace lcme{
  LCMEZZZ::LCMEZZZ(const char *name, const char *title,
		   Double_t polE,
		   Double_t polP)
  :LCMEBase(name,title,polE,polP),
   fZ1ModePtr (0),
   fZ2ModePtr (0),
   fZModePtr  (0),
   f1Ptr      (0),
   f2Ptr      (0),
   f3Ptr      (0),
   f4Ptr      (0),
   f5Ptr      (0),
   f6Ptr      (0),
   fZ1BosonPtr (0),
   fZ2BosonPtr (0),
   fZBosonPtr (0)
  {
    //  Constructor of bases.  Default parameter should be initialized here
    //
    cout << "Init LCMEZZZ " << endl;
    //  SetBeamPol(polE,polP);
    Initialize();
  }
  // --------------------------
  //  D-tor
  // --------------------------
  LCMEZZZ::~LCMEZZZ()
  {
    delete f1Ptr;
    delete f2Ptr;
    delete f3Ptr;
    delete f4Ptr;
    delete f5Ptr;
    delete f6Ptr;
    delete fZ1ModePtr;
    delete fZ2ModePtr;
    delete fZModePtr;
    delete fZ1BosonPtr;
    delete fZ2BosonPtr;
    delete fZBosonPtr;
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  Get Matrix Element Squared
  // --------------------------
  Double_t LCMEZZZ::GetMatrixElement2()
  {
    
    // sum final helicities combinations, average initial helicities combination
    Double_t weightElectron = (1.-fPolElectron)/2.;
    Double_t weightPositron = (1.+fPolPositron)/2.;

    Double_t sigma = 0.;
    for (Int_t i=-1;i<=1;i+=2) {
      Double_t wE = 0., wP = 0.;
      if (i==-1) {wE = weightElectron; wP = weightPositron;}        // (e-,e+)=(-1,1)
      if (i==1)  {wE = 1.-weightElectron; wP = 1.-weightPositron;}  // (e-,e+)=(1,-1)
      for (Int_t j=-1;j<=1;j+=2) {
	for (Int_t k=-1;k<=1;k+=2) {
	  for (Int_t l=-1;l<=1;l+=2) {
	    Int_t vHel[4] = {i,j,k,l};
	    sigma += wE*wP*GetMatrixElement2(vHel);
	  }
	}
      }
    }
    
    return (sigma);
  }
  Double_t LCMEZZZ::GetMatrixElement2(Int_t vHel[])
  {
    // with initial and final helicities combinations specified
    SetHelicities(vHel);
    Double_t sigma = DSigmaDX();
    return (sigma);
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  DSigmaDX
  // --------------------------
  Double_t LCMEZZZ::DSigmaDX()
  {
    // --------------------------------------------
    //  Phase space (&& calculate BetaBar)
    // --------------------------------------------
    // Z1 rest frame
    Double_t betaz1f = Beta(fM[0]*fM[0]/fQ2Z1,fM[1]*fM[1]/fQ2Z1);
    if (betaz1f <= 0.) return 0.;

    // Z2 rest frame
    Double_t betaz2f = Beta(fM[2]*fM[2]/fQ2Z2,fM[3]*fM[3]/fQ2Z2);
    if (betaz2f <= 0.) return 0.;
    
    // ZZ rest frame
    Double_t betaz1 = Beta(fQ2Z1/fQ2ZZ,fQ2Z2/fQ2ZZ);
    if (betaz1 <= 0.) return 0.;

    // Z rest frame
    Double_t betazf = Beta(fM[4]*fM[4]/fQ2Z,fM[5]*fM[5]/fQ2Z);
    if (betazf <= 0.) return 0.;

    // ZZZ rest frame
    Double_t betazz= Beta(fQ2ZZ/fQ2ZZZ,fQ2Z/fQ2ZZZ);
    if (betazz <= 0.) return 0.;
    
    // beam
    Double_t eb     = TMath::Sqrt(fQ2ZZZ)/2.;
    Double_t pb     = TMath::Sqrt((eb-kM_e)*(eb+kM_e));
    Double_t beta_e = pb/eb;
    fK[0].SetXYZT(0., 0., pb, eb);
    fK[1].SetXYZT(0., 0.,-pb, eb);
    
    // --------------------------------------------
    //  Calcuate differential cross section
    // --------------------------------------------
    Double_t s      = fQ2ZZZ;
#ifdef __DEBUG__
    for (Int_t i=0; i<6; i++) {
      cerr << " fP[" << i << "] = (" 
	   << fP[i].E () << ","
	   << fP[i].Px() << ","
	   << fP[i].Py() << ","
	   << fP[i].Pz() << ")" << endl;
    }
    ANL4DVector qz1 = fP[0] + fP[1];
    ANL4DVector qz2 = fP[2] + fP[3];
    ANL4DVector qz  = fP[4] + fP[5];
    cerr << " qz1 = (" 
	 << qz1.E () << ","
	 << qz1.Px() << ","
	 << qz1.Py() << ","
	 << qz1.Pz() << ")" << endl;
    cerr << " qz2 = (" 
	 << qz2.E () << ","
	 << qz2.Px() << ","
	 << qz2.Py() << ","
	 << qz2.Pz() << ")" << endl;
    cerr << " qz  = (" 
	 << qz.E () << ","
	 << qz.Px() << ","
	 << qz.Py() << ","
	 << qz.Pz() << ")" << endl;
    ANL4DVector pcm = qz1 + qz2 + qz;
    cerr << " mz1 = " << qz1.GetMass() << endl;
    cerr << " mz2 = " << qz2.GetMass() << endl;
    cerr << " mz  = " << qz .GetMass() << endl;
    cerr << " pcm = (" 
	 << pcm.E () << ","
	 << pcm.Px() << ","
	 << pcm.Py() << ","
	 << pcm.Pz() << ")" << endl;
#endif
    
    // -------------------
    //  Amplitude squared
    // -------------------
    Double_t  color = f1Ptr->GetColor() * f3Ptr->GetColor() * f5Ptr->GetColor();
    Complex_t amp   = FullAmplitude();
    Double_t  amp2  = TMath::Power(abs(amp),2) * color;
    
    // -------------------
    //  Put them together
    // -------------------
    static const Int_t    kNbr  = 5;
    static const Double_t kFact = k2Pi/(TMath::Power(k4Pi,3*kNbr));
    
    Double_t identp = 1./6.;                         // identical particle factor
    Double_t dPhase = kFact * betazz * betaz1 * betaz1f * betaz2f * betazf; // phase space factor
    Double_t flux   = 1./(2.* s * beta_e);           // beam flux factor
    Double_t spin   = 1.;                            // spin average for e+
    
    Double_t sigma  = identp * flux * spin * amp2 * dPhase; // in [1/GeV^2]
    sigma *= kGeV2fb;                              // now in [fb]
    
    return sigma;
  }
  
//_____________________________________________________________________________
// --------------------------
//  FullAmplitude()
// --------------------------
  Complex_t LCMEZZZ::FullAmplitude()
  {
    Double_t mz     = fZBosonPtr->GetMass();
    Double_t gamz   = fZBosonPtr->GetWidth();
    
    Double_t qf1    = f1Ptr->GetCharge();
    Double_t t3f1   = f1Ptr->GetISpin();
    Double_t glzf1  = -kGz*(t3f1 - qf1*kSin2W);
    Double_t grzf1  = -kGz*(     - qf1*kSin2W);
    
    Double_t qf3    = f3Ptr->GetCharge();
    Double_t t3f3   = f3Ptr->GetISpin();
    Double_t glzf3  = -kGz*(t3f3 - qf3*kSin2W);
    Double_t grzf3  = -kGz*(     - qf3*kSin2W);
    
    Double_t qf5    = f5Ptr->GetCharge();
    Double_t t3f5   = f5Ptr->GetISpin();
    Double_t glzf5  = -kGz*(t3f5 - qf5*kSin2W);
    Double_t grzf5  = -kGz*(     - qf5*kSin2W);
    
    static const Bool_t kIsIncoming = kTRUE;
    static const Bool_t kIsOutgoing = kFALSE;
    
    HELFermion em(fK[0], kM_e, fHelInitial[0], +1, kIsIncoming); // e-
    HELFermion ep(fK[1], kM_e, fHelInitial[1], -1, kIsOutgoing); // e+
    
    HELFermion f1 (fP[0], fM[0], fHelFinal [0], +1, kIsOutgoing); // f1
    HELFermion f2b(fP[1], fM[1], fHelFinal [1], -1, kIsIncoming); // f2b
    HELVector  z1(f2b, f1, glzf1, grzf1, mz, gamz);               // Z1
    
    HELFermion f3 (fP[2], fM[2], fHelFinal [2], +1, kIsOutgoing); // f3
    HELFermion f4b(fP[3], fM[3], fHelFinal [3], -1, kIsIncoming); // f4b
    HELVector  z2(f4b, f3, glzf3, grzf3, mz, gamz);               // Z2
    
    HELFermion f5 (fP[4], fM[4], fHelFinal [4], +1, kIsOutgoing); // f5
    HELFermion f6b(fP[5], fM[5], fHelFinal [5], -1, kIsIncoming); // f6b
    HELVector  zf(f6b, f5, glzf5, grzf5, mz, gamz);               // Z
    
    Complex_t amp = AmpEEtoZZZ(em, ep, z1, z2, zf);
    
    return amp;
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  AmpEEtoZZZ()
  // --------------------------
  Complex_t LCMEZZZ::AmpEEtoZZZ(const HELFermion &em,
				const HELFermion &ep,
				const HELVector  &z1,
				const HELVector  &z2,
				const HELVector  &zf)
  {
    Double_t  qe    = -1.;
    Double_t  t3e   = -1./2.;
    Double_t  glze  = -kGz*(t3e - qe*kSin2W);
    Double_t  grze  = -kGz*(    - qe*kSin2W);
    
    Double_t  mz    = fZBosonPtr->GetMass();
    Double_t  gamz  = fZBosonPtr->GetWidth();
    
    Double_t  me    = kMass[0][1][0];
    Double_t  game  = 0.;
    
    Double_t  gzzh  = kGz*kM_z;
    Double_t  mh    = 9999.;//125.;
    Double_t  gamh  = 4.e-3;
    
    //--------------------------------------------------
    // Calculate Amplitudes
    //--------------------------------------------------
#ifdef __PHASESPACE__
    //---------------------------
    // Just Phase Space
    //---------------------------
    Complex_t amp = 1;
#else
    //---------------------------
    // ZZZ Production Amplitude
    //---------------------------
#if 1
#if 1
    // (1)
    HELFermion emz1 (em, z1, glze, grze, me, game);
    HELFermion epz2 (ep, z2, glze, grze, me, game);
    HELVertex  amp01(emz1, epz2, zf, glze, grze);
    
    // (2)
    HELFermion epzf (ep, zf, glze, grze, me, game);
    HELVertex  amp02(emz1, epzf, z2, glze, grze);
    
    // (3)
    HELFermion emz2 (em, z2, glze, grze, me, game);
    HELFermion epz1 (ep, z1, glze, grze, me, game);
    HELVertex  amp03(emz2, epz1, zf, glze, grze);
    
    // (4)
    HELVertex  amp04(emz2, epzf, z1, glze, grze);
    
    // (5)
    HELFermion emzf (em, zf, glze, grze, me, game);
    HELVertex  amp05(emzf, epz1, z2, glze, grze);
    
    // (6)
    HELVertex  amp06(emzf, epz2, z1, glze, grze);
#else
    // (1)
    HELFermion emz1   (em  , z1, glze, grze, me, game);
    HELFermion emz1zf (emz1, zf, glze, grze, me, game);
    HELVertex  amp01(emz1zf, ep, z2, glze, grze);
    
    // (2)
    HELFermion emz1z2 (emz1, z2, glze, grze, me, game);
    HELVertex  amp02(emz1z2, ep, zf, glze, grze);
    
    // (3)
    HELFermion emz2   (em  , z2, glze, grze, me, game);
    HELFermion emz2zf (emz2, zf, glze, grze, me, game);
    HELVertex  amp03(emz2zf, ep, z1, glze, grze);
    
    // (4)
    HELFermion emz2z1 (emz2, z1, glze, grze, me, game);
    HELVertex  amp04(emz2z1, ep, zf, glze, grze);
    
    // (5)
    HELFermion emzf   (em  , zf, glze, grze, me, game);
    HELFermion emzfz2 (emzf, z2, glze, grze, me, game);
    HELVertex  amp05(emzfz2, ep, z1, glze, grze);
    
    // (6)
    HELFermion emzfz1 (emzf, z1, glze, grze, me, game);
    HELVertex  amp06(emzfz1, ep, z2, glze, grze);
#endif
#else
    HELFermion emz1 (em, z1, glze, grze, me, game);
    HELFermion emz2 (em, z2, glze, grze, me, game);
    HELFermion emzf (em, zf, glze, grze, me, game);
    
    Complex_t  amp01 = AmpEEtoZZ(emz1, ep, z2, zf);
    Complex_t  amp02 = AmpEEtoZZ(emz2, ep, z1, zf);
    Complex_t  amp03 = AmpEEtoZZ(emzf, ep, z1, z2);
    Complex_t  amp04 = 0.;
    Complex_t  amp05 = 0.;
    Complex_t  amp06 = 0.;
#endif
    
    //--
    // Higgs
    //--
    HELVector  zs (em, ep, glze, grze, mz, gamz);
    // (7)
    HELScalar  h12(z1, z2, gzzh, mh, gamh);
    HELVertex  amp07(zf, zs, h12, gzzh);
    
    // (8)
    HELScalar  h1f(z1, zf, gzzh, mh, gamh);
    HELVertex  amp08(z2, zs, h1f, gzzh);
    
    // (9)
    HELScalar  h2f(z2, zf, gzzh, mh, gamh);
    HELVertex  amp09(z1, zs, h2f, gzzh);
    
    //--
    // Sum
    //--
#if 1
    // without Higgs
    Complex_t amp  = amp01 + amp02 + amp03 + amp04 + amp05 + amp06;
#else
    Complex_t amp  = amp01 + amp02 + amp03 + amp04 + amp05 + amp06
      + amp07 + amp08 + amp09;
#endif

#endif /* end __PHASESPACE__ */

    return amp;
  }

  //_____________________________________________________________________________
  // --------------------------
  //  AmpEEtoZZ()
  // --------------------------
  Complex_t LCMEZZZ::AmpEEtoZZ(const HELFermion &em,
                               const HELFermion &ep,
                               const HELVector  &z1,
                               const HELVector  &z2)
  {
    Double_t  qe    = -1.;
    Double_t  t3e   = -1./2.;
    Double_t  glze  = -kGz*(t3e - qe*kSin2W);
    Double_t  grze  = -kGz*(    - qe*kSin2W);
    
    //    Double_t  gamz  = fZ1BosonPtr->GetWidth();
    
    Double_t  me    = kMass[0][1][0];
    Double_t  game  = 0.;
    
    //--------------------------------------------------
    // Calculate Amplitudes
    //--------------------------------------------------
#ifdef __PHASESPACE__
    //---------------------------
    // Just Phase Space
    //---------------------------
    Complex_t amp = 1;
#else
    //---------------------------
    // ZZ Production Amplitude
    //---------------------------
    //--
    // T-channel
    //--
    HELFermion emz1 (em, z1, glze, grze, me, game);
    HELVertex  amp01(emz1, ep, z2, glze, grze);
    
    //--
    // U-channel ne
    //--
    HELFermion emz2 (em, z2, glze, grze, me, game);
    HELVertex  amp02(emz2, ep, z1, glze, grze);
    
    //--
    // Sum
    //--
    Complex_t amp = amp01 + amp02;
#endif /* end __PHASESPACE__ */
    
    return amp;
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  Initialization
  // --------------------------
  void LCMEZZZ::Initialize()
  {
    // --------------------------------------------
    //  Initialize Z decay table
    // --------------------------------------------
    if (!fZ1BosonPtr) fZ1BosonPtr = new GENPDTZBoson();
    if (!fZ2BosonPtr) fZ2BosonPtr = new GENPDTZBoson();
    if (!fZBosonPtr)  fZBosonPtr  = new GENPDTZBoson();
    
    // Set Beam Polarisations
    //  SetBeamPol(0.,0.);
    
    //  cerr << "LCMEZZZ initialized!" << endl;
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  Set Z decay mode
  // --------------------------
  void LCMEZZZ::SetZDecayMode(Int_t iDecayMode1, Int_t iDecayMode2, Int_t iDecayMode)
  {
    std::cout << "set Z1 decay mode to be " << iDecayMode1 << std::endl;
    fZ1DecayMode = iDecayMode1;
    for (Int_t m=1; m<=fZ1BosonPtr->GetEntries(); m++) {
      GENDecayMode *mp = fZ1BosonPtr->GetMode(m); 
      if (mp && (m != fZ1DecayMode)) {
        mp->Lock();
      }
    }
    fZ1BosonPtr->DebugPrint();
    fZ1ModePtr = fZ1BosonPtr->GetMode(fZ1DecayMode); 
    f1Ptr = static_cast<GENPDTEntry *>(fZ1ModePtr->At(0));
    f2Ptr = static_cast<GENPDTEntry *>(fZ1ModePtr->At(1));

    std::cout << "set Z2 decay mode to be " << iDecayMode2 << std::endl;
    fZ2DecayMode = iDecayMode2;
    for (Int_t m=1; m<=fZ2BosonPtr->GetEntries(); m++) {
      GENDecayMode *mp = fZ2BosonPtr->GetMode(m); 
      if (mp && (m != fZ2DecayMode)) {
        mp->Lock();
      }
    }
    fZ2BosonPtr->DebugPrint();
    fZ2ModePtr = fZ2BosonPtr->GetMode(fZ2DecayMode); 
    f3Ptr = static_cast<GENPDTEntry *>(fZ2ModePtr->At(0));
    f4Ptr = static_cast<GENPDTEntry *>(fZ2ModePtr->At(1));
    

    std::cout << "set Z decay mode to be " << iDecayMode << std::endl;
    fZDecayMode = iDecayMode;
    for (Int_t m=1; m<=fZBosonPtr->GetEntries(); m++) {
      GENDecayMode *mp = fZBosonPtr->GetMode(m); 
      if (mp && (m != fZDecayMode)) {
        mp->Lock();
      }
    }
    fZBosonPtr->DebugPrint();
    fZModePtr = fZBosonPtr->GetMode(fZDecayMode); 
    f5Ptr = static_cast<GENPDTEntry *>(fZModePtr->At(0));
    f6Ptr = static_cast<GENPDTEntry *>(fZModePtr->At(1));
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  Set four-momentum of final states
  // --------------------------
  void LCMEZZZ::SetMomentumFinal(TLorentzVector vLortz[])
  {
    fLortzZ1f1 = vLortz[0];
    fLortzZ1f2 = vLortz[1];
    fLortzZ2f1 = vLortz[2];
    fLortzZ2f2 = vLortz[3];
    fLortzZf1  = vLortz[4];
    fLortzZf2  = vLortz[5];
    
    // --------------------------
    //  Calculate variables in rest frame based on input kinematics
    // --------------------------
    // just providing information for verification, not actually being used by calculation of ME
    // Z1 rest frame
    fLortzZ1 = fLortzZ1f1 + fLortzZ1f2;
    GetVariablesInRestFrame(fLortzZ1f1,fLortzZ1f2,fQ2Z1,fCosThetaZ1F,fPhiZ1F);
    // Z2 rest frame    
    fLortzZ2 = fLortzZ2f1 + fLortzZ2f2;
    GetVariablesInRestFrame(fLortzZ2f1,fLortzZ2f2,fQ2Z2,fCosThetaZ2F,fPhiZ2F);
    // ZZ rest frame
    TLorentzVector fLortzZZ = fLortzZ1 + fLortzZ2;
    GetVariablesInRestFrame(fLortzZ1,fLortzZ2,fQ2ZZ,fCosThetaZ1,fPhiZ1);    
    // Z rest frame
    TLorentzVector fLortzZ = fLortzZf1 + fLortzZf2;
    GetVariablesInRestFrame(fLortzZf1,fLortzZf2,fQ2Z,fCosThetaZF,fPhiZF);
    // ZZZ rest frame
    TLorentzVector fLortzZZZ = fLortzZZ + fLortzZ;
    GetVariablesInRestFrame(fLortzZZ,fLortzZ,fQ2ZZZ,fCosTheta,fPhi);    

    // feed to internal variables
    ANL4DVector pf1(fLortzZ1f1);
    ANL4DVector pf2(fLortzZ1f2);
    fP[0] = pf1;
    fP[1] = pf2;
    fM[0] = pf1.GetMass();
    fM[1] = pf2.GetMass();

    ANL4DVector pf3(fLortzZ2f1);
    ANL4DVector pf4(fLortzZ2f2);
    fP[2] = pf3;
    fP[3] = pf4;
    fM[2] = pf3.GetMass();
    fM[3] = pf4.GetMass();

    ANL4DVector pf5(fLortzZf1);
    ANL4DVector pf6(fLortzZf2);
    fP[4] = pf5;
    fP[5] = pf6;
    fM[4] = pf5.GetMass();
    fM[5] = pf6.GetMass();

    // beam
    Double_t eb     = TMath::Sqrt(fQ2ZZZ)/2.;
    Double_t pb     = TMath::Sqrt((eb-kM_e)*(eb+kM_e));
    fK[0].SetXYZT(0., 0., pb, eb);
    fK[1].SetXYZT(0., 0.,-pb, eb);
    
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  SetHelicities
  // --------------------------
  void LCMEZZZ::SetHelicities(Int_t vHel[])
  {
    Int_t iHelI  = vHel[0];
    Int_t iHelF1 = vHel[1];
    Int_t iHelF2 = vHel[2];
    Int_t iHelF3 = vHel[3];
    
    static const Int_t kNi = 2;
    static const Int_t kIHelComb[kNi][2] = {{-1, +1},
					    {+1, -1}};
    static const Int_t kNf = 8;
    static const Int_t kFHelComb[kNf][6] = {{-1, +1, -1, +1, -1, +1},
					    {-1, +1, -1, +1, +1, -1},
					    {-1, +1, +1, -1, -1, +1},
					    {-1, +1, +1, -1, +1, -1},
					    {+1, -1, -1, +1, -1, +1},
					    {+1, -1, -1, +1, +1, -1},
					    {+1, -1, +1, -1, -1, +1},
					    {+1, -1, +1, -1, +1, -1}};

    Int_t iJCombI = 0, iJCombF = 0;
    if (iHelI == -1) {
      iJCombI = 0;
    }
    else if (iHelI == 1) {
      iJCombI = 1;
    }
    if (iHelF1 == -1 && iHelF2 == -1 && iHelF3 == -1) {
      iJCombF = 0;
    }
    else if (iHelF1 == -1 && iHelF2 == -1 && iHelF3 == 1) {
      iJCombF = 1;
    }
    else if (iHelF1 == -1 && iHelF2 == 1 && iHelF3 == -1) {
      iJCombF = 2;
    }
    else if (iHelF1 == -1 && iHelF2 == 1 && iHelF3 == 1) {
      iJCombF = 3;
    }
    else if (iHelF1 == 1 && iHelF2 == -1 && iHelF3 == -1) {
      iJCombF = 4;
    }
    else if (iHelF1 == 1 && iHelF2 == -1 && iHelF3 == 1) {
      iJCombF = 5;
    }
    else if (iHelF1 == 1 && iHelF2 == 1 && iHelF3 == -1) {
      iJCombF = 6;
    }
    else if (iHelF1 == 1 && iHelF2 == 1 && iHelF3 == 1) {
      iJCombF = 7;
    }
    
    fHelInitial[0] = kIHelComb[iJCombI][0];
    fHelInitial[1] = kIHelComb[iJCombI][1];
    fHelFinal  [0] = kFHelComb[iJCombF][0];
    fHelFinal  [1] = kFHelComb[iJCombF][1];
    fHelFinal  [2] = kFHelComb[iJCombF][2];
    fHelFinal  [3] = kFHelComb[iJCombF][3];
    fHelFinal  [4] = kFHelComb[iJCombF][4];
    fHelFinal  [5] = kFHelComb[iJCombF][5];
  }
}
