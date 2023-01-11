//*****************************************************************************
//* =====================
//*  LCMEZZH
//* =====================
//*  
//* (Description)
//*    e+e- --> ZZH Matrix Element
//*
//* (Update Record)
//*    2007/03/29  K.Fujii	Original version in Physsim as generator.
//*    2014/12/08  J.Tian	Modified to calculate Matrix Element only.
//*
//*****************************************************************************

#include "physsim/LCMEZZH.h"

#include <sstream>
#include <iomanip>
#ifdef __PHASESPACE__
#define __NODECAY__
#endif

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "physsim/GENNumCon.h"

ClassImp(lcme::LCMEZZH)

//-----------------------------------------------------------------------------
// ==============================
//  class LCMEZZH
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
namespace lcme{
  LCMEZZH::LCMEZZH(const char *name, const char *title,
		   Double_t massHiggs,
		   Double_t polE,
		   Double_t polP)
  :LCMEBase(name,title,polE,polP),
   fZ1ModePtr (0),
   fZ2ModePtr (0),
   f1Ptr      (0),
   f2Ptr      (0),
   f3Ptr      (0),
   f4Ptr      (0),
   fZ1BosonPtr (0),
   fZ2BosonPtr (0),
   fPropagator(0)
  {
    //  Constructor of bases.  Default parameter should be initialized here
    //
    cout << "Init LCMEZZH " << endl;
    SetMass(massHiggs);
    //  SetBeamPol(polE,polP);
    Initialize();
  }
  // --------------------------
  //  D-tor
  // --------------------------
  LCMEZZH::~LCMEZZH()
  {
    delete f1Ptr;
    delete f2Ptr;
    delete f3Ptr;
    delete f4Ptr;
    delete fZ1ModePtr;
    delete fZ2ModePtr;
    delete fZ1BosonPtr;
    delete fZ2BosonPtr;
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  Get Matrix Element Squared
  // --------------------------
  Double_t LCMEZZH::GetMatrixElement2()
  {
    
    // sum final helicities combinations, average initial helicities combination
    Int_t vHelLLL[3] = {-1,-1,-1};
    Int_t vHelLLR[3] = {-1,-1,1};
    Int_t vHelLRL[3] = {-1,1,-1};
    Int_t vHelLRR[3] = {-1,1,1};
    Int_t vHelRLL[3] = {1,-1,-1};
    Int_t vHelRLR[3] = {1,-1,1};
    Int_t vHelRRL[3] = {1,1,-1};
    Int_t vHelRRR[3] = {1,1,1};
    Double_t sigmaLLL = GetMatrixElement2(vHelLLL);
    Double_t sigmaLLR = GetMatrixElement2(vHelLLR);
    Double_t sigmaLRL = GetMatrixElement2(vHelLRL);
    Double_t sigmaLRR = GetMatrixElement2(vHelLRR);
    Double_t sigmaRLL = GetMatrixElement2(vHelRLL);
    Double_t sigmaRLR = GetMatrixElement2(vHelRLR);
    Double_t sigmaRRL = GetMatrixElement2(vHelRRL);
    Double_t sigmaRRR = GetMatrixElement2(vHelRRR);
    
    Double_t weightElectron = (1.-fPolElectron)/2.;
    Double_t weightPositron = (1.+fPolPositron)/2.;
    
    Double_t sigma = 0.;
    sigma += (sigmaLLL+sigmaLLR+sigmaLRL+sigmaLRR)*weightElectron*weightPositron;
    sigma += (sigmaRLL+sigmaRLR+sigmaRRL+sigmaRRR)*(1.-weightElectron)*(1.-weightPositron);
    return (sigma);
  }
  Double_t LCMEZZH::GetMatrixElement2(Int_t vHel[])
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
  Double_t LCMEZZH::DSigmaDX()
  {
    // --------------------------------------------
    //  Phase space (&& calculate BetaBar)
    // --------------------------------------------
    // Z1 rest frame
    Double_t betaz1f = Beta(fM[1]*fM[1]/fQ2Z1,fM[2]*fM[2]/fQ2Z1);
    if (betaz1f <= 0.) return 0.;

    // Z2 rest frame
    Double_t betaz2f = Beta(fM[3]*fM[3]/fQ2Z2,fM[4]*fM[4]/fQ2Z2);
    if (betaz2f <= 0.) return 0.;
    
    // ZZ rest frame
    Double_t betaz = Beta(fQ2Z1/fQ2ZZ,fQ2Z2/fQ2ZZ);
    if (betaz <= 0.) return 0.;

    // ZZH rest frame
    Double_t betah = Beta(fQ2H/fQ2ZZH,fQ2ZZ/fQ2ZZH);
    if (betah <= 0.) return 0.;
    
    // beam
    Double_t eb     = TMath::Sqrt(fQ2ZZH)/2.;
    Double_t pb     = TMath::Sqrt((eb-kM_e)*(eb+kM_e));
    Double_t beta_e = pb/eb;
    fK[0].SetXYZT(0., 0., pb, eb);
    fK[1].SetXYZT(0., 0.,-pb, eb);
    
    // --------------------------------------------
    //  Calcuate differential cross section
    // --------------------------------------------
    Double_t s      = fQ2ZZH;
#ifdef __DEBUG__
    cerr << " fP[0] = (" 
	 << fP[0].E () << ","
	 << fP[0].Px() << ","
	 << fP[0].Py() << ","
	 << fP[0].Pz() << ")" << endl;
    cerr << " fP[1] = (" 
	 << fP[1].E () << ","
	 << fP[1].Px() << ","
	 << fP[1].Py() << ","
	 << fP[1].Pz() << ")" << endl;
    cerr << " fP[2] = (" 
	 << fP[2].E () << ","
	 << fP[2].Px() << ","
	 << fP[2].Py() << ","
	 << fP[2].Pz() << ")" << endl;
    cerr << " fP[3] = (" 
	 << fP[3].E () << ","
	 << fP[3].Px() << ","
	 << fP[3].Py() << ","
	 << fP[3].Pz() << ")" << endl;
    ANL4DVector qz = fP[2] + fP[3];
    cerr << " qz = (" 
	 << qz.E () << ","
	 << qz.Px() << ","
	 << qz.Py() << ","
	 << qz.Pz() << ")" << endl;
    ANL4DVector pcm = fP[0] + fP[1] + qz;
    cerr << " pcm = (" 
	 << pcm.E () << ","
	 << pcm.Px() << ","
	 << pcm.Py() << ","
	 << pcm.Pz() << ")" << endl;
#endif
    
    // -------------------
    //  Amplitude squared
    // -------------------
    Double_t  color = f1Ptr->GetColor() * f3Ptr->GetColor();
    Complex_t amp   = FullAmplitude();
    if (GetPropagator()) amp   *= GetHiggsPropagator(fM[0]*fM[0]);
    Double_t  amp2  = TMath::Power(abs(amp),2) * color;
    
    // -------------------
    //  Put them together
    // -------------------
    static const Int_t    kNbr  = 4;
    static const Double_t kFact = k2Pi/(TMath::Power(k4Pi,3*kNbr));
    
    Double_t identp = 1./2.;                         // identical particle factor
    Double_t dPhase = kFact * betah * betaz * betaz1f * betaz2f; // phase space factor
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
  Complex_t LCMEZZH::FullAmplitude()
  {
    Double_t gamz   = fZ1BosonPtr->GetWidth();
    
    Double_t qf1    = f1Ptr->GetCharge();
    Double_t t3f1   = f1Ptr->GetISpin();
    Double_t glfz1  = -kGz*(t3f1 - qf1*kSin2W);
    Double_t grfz1  = -kGz*(     - qf1*kSin2W);
    
    Double_t qf3    = f3Ptr->GetCharge();
    Double_t t3f3   = f3Ptr->GetISpin();
    Double_t glfz3  = -kGz*(t3f3 - qf3*kSin2W);
    Double_t grfz3  = -kGz*(     - qf3*kSin2W);
    
    static const Bool_t kIsIncoming = kTRUE;
    static const Bool_t kIsOutgoing = kFALSE;
    
    HELFermion em(fK[0], kM_e, fHelInitial[0], +1, kIsIncoming); // e-
    HELFermion ep(fK[1], kM_e, fHelInitial[1], -1, kIsOutgoing); // e+
    
    HELScalar  hs(fP[0]); // higgs
    
    HELFermion f1 (fP[1], fM[1], fHelFinal [1], +1, kIsOutgoing); // f1
    HELFermion f2b(fP[2], fM[2], fHelFinal [2], -1, kIsIncoming); // f2b
    HELVector  z1 (f2b, f1, glfz1, grfz1, kM_z, gamz);            // Z1
    
    HELFermion f3 (fP[3], fM[3], fHelFinal [3], +1, kIsOutgoing); // f3
    HELFermion f4b(fP[4], fM[4], fHelFinal [4], -1, kIsIncoming); // f4d
    HELVector  z2 (f4b, f3, glfz3, grfz3, kM_z, gamz);            // Z2
    
    Complex_t amp = AmpEEtoZZH(em, ep, z1, z2, hs);
    
    return amp;
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  AmpEEtoZZH()
  // --------------------------
  Complex_t LCMEZZH::AmpEEtoZZH(const HELFermion &em,
				const HELFermion &ep,
				const HELVector  &z1,
				const HELVector  &z2,
				const HELScalar  &hs)
  {
    Double_t  qe    = -1.;
    Double_t  game  = 0.;
    
    Double_t  t3e   = -1./2.;
    Double_t  glze  = -kGz*(t3e - qe*kSin2W);
    Double_t  grze  = -kGz*(    - qe*kSin2W);
    
    Double_t  gamz  = fZ1BosonPtr->GetWidth();
    
    Double_t  gzzh  = kGz*kM_z;
    
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
    // ZZH Production Amplitude
    //---------------------------
    HELFermion emz1  (em, z1, glze, grze, kM_e, game);
    HELVector  z2star(z2, hs, gzzh, kM_z, gamz);
    HELVertex  ampz1z2s(emz1, ep, z2star, glze, grze);
    
    HELFermion emz2(em, z2, glze, grze, kM_e, game);
    HELVector  z1star(z1, hs, gzzh, kM_z, gamz);
    HELVertex  ampz2z1s(emz2, ep, z1star, glze, grze);
    
    HELFermion emz1s (em, z1star, glze, grze, kM_e, game);
    HELVertex  ampz1sz2(emz1s, ep, z2, glze, grze);
    
    HELFermion emz2s (em, z2star, glze, grze, kM_e, game);
    HELVertex  ampz2sz1(emz2s, ep, z1, glze, grze);
    
    Complex_t amp = ampz1z2s + ampz2z1s + ampz1sz2 + ampz2sz1;
#endif /* end __PHASESPACE__ */
    
    return amp;
  }
  
  
  //_____________________________________________________________________________
  // --------------------------
  //  Initialization
  // --------------------------
  void LCMEZZH::Initialize()
  {
    // --------------------------------------------
    //  Initialize Z decay table
    // --------------------------------------------
    if (!fZ1BosonPtr) fZ1BosonPtr = new GENPDTZBoson();
    if (!fZ2BosonPtr) fZ2BosonPtr = new GENPDTZBoson();
    
    
    // Set Higgs Mass
    //  SetMass(125.);
    
    
    // Set Beam Polarisations
    //  SetBeamPol(0.,0.);
    
    //  cerr << "LCMEZZH initialized!" << endl;
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  Set Z decay mode
  // --------------------------
  void LCMEZZH::SetZDecayMode(Int_t iDecayMode1, Int_t iDecayMode2)
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
    
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  Set four-momentum of final states
  // --------------------------
  void LCMEZZH::SetMomentumFinal(TLorentzVector vLortz[])
  {
    fLortzZ1f1 = vLortz[0];
    fLortzZ1f2 = vLortz[1];
    fLortzZ2f1 = vLortz[2];
    fLortzZ2f2 = vLortz[3];
    fLortzH    = vLortz[4];
    
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
    GetVariablesInRestFrame(fLortzZ1,fLortzZ2,fQ2ZZ,fCosThetaZ,fPhiZ);    
    // ZZH rest frame
    TLorentzVector fLortzZZH = fLortzZZ + fLortzH;
    GetVariablesInRestFrame(fLortzH,fLortzZZ,fQ2ZZH,fCosTheta,fPhi);    
    // H
    fQ2H = fLortzH.M2();

    // feed to internal variables
    ANL4DVector pf1(fLortzZ1f1);
    ANL4DVector pf2(fLortzZ1f2);
    fP[1] = pf1;
    fP[2] = pf2;
    fM[1] = pf1.GetMass();
    fM[2] = pf2.GetMass();

    ANL4DVector pf3(fLortzZ2f1);
    ANL4DVector pf4(fLortzZ2f2);
    fP[3] = pf3;
    fP[4] = pf4;
    fM[3] = pf3.GetMass();
    fM[4] = pf4.GetMass();

    ANL4DVector pfh(fLortzH);
    fP[0] = pfh;
    fM[0] = pfh.GetMass();

    // beam
    Double_t eb     = TMath::Sqrt(fQ2ZZH)/2.;
    Double_t pb     = TMath::Sqrt((eb-kM_e)*(eb+kM_e));
    fK[0].SetXYZT(0., 0., pb, eb);
    fK[1].SetXYZT(0., 0.,-pb, eb);
    
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  SetHelicities
  // --------------------------
  void LCMEZZH::SetHelicities(Int_t vHel[])
  {
    Int_t iHelI  = vHel[0];
    Int_t iHelF1 = vHel[1];
    Int_t iHelF2 = vHel[2];
    
    static const Int_t kNi = 2;
    static const Int_t kIHelComb[kNi][2] = {{-1, +1},
					    {+1, -1}};
    static const Int_t kNf = 4;
    static const Int_t kFHelComb[kNf][5] = {{0, -1, +1, -1, +1},
					    {0, -1, +1, +1, -1},
					    {0, +1, -1, -1, +1},
					    {0, +1, -1, +1, -1}};

    Int_t iJCombI = 0, iJCombF = 0;
    if (iHelI == -1) {
      iJCombI = 0;
    }
    else if (iHelI == 1) {
      iJCombI = 1;
    }
    if (iHelF1 == -1 && iHelF2 == -1) {
      iJCombF = 0;
    }
    else if (iHelF1 == -1 && iHelF2 == 1) {
      iJCombF = 1;
    }
    else if (iHelF1 == 1 && iHelF2 == -1) {
      iJCombF = 2;
    }
    else if (iHelF1 == 1 && iHelF2 == 1) {
      iJCombF = 3;
    }
    
    fHelInitial[0] = kIHelComb[iJCombI][0];
    fHelInitial[1] = kIHelComb[iJCombI][1];
    fHelFinal  [0] = kFHelComb[iJCombF][0];
    fHelFinal  [1] = kFHelComb[iJCombF][1];
    fHelFinal  [2] = kFHelComb[iJCombF][2];
    fHelFinal  [3] = kFHelComb[iJCombF][3];
    fHelFinal  [4] = kFHelComb[iJCombF][4];
  }
}
