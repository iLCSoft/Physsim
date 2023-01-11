//*****************************************************************************
//* =====================
//*  LCMEEEZ
//* =====================
//*  
//* (Description)
//*    e+e- --> EEZ Matrix Element
//*
//* (Update Record)
//*    2012/03/30  K.Fujii	Original version in Physsim as generator.
//*    2014/04/08  J.Tian	Modified to calculate Matrix Element only.
//*
//*****************************************************************************

#include "physsim/LCMEEEZ.h"

#include <sstream>
#include <iomanip>

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "physsim/GENNumCon.h"

ClassImp(lcme::LCMEEEZ)

//-----------------------------------------------------------------------------
// ==============================
//  class LCMEEEZ
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
namespace lcme{
  LCMEEEZ::LCMEEEZ(const char *name, const char *title,
		   Double_t polE,
		   Double_t polP,
		   Bool_t   iNoDecay)
  :LCMEBase(name,title,polE,polP),
   fZBosonPtr (0),
   fZModePtr  (0),
   f3Ptr      (0),
   f4Ptr      (0)
  {
    //  Constructor of bases.  Default parameter should be initialized here
    //
    cout << "Init LCMEEEZ without Z decay" << endl;
    //  SetBeamPol(polE,polP);
    fNoDecay = iNoDecay;
    Initialize();
  }
  // --------------------------
  //  D-tor
  // --------------------------
  LCMEEEZ::~LCMEEEZ()
  {
    delete f3Ptr;
    delete f4Ptr;
    delete fZModePtr;
    delete fZBosonPtr;
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  Get Matrix Element Squared
  // --------------------------
  Double_t LCMEEEZ::GetMatrixElement2()
  {
    
    // sum final helicities combinations, average initial helicities combination
    Double_t weightElectron = (1.-fPolElectron)/2.;
    Double_t weightPositron = (1.+fPolPositron)/2.;

    Double_t sigma = 0.;
    for (Int_t i=-2;i<=2&&i!=0;i++) {
      if (i==0) continue;
      Double_t wE = 0., wP = 0.;
      if (i==-2) {wE = weightElectron; wP = 1.-weightPositron;}     // (e-,e+)=(-1,-1)
      if (i==-1) {wE = weightElectron; wP = weightPositron;}        // (e-,e+)=(-1,1)
      if (i==1)  {wE = 1.-weightElectron; wP = 1.-weightPositron;}  // (e-,e+)=(1,-1)
      if (i==2)  {wE = 1.-weightElectron; wP = weightPositron;}     // (e-,e+)=(1,1)
      for (Int_t j=-2;j<=2&&j!=0;j++) {
	if (j==0) continue;
	for (Int_t k=-1;k<=1;k++) {
	  if (k==0 && !fNoDecay) continue;
	  Int_t vHel[3] = {i,j,k};
	  sigma += wE*wP*GetMatrixElement2(vHel);
	}
      }
    }
    return (sigma);
  }
  Double_t LCMEEEZ::GetMatrixElement2(Int_t vHel[])
  {
    // with initial and final helicities combinations specified
    if (!fNoDecay) {
      SetHelicities(vHel);
    }
    else {
      SetHelicitiesNoDecay(vHel);
    }
    Double_t sigma = DSigmaDX();
    return (sigma);
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  DSigmaDX
  // --------------------------
  Double_t LCMEEEZ::DSigmaDX()
  {
    // --------------------------------------------
    //  Phase space (&& calculate BetaBar)
    // --------------------------------------------
    // Z rest frame
    Double_t betaf = 1.;
    if (!fNoDecay) {
      ANL4DVector pf1(fLortzZf1);
      ANL4DVector pf2(fLortzZf2);
      fP[2] = pf1;
      fP[3] = pf2;
      fM[2] = pf1.GetMass();
      fM[3] = pf2.GetMass();
      betaf = Beta(fM[2]*fM[2]/fQ2Z,fM[3]*fM[3]/fQ2Z);
      if (betaf <= 0.) return 0.;
    }
    else {
      ANL4DVector pz(fLortzZ);
      fP[2] = pz;
      fM[2] = pz.GetMass();
    }
    
    // EE rest frame
    ANL4DVector pe1(fLortzE1);
    ANL4DVector pe2(fLortzE2);
    fP[0] = pe1;
    fP[1] = pe2;
    fM[0] = kM_e;
    fM[1] = kM_e;
    Double_t betaE = Beta(fM[0]*fM[0]/fQ2EE,fM[1]*fM[1]/fQ2EE);
    if (betaE <= 0.) return 0.;
    
    // EEZ rest frame
    Double_t betax = fNoDecay? Beta(fM[2]*fM[2]/fQ2EEZ,fQ2EE/fQ2EEZ):Beta(fQ2Z/fQ2EEZ,fQ2EE/fQ2EEZ);
    if (betax <= 0.) return 0.;
    
    // beam
    Double_t eb     = TMath::Sqrt(fQ2EEZ)/2.;
    Double_t pb     = TMath::Sqrt((eb-kM_e)*(eb+kM_e));
    Double_t beta_e = pb/eb;
    fK[0].SetXYZT(0., 0., pb, eb);
    fK[1].SetXYZT(0., 0.,-pb, eb);
    
    // --------------------------------------------
    //  Calcuate differential cross section
    // --------------------------------------------
    Double_t s      = fQ2EEZ;

    // -------------------
    //  Amplitude squared
    // -------------------
    Double_t  color = fNoDecay? 1:f3Ptr->GetColor();
    Complex_t amp   = FullAmplitude();
    Double_t  amp2  = TMath::Power(abs(amp),2) * color;
    
    // -------------------
    //  Put them together
    // -------------------
    static const Int_t    kNbr  = fNoDecay? 2:3;
    static const Double_t kFact = k2Pi/(TMath::Power(k4Pi,3*kNbr));
    
    Double_t identp = 1.;                            // identical particle factor
    Double_t dPhase = kFact * betax * betaf * betaE; // phase space factor
    Double_t flux   = 1./(2.* s * beta_e);           // beam flux factor
    
    Double_t sigma  = identp * flux * amp2 * dPhase; // in [1/GeV^2]
    sigma *= kGeV2fb;                                // now in [fb]
    
    return sigma;
  }
  
//_____________________________________________________________________________
// --------------------------
//  FullAmplitude()
// --------------------------
  Complex_t LCMEEEZ::FullAmplitude()
  {
    static const Bool_t kIsIncoming = kTRUE;
    static const Bool_t kIsOutgoing = kFALSE;
    
    HELFermion em(fK[0], kM_e,  fHelInitial[0], +1, kIsIncoming);
    HELFermion ep(fK[1], kM_e,  fHelInitial[1], -1, kIsOutgoing);
    
    HELFermion e (fP[0], fM[0], fHelFinal  [0], +1, kIsOutgoing);
    HELFermion eb(fP[1], fM[1], fHelFinal  [1], -1, kIsIncoming);

    Complex_t amp;
    if (!fNoDecay) {
      Double_t mz     = fZBosonPtr->GetMass();
      Double_t gamz   = fZBosonPtr->GetWidth();
      Double_t qf     = f3Ptr->GetCharge();
      Double_t t3f    = f3Ptr->GetISpin();
      Double_t glz    = -kGz*(t3f - qf*kSin2W);
      Double_t grz    = -kGz*(    - qf*kSin2W);
      HELFermion f (fP[2], fM[2], fHelFinal  [2], +1, kIsOutgoing);
      HELFermion fb(fP[3], fM[3], fHelFinal  [3], -1, kIsIncoming);
      HELVector  zf(fb, f, glz, grz, mz, gamz);
      amp = AmpEEtoEEZ(em, ep, e, eb, zf);
    }
    else {
      HELVector  zf(fP[2], fM[2], fHelFinal[2],+1);
      amp = AmpEEtoEEZ(em, ep, e, eb, zf);
    }
    
    return amp;
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  AmpEEtoEEZ()
  // --------------------------
  Complex_t LCMEEEZ::AmpEEtoEEZ(const HELFermion &em,
				const HELFermion &ep,
				const HELFermion &e,
				const HELFermion &eb,
				const HELVector  &zf)
  {
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
    // Z Production Amplitude
    //---------------------------
    Double_t  me    = kM_e;
    Double_t  game  = 0.;
    Double_t  mz    = fZBosonPtr->GetMass();
    Double_t  gamz  = fZBosonPtr->GetWidth();
    Double_t  qe    = -1.;
    Double_t  t3e   = -1./2.;
    Double_t  glze  = -kGz*(t3e - qe*kSin2W);
    Double_t  grze  = -kGz*(    - qe*kSin2W);
    Double_t  glae  = -kGe*qe;
    Double_t  grae  = -kGe*qe;
    
    Double_t  ebm   = TMath::Abs(em.GetFourMomentum()(0));
    Double_t  e1    = TMath::Abs(e .GetFourMomentum()(0));
    Double_t  e2    = TMath::Abs(eb.GetFourMomentum()(0));
    
    HELVector a1(ebm, e1, fSh1, fCh1, fPhi1, 
		 fHelInitial[0], fHelFinal[0],+1, kGe, me);
    HELVector z1(em, e , glze, grze, mz, gamz);
    
    HELVector a2(ebm, e2, fSh2, fCh2, fPhi2, 
		 fHelInitial[1], fHelFinal[1],-1, kGe, me);
    HELVector z2(eb, ep, glze, grze, mz, gamz);
    
    HELFermion ebv(eb, zf, glze, grze, me, game);
    Complex_t amp1 = HELVertex(ebv, ep, a1, glae, grae);
    Complex_t amp2 = HELVertex(ebv, ep, z1, glze, grze);
    
    HELFermion epv(ep, zf, glze, grze, me, game);
    Complex_t amp3 = HELVertex(eb, epv, a1, glae, grae);
    Complex_t amp4 = HELVertex(eb, epv, z1, glze, grze);
    
    HELFermion ev(e, zf, glze, grze, me, game);
    Complex_t amp5 = HELVertex(em, ev, a2, glae, grae);
    Complex_t amp6 = HELVertex(em, ev, z2, glze, grze);
    
    HELFermion emv(em, zf, glze, grze, me, game);
    Complex_t amp7 = HELVertex(emv, e, a2, glae, grae);
    Complex_t amp8 = HELVertex(emv, e, z2, glze, grze);
    
#if 1
    Complex_t amp = amp1 + amp2 + amp3 + amp4
      + amp5 + amp6 + amp7 + amp8; 
#else
    Complex_t amp = amp1 + amp3 + amp5 + amp7; 
#endif
#endif /* end __PHASESPACE__ */
    
    return amp;
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  Initialization
  // --------------------------
  void LCMEEEZ::Initialize()
  {
    // --------------------------------------------
    //  Initialize Z decay table
    // --------------------------------------------
    if (! fZBosonPtr) fZBosonPtr = new GENPDTZBoson();
    
    // Set Beam Polarisations
    //  SetBeamPol(0.,0.);
    
    //  cerr << "LCMEEEZ initialized!" << endl;
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  Set Z decay mode
  // --------------------------
  void LCMEEEZ::SetZDecayMode(Int_t iDecayMode)
  {
    if (fNoDecay) {
      std::cerr << "nothing to do with Z decay" << std::endl;
      return;
    }
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
    f3Ptr = static_cast<GENPDTEntry *>(fZModePtr->At(0));
    f4Ptr = static_cast<GENPDTEntry *>(fZModePtr->At(1));
    
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  Set four-momentum of final states
  // --------------------------
  void LCMEEEZ::SetMomentumFinal(TLorentzVector vLortz[])
  {
    fLortzE1  = vLortz[0];
    fLortzE2  = vLortz[1];
    fLortzZf1 = vLortz[2];
    fLortzZf2 = vLortz[3];
    
    // --------------------------
    //  Calculate variables in rest frame based on input kinematics
    // --------------------------
    // just providing information for verification, not actually being used by calculation of ME
    // EE rest frame
    TLorentzVector fLortzEE = fLortzE1 + fLortzE2;
    GetVariablesInRestFrame(fLortzE1,fLortzE2,fQ2EE,fCosThetaE,fPhiE);
    // Z rest frame    
    fLortzZ = fLortzZf1 + fLortzZf2;
    GetVariablesInRestFrame(fLortzZf1,fLortzZf2,fQ2Z,fCosThetaF,fPhiF);
    // EEZ rest frame
    //    TLorentzVector fLortzEEZ = fLortzEE + fLortzH;
    GetVariablesInRestFrame(fLortzEE,fLortzZ,fQ2EEZ,fCosTheta,fPhi);    

    // variables for co-linear divergence
    Double_t cosThetaE1 = fLortzE1.CosTheta();
    fSh1  = TMath::Sqrt((1.-cosThetaE1)/2);
    fCh1  = TMath::Sqrt((1.+cosThetaE1)/2);
    fPhi1 = fLortzE1.Phi();
    Double_t cosThetaE2 = fLortzE2.CosTheta();
    fSh2  = TMath::Sqrt((1.-cosThetaE2)/2);
    fCh2  = TMath::Sqrt((1.+cosThetaE2)/2);
    fPhi2 = fLortzE2.Phi();
  }
  void LCMEEEZ::SetMomentumFinalNoDecay(TLorentzVector vLortz[])
  {
    fLortzE1  = vLortz[0];
    fLortzE2  = vLortz[1];
    fLortzZ   = vLortz[2];
    
    // --------------------------
    //  Calculate variables in rest frame based on input kinematics
    // --------------------------
    // just providing information for verification, not actually being used by calculation of ME
    // EE rest frame
    TLorentzVector fLortzEE = fLortzE1 + fLortzE2;
    GetVariablesInRestFrame(fLortzE1,fLortzE2,fQ2EE,fCosThetaE,fPhiE);

    // EEZ rest frame
    //    TLorentzVector fLortzEEZ = fLortzEE + fLortzH;
    GetVariablesInRestFrame(fLortzEE,fLortzZ,fQ2EEZ,fCosTheta,fPhi);    

    // variables for co-linear divergence
    Double_t cosThetaE1 = fLortzE1.CosTheta();
    fSh1  = TMath::Sqrt((1.-cosThetaE1)/2);
    fCh1  = TMath::Sqrt((1.+cosThetaE1)/2);
    fPhi1 = fLortzE1.Phi();
    Double_t cosThetaE2 = fLortzE2.CosTheta();
    fSh2  = TMath::Sqrt((1.-cosThetaE2)/2);
    fCh2  = TMath::Sqrt((1.+cosThetaE2)/2);
    fPhi2 = fLortzE2.Phi();
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  SetHelicities
  // --------------------------
  void LCMEEEZ::SetHelicities(Int_t vHel[])
  {
    Int_t iHelI  = vHel[0];
    Int_t iHelF1 = vHel[1];
    Int_t iHelF2 = vHel[2];
    
    static const Int_t kNi = 4;
    static const Int_t kIHelComb[kNi][2] = {{-1, -1},
					    {-1, +1},
					    {+1, -1},
					    {+1, +1}};
   static const Int_t kNf = 8;
   static const Int_t kFHelComb[kNf][4] = {{-1, -1, -1, +1},
                                           {-1, -1, +1, -1},
                                           {-1, +1, -1, +1},
                                           {-1, +1, +1, -1},
                                           {+1, -1, -1, +1},
                                           {+1, -1, +1, -1},
                                           {+1, +1, -1, +1},
                                           {+1, +1, +1, -1}};
    
    Int_t iJCombI = 0, iJCombF = 0;
    if (iHelI == -2) {
      iJCombI = 0;
    }
    else if (iHelI == -1) {
      iJCombI = 1;
    }
    else if (iHelI == 1) {
      iJCombI = 2;
    }
    else if (iHelI == 2) {
      iJCombI = 3;
    }
    if (iHelF1 == -2 && iHelF2 == -1) {
      iJCombF = 0;
    }
    else if (iHelF1 == -2 && iHelF2 == 1) {
      iJCombF = 1;
    }
    else if (iHelF1 == -1 && iHelF2 == -1) {
      iJCombF = 2;
    }
    else if (iHelF1 == -1 && iHelF2 == 1) {
      iJCombF = 3;
    }
    else if (iHelF1 == 1  && iHelF2 == -1) {
      iJCombF = 4;
    }
    else if (iHelF1 == 1  && iHelF2 == 1) {
      iJCombF = 5;
    }
    else if (iHelF1 == 2  && iHelF2 == -1) {
      iJCombF = 6;
    }
    else if (iHelF1 == 2  && iHelF2 == 1) {
      iJCombF = 7;
    }
    
    fHelInitial[0] = kIHelComb[iJCombI][0];
    fHelInitial[1] = kIHelComb[iJCombI][1];
    fHelFinal  [0] = kFHelComb[iJCombF][0];
    fHelFinal  [1] = kFHelComb[iJCombF][1];
    fHelFinal  [2] = kFHelComb[iJCombF][2];
    fHelFinal  [3] = kFHelComb[iJCombF][3];
  }
  void LCMEEEZ::SetHelicitiesNoDecay(Int_t vHel[])
  {
    Int_t iHelI  = vHel[0];
    Int_t iHelF1 = vHel[1];
    Int_t iHelF2 = vHel[2];
    
    static const Int_t kNi = 4;
    static const Int_t kIHelComb[kNi][2] = {{-1, -1},
					    {-1, +1},
					    {+1, -1},
					    {+1, +1}};
   static const Int_t kNf = 12;
   static const Int_t kFHelComb[kNf][3] = {{-1, -1, -1},
                                           {-1, -1,  0},
                                           {-1, -1, +1},
                                           {-1, +1, -1},
                                           {-1, +1,  0},
                                           {-1, +1, +1},
                                           {+1, -1, -1},
                                           {+1, -1,  0},
                                           {+1, -1, +1},
                                           {+1, +1, -1},
                                           {+1, +1,  0},
                                           {+1, +1, +1}};
    
    Int_t iJCombI = 0, iJCombF = 0;
    if (iHelI == -2) {
      iJCombI = 0;
    }
    else if (iHelI == -1) {
      iJCombI = 1;
    }
    else if (iHelI == 1) {
      iJCombI = 2;
    }
    else if (iHelI == 2) {
      iJCombI = 3;
    }
    if (iHelF1 == -2 && iHelF2 == -1) {
      iJCombF = 0;
    }
    else if (iHelF1 == -2 && iHelF2 == 0) {
      iJCombF = 1;
    }
    else if (iHelF1 == -2 && iHelF2 == +1) {
      iJCombF = 2;
    }
    else if (iHelF1 == -1 && iHelF2 == -1) {
      iJCombF = 3;
    }
    else if (iHelF1 == -1 && iHelF2 == 0) {
      iJCombF = 4;
    }
    else if (iHelF1 == -1 && iHelF2 == 1) {
      iJCombF = 5;
    }
    else if (iHelF1 == 1  && iHelF2 == -1) {
      iJCombF = 6;
    }
    else if (iHelF1 == 1  && iHelF2 == 0) {
      iJCombF = 7;
    }
    else if (iHelF1 == 1  && iHelF2 == +1) {
      iJCombF = 8;
    }
    else if (iHelF1 == 2  && iHelF2 == -1) {
      iJCombF = 9;
    }
    else if (iHelF1 == 2  && iHelF2 == 0) {
      iJCombF = 10;
    }
    else if (iHelF1 == 2  && iHelF2 == +1) {
      iJCombF = 11;
    }
    
    fHelInitial[0] = kIHelComb[iJCombI][0];
    fHelInitial[1] = kIHelComb[iJCombI][1];
    fHelFinal  [0] = kFHelComb[iJCombF][0];
    fHelFinal  [1] = kFHelComb[iJCombF][1];
    fHelFinal  [2] = kFHelComb[iJCombF][2];
  }
}
