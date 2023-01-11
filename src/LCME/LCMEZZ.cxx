//*****************************************************************************
//* =====================
//*  LCMEZZ
//* =====================
//*  
//* (Description)
//*    e+e- --> ZZ Matrix Element
//*
//* (Update Record)
//*    2012/03/30  K.Fujii	Original version in Physsim as generator.
//*    2014/04/08  J.Tian	Modified to calculate Matrix Element only.
//*
//*****************************************************************************

#include "physsim/LCMEZZ.h"

#include <sstream>
#include <iomanip>

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "physsim/GENNumCon.h"

ClassImp(lcme::LCMEZZ)

//-----------------------------------------------------------------------------
// ==============================
//  class LCMEZZ
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
namespace lcme{
  LCMEZZ::LCMEZZ(const char *name, const char *title,
		   Double_t polE,
		   Double_t polP,
		   Int_t    iNoDecay)
  :LCMEBase(name,title,polE,polP),
   fZ1BosonPtr (0),
   fZ1ModePtr  (0),
   f1Ptr       (0),
   f2Ptr       (0),
   fZ2BosonPtr (0),
   fZ2ModePtr  (0),
   f3Ptr       (0),
   f4Ptr       (0)
  {
    //  Constructor of bases.  Default parameter should be initialized here
    //
    cout << "Init LCMEZZ" << endl;
    //  SetBeamPol(polE,polP);
    fNoDecay = iNoDecay;
    fPropagator = 0;
    Initialize();
  }
  // --------------------------
  //  D-tor
  // --------------------------
  LCMEZZ::~LCMEZZ()
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
  Double_t LCMEZZ::GetMatrixElement2()
  {
    
    // sum final helicities combinations, average initial helicities combination
    //    Double_t weightElectron = (1.-fPolElectron)/2.;
    Double_t weightElectron = (1.-GetBeamPolE())/2.;
    //    Double_t weightPositron = (1.+fPolPositron)/2.;
    Double_t weightPositron = (1.+GetBeamPolP())/2.;

    Double_t sigma = 0.;
    for (Int_t i=-1;i<=1;i++) {
      if (i==0) continue;
      Double_t wE = 0., wP = 0.;
      if (i==-1) {wE = weightElectron; wP = weightPositron;}        // (e-,e+)=(-1,1)
      if (i==1)  {wE = 1.-weightElectron; wP = 1.-weightPositron;}  // (e-,e+)=(1,-1)
      for (Int_t j=-1;j<=1;j++) {
	if (j==0 && fNoDecay!=2) continue;
	for (Int_t k=-1;k<=1;k++) {
	  if (k==0 && fNoDecay==0) continue;
	  Int_t vHel[3] = {i,j,k};
	  sigma += wE*wP*GetMatrixElement2(vHel);
	}
      }
    }
    return (sigma);
  }
  Double_t LCMEZZ::GetMatrixElement2(Int_t vHel[])
  {
    // with initial and final helicities combinations specified
    if (!fNoDecay) {
      SetHelicities(vHel);
    }
    else if (fNoDecay == 1) {
      SetHelicitiesNoDecay(vHel);  // one Z not decay
    }
    else if (fNoDecay == 2) {
      SetHelicitiesNoDecay2(vHel); // both Z not decay
    }
    Double_t sigma = 0;
    if (GetMEType() == 1) {
      sigma = DSigmaDX();   // differential cross section
    }
    else if (GetMEType() == 2) {
      sigma = TMath::Power(abs(FullAmplitude()),2);  // squared matrix element
    }
    return (sigma);
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  DSigmaDX
  // --------------------------
  Double_t LCMEZZ::DSigmaDX()
  {
#if 0 // for debug
    for (Int_t i=0; i<3; i++) {
      cerr << " fP[" << i << "] = (" 
	   << fP[i].E () << ","
	   << fP[i].Px() << ","
	   << fP[i].Py() << ","
	   << fP[i].Pz() << ")" << endl;
    }
    ANL4DVector qz1 = fP[0] + fP[1];
    ANL4DVector qz2 = fP[2];
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
    ANL4DVector pcm = qz1 + qz2;
    cerr << " pcm = (" 
	 << pcm.E () << ","
	 << pcm.Px() << ","
	 << pcm.Py() << ","
	 << pcm.Pz() << ")" << endl;
#endif
    
    // --------------------------------------------
    //  Phase space (&& calculate BetaBar)
    // --------------------------------------------
    // Z1 rest frame
    Double_t betaz1f = 1.;
    if (fNoDecay != 2) {
      betaz1f = Beta(fM[0]*fM[0]/fQ2Z1,fM[1]*fM[1]/fQ2Z1);
      if (betaz1f <= 0.) return 0.;
    }

    // Z2 rest frame
    Double_t betaz2f = 1.;
    if (!fNoDecay) {
      betaz2f = Beta(fM[2]*fM[2]/fQ2Z2,fM[3]*fM[3]/fQ2Z2);
      if (betaz2f <= 0.) return 0.;
    }
    
    // ZZ rest frame
    Double_t betaz = Beta(fQ2Z1/fQ2ZZ,fQ2Z2/fQ2ZZ);
    if (betaz <= 0.) return 0.;

    // beam
    Double_t eb     = TMath::Sqrt(fQ2ZZ)/2.;
    Double_t pb     = TMath::Sqrt((eb-kM_e)*(eb+kM_e));
    Double_t beta_e = pb/eb;
    
    // --------------------------------------------
    //  Calcuate differential cross section
    // --------------------------------------------
    Double_t s      = fQ2ZZ;

    // -------------------
    //  Amplitude squared
    // -------------------
    Double_t  color = fNoDecay>1? 1: fNoDecay>0? f1Ptr->GetColor() : f1Ptr->GetColor()*f3Ptr->GetColor();
    Complex_t amp   = FullAmplitude();
    if (GetPropagator()) {
      if (fNoDecay > 0) {
	amp *= GetBosonPropagator(fQ2Z2,kM_z,GetZ2Width());
      }
      if (fNoDecay > 1) {
	amp *= GetBosonPropagator(fQ2Z1,kM_z,GetZ1Width());
      }
    }
    Double_t  amp2  = TMath::Power(abs(amp),2) * color;
    
    // -------------------
    //  Put them together
    // -------------------
    static const Int_t    kNbr  = fNoDecay>1? 1 : fNoDecay>0? 2:3;
    static const Double_t kFact = k2Pi/(TMath::Power(k4Pi,3*kNbr));
    
    Double_t identp = 1./2;                              // identical particle factor
    Double_t dPhase = kFact * betaz * betaz1f * betaz2f; // phase space factor
    Double_t flux   = 1./(2.* s * beta_e);               // beam flux factor
    
    Double_t sigma  = identp * flux * amp2 * dPhase; // in [1/GeV^2]
    sigma *= kGeV2fb;                                // now in [fb]

    return sigma;
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  FullAmplitude()
  // --------------------------
  Complex_t LCMEZZ::FullAmplitude()
  {
    Double_t mz     = fZ1BosonPtr->GetMass();
    Double_t gamz   = fZ1BosonPtr->GetWidth();
    
    
    static const Bool_t kIsIncoming = kTRUE;
    static const Bool_t kIsOutgoing = kFALSE;
    
    HELFermion em(fK[0], kM_e, fHelInitial[0], +1, kIsIncoming); // e-
    HELFermion ep(fK[1], kM_e, fHelInitial[1], -1, kIsOutgoing); // e+

    Complex_t amp;
    if (fNoDecay == 0) {
      Double_t qf1    = f1Ptr->GetCharge();
      Double_t t3f1   = f1Ptr->GetISpin();
      Double_t glzf1  = -kGz*(t3f1 - qf1*kSin2W);
      Double_t grzf1  = -kGz*(     - qf1*kSin2W);
      
      Double_t qf3    = f3Ptr->GetCharge();
      Double_t t3f3   = f3Ptr->GetISpin();
      Double_t glzf3  = -kGz*(t3f3 - qf3*kSin2W);
      Double_t grzf3  = -kGz*(     - qf3*kSin2W);

      HELFermion f1 (fP[0], fM[0], fHelFinal [0], +1, kIsOutgoing); // f1
      HELFermion f2b(fP[1], fM[1], fHelFinal [1], -1, kIsIncoming); // f2b
      HELVector  z1(f2b, f1, glzf1, grzf1, mz, gamz);               // Z1
      
      HELFermion f3 (fP[2], fM[2], fHelFinal [2], +1, kIsOutgoing); // fu
      HELFermion f4b(fP[3], fM[3], fHelFinal [3], -1, kIsIncoming); // fdbar
      HELVector  z2(f4b, f3, glzf3, grzf3, mz, gamz);               // Z2
      
      amp = AmpEEtoZZ(em, ep, z1, z2);
    }
    else if (fNoDecay == 1) {
      Double_t qf1    = f1Ptr->GetCharge();
      Double_t t3f1   = f1Ptr->GetISpin();
      Double_t glzf1  = -kGz*(t3f1 - qf1*kSin2W);
      Double_t grzf1  = -kGz*(     - qf1*kSin2W);
      
      HELFermion f1 (fP[0], fM[0], fHelFinal [0], +1, kIsOutgoing); // f1
      HELFermion f2b(fP[1], fM[1], fHelFinal [1], -1, kIsIncoming); // f2b
      HELVector  z1(f2b, f1, glzf1, grzf1, mz, gamz);               // Z1
      
      HELVector  z2(fP[2], fM[2], fHelFinal[2],+1);                 // Z2
      
      amp = AmpEEtoZZ(em, ep, z1, z2);
    }
    else if (fNoDecay == 2) {
      HELVector  z1(fP[0], fM[0], fHelFinal[0],+1);                 // Z1
      HELVector  z2(fP[1], fM[1], fHelFinal[1],+1);                 // Z2
      
      amp = AmpEEtoZZ(em, ep, z1, z2);
    }
    
    return amp;
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  AmpEEtoZZ()
  // --------------------------
  Complex_t LCMEZZ::AmpEEtoZZ(const HELFermion &em,
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
  void LCMEZZ::Initialize()
  {
    // --------------------------------------------
    //  Initialize Z decay table
    // --------------------------------------------
    if (! fZ1BosonPtr) fZ1BosonPtr = new GENPDTZBoson();
    if (! fZ2BosonPtr) fZ2BosonPtr = new GENPDTZBoson();
    
    // Set Beam Polarisations
    //  SetBeamPol(0.,0.);
    
    //  cerr << "LCMEZZ initialized!" << endl;
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  Set Z decay mode
  // --------------------------
  void LCMEZZ::SetZDecayMode(Int_t iDecayMode1, Int_t iDecayMode2)
  {
    if (fNoDecay == 2) {
      std::cerr << "nothing to do with Z decay" << std::endl;
      return;
    }
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
  // --------------------------
  void LCMEZZ::SetZDecayMode(Int_t iDecayMode1)
  {
    if (fNoDecay == 2) {
      std::cerr << "nothing to do with Z decay" << std::endl;
      return;
    }
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

    std::cout << "Z2 doesn't decay!"<< std::endl;
    
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  Set four-momentum of final states
  // --------------------------
  void LCMEZZ::SetMomentumFinal(TLorentzVector vLortz[])
  {
    fLortzZ1f1 = vLortz[0];
    fLortzZ1f2 = vLortz[1];
    fLortzZ2f1 = vLortz[2];
    fLortzZ2f2 = vLortz[3];
    
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
    GetVariablesInRestFrame(fLortzZ1,fLortzZ2,fQ2ZZ,fCosTheta,fPhi);    

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

    // beam
    Double_t eb     = TMath::Sqrt(fQ2ZZ)/2.;
    Double_t pb     = TMath::Sqrt((eb-kM_e)*(eb+kM_e));
    fK[0].SetXYZT(0., 0., pb, eb);
    fK[1].SetXYZT(0., 0.,-pb, eb);
    
  }
  void LCMEZZ::SetMomentumFinalNoDecay(TLorentzVector vLortz[])
  {
    fLortzZ1f1 = vLortz[0];
    fLortzZ1f2 = vLortz[1];
    fLortzZ2   = vLortz[2];
    
    // --------------------------
    //  Calculate variables in rest frame based on input kinematics
    // --------------------------
    // just providing information for verification, not actually being used by calculation of ME
    // Z1 rest frame
    fLortzZ1 = fLortzZ1f1 + fLortzZ1f2;
    GetVariablesInRestFrame(fLortzZ1f1,fLortzZ1f2,fQ2Z1,fCosThetaZ1F,fPhiZ1F);
    // Z2 rest frame
    fQ2Z2 = fLortzZ2.M2();
    // ZZ rest frame
    TLorentzVector fLortzZZ = fLortzZ1 + fLortzZ2;
    GetVariablesInRestFrame(fLortzZ1,fLortzZ2,fQ2ZZ,fCosTheta,fPhi);    

    // feed to internal variables
    ANL4DVector pf1(fLortzZ1f1);
    ANL4DVector pf2(fLortzZ1f2);
    fP[0] = pf1;
    fP[1] = pf2;
    fM[0] = pf1.GetMass();
    fM[1] = pf2.GetMass();

    ANL4DVector pz2(fLortzZ2);
    fP[2] = pz2;
    fM[2] = pz2.GetMass();

    // beam
    Double_t eb     = TMath::Sqrt(fQ2ZZ)/2.;
    Double_t pb     = TMath::Sqrt((eb-kM_e)*(eb+kM_e));
    fK[0].SetXYZT(0., 0., pb, eb);
    fK[1].SetXYZT(0., 0.,-pb, eb);
  }
  void LCMEZZ::SetMomentumFinalNoDecay2(TLorentzVector vLortz[])
  {
    fLortzZ1   = vLortz[0];
    fLortzZ2   = vLortz[1];
    
    // --------------------------
    //  Calculate variables in rest frame based on input kinematics
    // --------------------------
    // just providing information for verification, not actually being used by calculation of ME
    // Z1 rest frame
    fQ2Z1 = fLortzZ1.M2();
    // Z2 rest frame
    fQ2Z2 = fLortzZ2.M2();
    // ZZ rest frame
    TLorentzVector fLortzZZ = fLortzZ1 + fLortzZ2;
    GetVariablesInRestFrame(fLortzZ1,fLortzZ2,fQ2ZZ,fCosTheta,fPhi);    

    // feed to internal variables
    ANL4DVector pz1(fLortzZ1);
    fP[0] = pz1;
    fM[0] = pz1.GetMass();

    ANL4DVector pz2(fLortzZ2);
    fP[1] = pz2;
    fM[1] = pz2.GetMass();

    // beam
    Double_t eb     = TMath::Sqrt(fQ2ZZ)/2.;
    Double_t pb     = TMath::Sqrt((eb-kM_e)*(eb+kM_e));
    fK[0].SetXYZT(0., 0., pb, eb);
    fK[1].SetXYZT(0., 0.,-pb, eb);
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  SetHelicities
  // --------------------------
  void LCMEZZ::SetHelicities(Int_t vHel[])
  {
    Int_t iHelI  = vHel[0];
    Int_t iHelF1 = vHel[1];
    Int_t iHelF2 = vHel[2];
    
   static const Int_t kNi = 2;
   static const Int_t kIHelComb[kNi][2] = {{-1, +1},
                                           {+1, -1}};
   static const Int_t kNf = 4;
   static const Int_t kFHelComb[kNf][4] = {{-1, +1, -1, +1},
                                           {-1, +1, +1, -1},
                                           {+1, -1, -1, +1},
                                           {+1, -1, +1, -1}};
    
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
  }
  // --------------------------
  void LCMEZZ::SetHelicitiesNoDecay(Int_t vHel[])
  {
    Int_t iHelI  = vHel[0];
    Int_t iHelF1 = vHel[1];
    Int_t iHelF2 = vHel[2];
    
   static const Int_t kNi = 2;
   static const Int_t kIHelComb[kNi][2] = {{-1, +1},
                                           {+1, -1}};
   static const Int_t kNf = 6;
   static const Int_t kFHelComb[kNf][3] = {{-1, +1, -1},
                                           {-1, +1,  0},
                                           {-1, +1, +1},
                                           {+1, -1, -1},
                                           {+1, -1,  0},
                                           {+1, -1, +1}};
    
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
    else if (iHelF1 == -1 && iHelF2 == 0) {
      iJCombF = 1;
    }
    else if (iHelF1 == -1 && iHelF2 == 1) {
      iJCombF = 2;
    }
    else if (iHelF1 == 1 && iHelF2 == -1) {
      iJCombF = 3;
    }
    else if (iHelF1 == 1 && iHelF2 == 0) {
      iJCombF = 4;
    }
    else if (iHelF1 == 1 && iHelF2 == 1) {
      iJCombF = 5;
    }
    
    fHelInitial[0] = kIHelComb[iJCombI][0];
    fHelInitial[1] = kIHelComb[iJCombI][1];
    fHelFinal  [0] = kFHelComb[iJCombF][0];
    fHelFinal  [1] = kFHelComb[iJCombF][1];
    fHelFinal  [2] = kFHelComb[iJCombF][2];
  }
  // --------------------------
  void LCMEZZ::SetHelicitiesNoDecay2(Int_t vHel[])
  {
    Int_t iHelI  = vHel[0];
    Int_t iHelF1 = vHel[1];
    Int_t iHelF2 = vHel[2];
    
   static const Int_t kNi = 2;
   static const Int_t kIHelComb[kNi][2] = {{-1, +1},
                                           {+1, -1}};
   static const Int_t kNf = 9;
   static const Int_t kFHelComb[kNf][2] = {{-1,  -1},
                                           {-1,   0},
                                           {-1,  +1},
                                           { 0,  -1},
                                           { 0,   0},
                                           { 0,  +1},
                                           {+1,  -1},
                                           {+1,   0},
                                           {+1,  +1}};
    
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
    else if (iHelF1 == -1 && iHelF2 == 0) {
      iJCombF = 1;
    }
    else if (iHelF1 == -1 && iHelF2 == 1) {
      iJCombF = 2;
    }
    else if (iHelF1 == 0 && iHelF2 == -1) {
      iJCombF = 3;
    }
    else if (iHelF1 == 0 && iHelF2 == 0) {
      iJCombF = 4;
    }
    else if (iHelF1 == 0 && iHelF2 == 1) {
      iJCombF = 5;
    }
    else if (iHelF1 == 1 && iHelF2 == -1) {
      iJCombF = 6;
    }
    else if (iHelF1 == 1 && iHelF2 == 0) {
      iJCombF = 7;
    }
    else if (iHelF1 == 1 && iHelF2 == 1) {
      iJCombF = 8;
    }
    
    fHelInitial[0] = kIHelComb[iJCombI][0];
    fHelInitial[1] = kIHelComb[iJCombI][1];
    fHelFinal  [0] = kFHelComb[iJCombF][0];
    fHelFinal  [1] = kFHelComb[iJCombF][1];

  }
}
