//*****************************************************************************
//* =====================
//*  LCMEWW
//* =====================
//*  
//* (Description)
//*    e+e- --> WW Matrix Element
//*
//* (Update Record)
//*    2012/03/30  K.Fujii	Original version in Physsim as generator.
//*    2014/04/08  J.Tian	Modified to calculate Matrix Element only.
//*
//*****************************************************************************

#include "LCMEWW.h"

#include <sstream>
#include <iomanip>

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "GENNumCon.h"

ClassImp(lcme::LCMEWW)

//-----------------------------------------------------------------------------
// ==============================
//  class LCMEWW
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
namespace lcme{
  LCMEWW::LCMEWW(const char *name, const char *title,
		   Double_t polE,
		   Double_t polP,
		   Int_t    iNoDecay)
  :LCMEBase(name,title,polE,polP),
   fWmBosonPtr (0),
   fWmModePtr  (0),
   f1Ptr       (0),
   f2Ptr       (0),
   fWpBosonPtr (0),
   fWpModePtr  (0),
   f3Ptr       (0),
   f4Ptr       (0),
   fZBosonPtr  (0)
  {
    //  Constructor of bases.  Default parameter should be initialized here
    //
    cout << "Init LCMEWW" << endl;
    //  SetBeamPol(polE,polP);
    Initialize();
  }
  // --------------------------
  //  D-tor
  // --------------------------
  LCMEWW::~LCMEWW()
  {
    delete f1Ptr;
    delete f2Ptr;
    delete f3Ptr;
    delete f4Ptr;
    delete fWmModePtr;
    delete fWpModePtr;
    delete fWmBosonPtr;
    delete fWpBosonPtr;
    delete fZBosonPtr;
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  Get Matrix Element Squared
  // --------------------------
  Double_t LCMEWW::GetMatrixElement2()
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
      Int_t vHel[1] = {i};
      sigma += wE*wP*GetMatrixElement2(vHel);
    }
    return (sigma);
  }
  Double_t LCMEWW::GetMatrixElement2(Int_t vHel[])
  {
    // with initial and final helicities combinations specified
    SetHelicities(vHel);
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
  Double_t LCMEWW::DSigmaDX()
  {
#ifdef __DEBUG__
  for (Int_t i=0; i<4; i++) {
    cerr << " fP[" << i << "] = (" 
         << fP[i].E () << ","
         << fP[i].Px() << ","
         << fP[i].Py() << ","
         << fP[i].Pz() << ")" << endl;
  }
  ANL4DVector qwm = fP[0] + fP[1];
  ANL4DVector qwp = fP[2] + fP[3];
  cerr << " qwm = (" 
       << qwm.E () << ","
       << qwm.Px() << ","
       << qwm.Py() << ","
       << qwm.Pz() << ")" << endl;
  cerr << " qwp = (" 
       << qwp.E () << ","
       << qwp.Px() << ","
       << qwp.Py() << ","
       << qwp.Pz() << ")" << endl;
  ANL4DVector pcm = qwm + qwp;
  cerr << " pcm = (" 
       << pcm.E () << ","
       << pcm.Px() << ","
       << pcm.Py() << ","
       << pcm.Pz() << ")" << endl;
#endif
    
    // --------------------------------------------
    //  Phase space (&& calculate BetaBar)
    // --------------------------------------------
    // Wm rest frame
    Double_t betawmf = Beta(fM[0]*fM[0]/fQ2Wm,fM[1]*fM[1]/fQ2Wm);
    if (betawmf <= 0.) return 0.;

    // Wp rest frame
    Double_t betawpf = Beta(fM[2]*fM[2]/fQ2Wp,fM[3]*fM[3]/fQ2Wp);
    if (betawpf <= 0.) return 0.;
    
    // WW rest frame
    Double_t betaw = Beta(fQ2Wm/fQ2WW,fQ2Wp/fQ2WW);
    if (betaw <= 0.) return 0.;

    // beam
    Double_t eb     = TMath::Sqrt(fQ2WW)/2.;
    Double_t pb     = TMath::Sqrt((eb-kM_e)*(eb+kM_e));
    Double_t beta_e = pb/eb;
    
    // --------------------------------------------
    //  Calcuate differential cross section
    // --------------------------------------------
    Double_t s      = fQ2WW;

    // -------------------
    //  Amplitude squared
    // -------------------
    Double_t  color =f1Ptr->GetColor()*f3Ptr->GetColor();
    Complex_t amp   = FullAmplitude();
    Int_t     ig1   = f1Ptr->GetGenNo() - 1;
    Int_t     ig2   = f2Ptr->GetGenNo() - 1;
    Int_t     ig3   = f3Ptr->GetGenNo() - 1;
    Int_t     ig4   = f4Ptr->GetGenNo() - 1;
    Double_t  mix   = TMath::Power(kVkm[ig1][ig2]*kVkm[ig3][ig4],2);
    Double_t  amp2  = TMath::Power(abs(amp),2) * color * mix;
    
    // -------------------
    //  Put them together
    // -------------------
    static const Int_t    kNbr  = 3;
    static const Double_t kFact = k2Pi/(TMath::Power(k4Pi,3*kNbr));
    
    Double_t identp = 1.;                              // identical particle factor
    Double_t dPhase = kFact * betaw * betawmf * betawpf; // phase space factor
    Double_t flux   = 1./(2.* s * beta_e);               // beam flux factor
    
    Double_t sigma  = identp * flux * amp2 * dPhase; // in [1/GeV^2]
    sigma *= kGeV2fb;                                // now in [fb]

    return sigma;
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  FullAmplitude()
  // --------------------------
  Complex_t LCMEWW::FullAmplitude()
  {
    Double_t gamw   = fWmBosonPtr->GetWidth();
    Double_t glw    = -kGw*kSqh;
    Double_t grw    = 0.;
    
    static const Bool_t kIsIncoming = kTRUE;
    static const Bool_t kIsOutgoing = kFALSE;
    
    HELFermion em(fK[0], kM_e, fHelInitial[0], +1, kIsIncoming); // e-
    HELFermion ep(fK[1], kM_e, fHelInitial[1], -1, kIsOutgoing); // e+
    
    HELFermion f1b(fP[0], fM[0], fHelFinal [0], -1, kIsIncoming); // fubar
    HELFermion f2 (fP[1], fM[1], fHelFinal [1], +1, kIsOutgoing); // fd
    HELVector  wm(f1b, f2, glw, grw, kM_w, gamw);                 // W-
    
    HELFermion f3 (fP[2], fM[2], fHelFinal [2], +1, kIsOutgoing); // fu
    HELFermion f4b(fP[3], fM[3], fHelFinal [3], -1, kIsIncoming); // fdbar
    HELVector  wp(f4b, f3, glw, grw, kM_w, gamw);                 // W+
    
    Complex_t amp = AmpEEtoWW(em, ep, wm, wp);
    
    return amp;
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  AmpEEtoWW()
  // --------------------------
  Complex_t LCMEWW::AmpEEtoWW(const HELFermion &em,
			      const HELFermion &ep,
			      const HELVector  &wm,
			      const HELVector  &wp)
  {
    Double_t  qe    = -1.;
    Double_t  ge    = -qe*kGe;
    Double_t  glae  = ge;
    Double_t  grae  = ge;
    
    Double_t  t3e   = -1./2.;
    Double_t  glze  = -kGz*(t3e - qe*kSin2W);
    Double_t  grze  = -kGz*(    - qe*kSin2W);
    
    Double_t  gamz  = fZBosonPtr->GetWidth();
    
    Double_t  gww3  = kGw;
    Double_t  glw   = -kGw*kSqh;
    Double_t  grw   = 0.;
    
    Double_t  mne   = kMass[0][0][0];
    Double_t  gamne = 0.;
    
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
    // WW Production Amplitude
    //---------------------------
    //--
    // S-channel
    //--
    HELVector w3(em, ep, glae, grae, glze, grze, kM_z, gamz);
    HELVertex ampww3(wm, wp, w3, gww3);
    
    //--
    // T-channel ne
    //--
    HELFermion ne(em, wm, glw, grw, mne, gamne);
    HELVertex  amptee(ne, ep, wp, glw, grw);
    
    //--
    // Sum
    //--
    Complex_t amp = ampww3 + amptee;
#endif /* end __PHASESPACE__ */
    
    return amp;
  }
  

  //_____________________________________________________________________________
  // --------------------------
  //  Initialization
  // --------------------------
  void LCMEWW::Initialize()
  {
    // --------------------------------------------
    //  Initialize W decay table
    // --------------------------------------------
    if (! fWmBosonPtr) fWmBosonPtr = new GENPDTWBoson();
    if (! fWpBosonPtr) fWpBosonPtr = new GENPDTWBoson();
    if (! fZBosonPtr)  fZBosonPtr  = new GENPDTZBoson();
    
    // Set Beam Polarisations
    //  SetBeamPol(0.,0.);
    
    //  cerr << "LCMEWW initialized!" << endl;
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  Set W decay mode
  // --------------------------
  void LCMEWW::SetWDecayMode(Int_t iDecayMode1, Int_t iDecayMode2)
  {
    std::cout << "set W- decay mode to be " << iDecayMode1 << std::endl;
    fWmDecayMode = iDecayMode1;
    for (Int_t m=1; m<=fWmBosonPtr->GetEntries(); m++) {
      GENDecayMode *mp = fWmBosonPtr->GetMode(m); 
      if (mp && (m != fWmDecayMode)) {
        mp->Lock();
      }
    }
    fWmBosonPtr->DebugPrint();
    fWmModePtr = fWmBosonPtr->GetMode(fWmDecayMode); 
    f1Ptr = static_cast<GENPDTEntry *>(fWmModePtr->At(0));
    f2Ptr = static_cast<GENPDTEntry *>(fWmModePtr->At(1));

    std::cout << "set W+ decay mode to be " << iDecayMode2 << std::endl;
    fWpDecayMode = iDecayMode2;
    for (Int_t m=1; m<=fWpBosonPtr->GetEntries(); m++) {
      GENDecayMode *mp = fWpBosonPtr->GetMode(m); 
      if (mp && (m != fWpDecayMode)) {
        mp->Lock();
      }
    }
    fWpBosonPtr->DebugPrint();
    fWpModePtr = fWpBosonPtr->GetMode(fWpDecayMode); 
    f3Ptr = static_cast<GENPDTEntry *>(fWpModePtr->At(0));
    f4Ptr = static_cast<GENPDTEntry *>(fWpModePtr->At(1));
    
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  Set four-momentum of final states
  // --------------------------
  void LCMEWW::SetMomentumFinal(TLorentzVector vLortz[])
  {
    fLortzWmf1 = vLortz[0];
    fLortzWmf2 = vLortz[1];
    fLortzWpf1 = vLortz[2];
    fLortzWpf2 = vLortz[3];
    
    // --------------------------
    //  Calculate variables in rest frame based on input kinematics
    // --------------------------
    // just providing information for verification, not actually being used by calculation of ME
    // Wm rest frame
    fLortzWm = fLortzWmf1 + fLortzWmf2;
    GetVariablesInRestFrame(fLortzWmf1,fLortzWmf2,fQ2Wm,fCosThetaWmF,fPhiWmF);
    // Wp rest frame    
    fLortzWp = fLortzWpf1 + fLortzWpf2;
    GetVariablesInRestFrame(fLortzWpf1,fLortzWpf2,fQ2Wp,fCosThetaWpF,fPhiWpF);
    // WW rest frame
    TLorentzVector fLortzWW = fLortzWm + fLortzWp;
    GetVariablesInRestFrame(fLortzWm,fLortzWp,fQ2WW,fCosTheta,fPhi);    

    // feed to internal variables
    ANL4DVector pf1(fLortzWmf1);
    ANL4DVector pf2(fLortzWmf2);
    fP[0] = pf1;
    fP[1] = pf2;
    fM[0] = pf1.GetMass();
    fM[1] = pf2.GetMass();

    ANL4DVector pf3(fLortzWpf1);
    ANL4DVector pf4(fLortzWpf2);
    fP[2] = pf3;
    fP[3] = pf4;
    fM[2] = pf3.GetMass();
    fM[3] = pf4.GetMass();

    // beam
    Double_t eb     = TMath::Sqrt(fQ2WW)/2.;
    Double_t pb     = TMath::Sqrt((eb-kM_e)*(eb+kM_e));
    fK[0].SetXYZT(0., 0., pb, eb);
    fK[1].SetXYZT(0., 0.,-pb, eb);
    
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  SetHelicities
  // --------------------------
  void LCMEWW::SetHelicities(Int_t vHel[])
  {
    Int_t iHelI  = vHel[0];
    
   static const Int_t kNi = 2;
   static const Int_t kIHelComb[kNi][2] = {{-1, +1},
                                           {+1, -1}};
   static const Int_t kNf = 1;
   static const Int_t kFHelComb[kNf][4] = {{+1, -1, -1, +1}};
    
    Int_t iJCombI = 0, iJCombF = 0;
    if (iHelI == -1) {
      iJCombI = 0;
    }
    else if (iHelI == 1) {
      iJCombI = 1;
    }
    
    fHelInitial[0] = kIHelComb[iJCombI][0];
    fHelInitial[1] = kIHelComb[iJCombI][1];
    fHelFinal  [0] = kFHelComb[iJCombF][0];
    fHelFinal  [1] = kFHelComb[iJCombF][1];
    fHelFinal  [2] = kFHelComb[iJCombF][2];
    fHelFinal  [3] = kFHelComb[iJCombF][3];
  }
}
