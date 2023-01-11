//*****************************************************************************
//* =====================
//*  LCMEZHH
//* =====================
//*  
//* (Description)
//*    e+e- --> ZHH Matrix Element
//*
//* (Update Record)
//*    2007/03/29  K.Fujii	Original version in Physsim as generator.
//*    2014/01/14  J.Tian	Modified to calculate Matrix Element only.
//*
//*****************************************************************************

#include "physsim/LCMEZHH.h"

#include <sstream>
#include <iomanip>
#ifdef __PHASESPACE__
#define __NODECAY__
#endif

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "physsim/GENNumCon.h"

ClassImp(lcme::LCMEZHH)

//-----------------------------------------------------------------------------
// ==============================
//  class LCMEZHH
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
namespace lcme{
  LCMEZHH::LCMEZHH(const char *name, const char *title,
		   Double_t massHiggs,
		   Double_t polE,
		   Double_t polP)
  :LCMEBase(name,title,polE,polP),
   fZModePtr  (0),
   f3Ptr      (0),
   f4Ptr      (0),
   fZBosonPtr (0),
   fPropagator(0)
  {
    //  Constructor of bases.  Default parameter should be initialized here
    //
    cout << "Init LCMEZHH " << endl;
    SetMass(massHiggs);
    //  SetBeamPol(polE,polP);
    Initialize();
  }
  // --------------------------
  //  D-tor
  // --------------------------
  LCMEZHH::~LCMEZHH()
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
  Double_t LCMEZHH::GetMatrixElement2()
  {
    
    // sum final helicities combinations, average initial helicities combination
    Int_t vHelLL[2] = {-1,-1};
    Int_t vHelLR[2] = {-1,1};
    Int_t vHelRL[2] = {1,-1};
    Int_t vHelRR[2] = {1,1};
    Double_t sigmaLL = GetMatrixElement2(vHelLL);
    Double_t sigmaLR = GetMatrixElement2(vHelLR);
    Double_t sigmaRL = GetMatrixElement2(vHelRL);
    Double_t sigmaRR = GetMatrixElement2(vHelRR);
    
    Double_t weightElectron = (1.-fPolElectron)/2.;
    Double_t weightPositron = (1.+fPolPositron)/2.;
    
    Double_t sigma = 0.;
    sigma += (sigmaLL+sigmaLR)*weightElectron*weightPositron;
    sigma += (sigmaRL+sigmaRR)*(1.-weightElectron)*(1.-weightPositron);
    return (sigma);
  }
  Double_t LCMEZHH::GetMatrixElement2(Int_t vHel[])
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
  Double_t LCMEZHH::DSigmaDX()
  {
    // --------------------------------------------
    //  Phase space (&& calculate BetaBar)
    // --------------------------------------------
    // HH rest frame
    ANL4DVector ph1(fLortzH1);
    ANL4DVector ph2(fLortzH2);
    fP[0] = ph1;
    fM[0] = ph1.GetMass();
    fP[1] = ph2;
    fM[1] = ph2.GetMass();
    Double_t betah = Beta(fM[0]*fM[0]/fQ2HH,fM[1]*fM[1]/fQ2HH);
    if (betah <= 0.) return 0.;
    
    // Z rest frame
    ANL4DVector pf1(fLortzZf1);
    ANL4DVector pf2(fLortzZf2);
    fP[2] = pf1;
    fP[3] = pf2;
    fM[2] = pf1.GetMass();
    fM[3] = pf2.GetMass();
    Double_t betaf = Beta(fM[2]*fM[2]/fQ2Z,fM[3]*fM[3]/fQ2Z);
    if (betaf <= 0.) return 0.;
    
    // ZHH rest frame
    Double_t betax = Beta(fQ2HH/fQ2ZHH,fQ2Z/fQ2ZHH);
    if (betax <= 0.) return 0.;
    
    // beam
    Double_t eb     = TMath::Sqrt(fQ2ZHH)/2.;
    Double_t pb     = TMath::Sqrt((eb-kM_e)*(eb+kM_e));
    Double_t beta_e = pb/eb;
    fK[0].SetXYZT(0., 0., pb, eb);
    fK[1].SetXYZT(0., 0.,-pb, eb);
    
    // --------------------------------------------
    //  Calcuate differential cross section
    // --------------------------------------------
    Double_t s      = fQ2ZHH;
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
    Double_t  color = f3Ptr->GetColor();
    Complex_t amp   = FullAmplitude();
    if (GetPropagator()) amp   *= GetHiggsPropagator(fM[0]*fM[0]) *
			          GetHiggsPropagator(fM[1]*fM[1]);
    Double_t  amp2  = TMath::Power(abs(amp),2) * color;
    
    // -------------------
    //  Put them together
    // -------------------
    static const Int_t    kNbr  = 3;
    static const Double_t kFact = k2Pi/(TMath::Power(k4Pi,3*kNbr));
    
    Double_t identp = 1./2.;                         // identical particle factor
    Double_t dPhase = kFact * betax * betaf * betah; // phase space factor
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
  Complex_t LCMEZHH::FullAmplitude()
  {
    Double_t gamz   = fZBosonPtr->GetWidth();
    Double_t qf     = f3Ptr->GetCharge();
    Double_t t3f    = f3Ptr->GetISpin();
    Double_t glz    = -kGz*(t3f - qf*kSin2W);
    Double_t grz    = -kGz*(    - qf*kSin2W);
    
    static const Bool_t kIsIncoming = kTRUE;
    static const Bool_t kIsOutgoing = kFALSE;
    
    HELFermion em(fK[0], kM_e, fHelInitial[0], +1, kIsIncoming);
    HELFermion ep(fK[1], kM_e, fHelInitial[1], -1, kIsOutgoing);
    
    HELScalar  h1(fP[0]);
    HELScalar  h2(fP[1]);
    
    HELFermion f (fP[2], fM[2], fHelFinal [2], +1, kIsOutgoing);
    HELFermion fb(fP[3], fM[3], fHelFinal [3], -1, kIsIncoming);
    HELVector  zf(fb, f, glz, grz, kM_z, gamz);
    
    Complex_t amp = AmpEEtoZHH(em, ep, h1, h2, zf);
    
    return amp;
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  AmpEEtoZHH()
  // --------------------------
  Complex_t LCMEZHH::AmpEEtoZHH(const HELFermion &em,
				const HELFermion &ep,
				const HELScalar  &h1,
				const HELScalar  &h2,
				const HELVector  &zf)
  {
    Double_t  qe    = -1.;
    Double_t  t3e   = -1./2.;
    Double_t  glze  = -kGz*(t3e - qe*kSin2W);
    Double_t  grze  = -kGz*(    - qe*kSin2W);
    Double_t  gamz  = fZBosonPtr->GetWidth();
    
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
    // Higgs Production Amplitude
    //---------------------------
    HELVector zs(em, ep, glze, grze, kM_z, gamz);
    
    Double_t v     = 2.*kM_w/kGw;
    Double_t ghhh  = -TMath::Power(fMass,2)/v*3.; 
    Double_t gzzh  = kGz*kM_z;
    Double_t gzzhh = kGz*kGz/2.; 
    
    HELScalar hh(h1, h2, ghhh, fMass, 0.);
    HELVertex amp1(zs, zf, hh, gzzh);         // HHH self-coupling
    
    HELVertex amp2(zs, zf, h1, h2, gzzhh);    // ZZHH 4-point
    
    HELVector vz1(zf, h1, gzzh, kM_z, gamz);
    HELVertex amp3(zs, vz1, h2, gzzh);        // double H-strahlung
    
    HELVector vz2(zf, h2, gzzh, kM_z, gamz);
    HELVertex amp4(zs, vz2, h1, gzzh);        // double H-strahlung
    
    Complex_t amp  = amp1 + amp2 + amp3 + amp4;
    //   cerr << "Debug: AMP = " << amp << endl;
#endif /* end __PHASESPACE__ */
    
    return amp;
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  Initialization
  // --------------------------
  void LCMEZHH::Initialize()
  {
    // --------------------------------------------
    //  Initialize Z decay table
    // --------------------------------------------
    if (!fZBosonPtr) fZBosonPtr = new GENPDTZBoson();
    
    
    // Set Higgs Mass
    //  SetMass(125.);
    
    
    // Set Beam Polarisations
    //  SetBeamPol(0.,0.);
    
    //  cerr << "LCMEZHH initialized!" << endl;
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  Set Z decay mode
  // --------------------------
  void LCMEZHH::SetZDecayMode(Int_t iDecayMode)
  {
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
  void LCMEZHH::SetMomentumFinal(TLorentzVector vLortz[])
  {
    fLortzZf1 = vLortz[0];
    fLortzZf2 = vLortz[1];
    fLortzH1  = vLortz[2];
    fLortzH2  = vLortz[3];
    
    // --------------------------
    //  Calculate variables in rest frame based on input kinematics
    // --------------------------
    // just providing information for verification, not actually being used by calculation of ME
    // Z rest frame
    TLorentzVector fLortzZ = fLortzZf1 + fLortzZf2;
    GetVariablesInRestFrame(fLortzZf1,fLortzZf2,fQ2Z,fCosThetaF,fPhiF);
    
    // HH rest frame
    TLorentzVector fLortzHH = fLortzH1 + fLortzH2;
    GetVariablesInRestFrame(fLortzH1,fLortzH2,fQ2HH,fCosThetaH,fPhiH);
    
    // ZHH rest frame
    //    TLorentzVector fLortzZHH = fLortzZ + fLortzHH;
    GetVariablesInRestFrame(fLortzZ,fLortzHH,fQ2ZHH,fCosTheta,fPhi);    
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  SetHelicities
  // --------------------------
  void LCMEZHH::SetHelicities(Int_t vHel[])
  {
    Int_t iHelI = vHel[0];
    Int_t iHelF = vHel[1];
    
    static const Int_t kNi = 2;
    static const Int_t kIHelComb[kNi][2] = {{-1, +1},
					    {+1, -1}};
    static const Int_t kNf = 2;
    static const Int_t kFHelComb[kNf][4] = {{0, 0, -1, +1},
					    {0, 0, +1, -1}};
    
    Int_t iJCombI = 0, iJCombF = 0;
    if (iHelI == -1) {
      iJCombI = 0;
    }
    else if (iHelI == 1) {
      iJCombI = 1;
    }
    if (iHelF == -1) {
      iJCombF = 0;
    }
    else if (iHelF == 1) {
      iJCombF = 1;
    }
    
    fHelInitial[0] = kIHelComb[iJCombI][0];
    fHelInitial[1] = kIHelComb[iJCombI][1];
    fHelFinal  [0] = kFHelComb[iJCombF][0];
    fHelFinal  [1] = kFHelComb[iJCombF][1];
    fHelFinal  [2] = kFHelComb[iJCombF][2];
    fHelFinal  [3] = kFHelComb[iJCombF][3];
  }
}
