//*****************************************************************************
//* =====================
//*  LCMEZH
//* =====================
//*  
//* (Description)
//*    e+e- --> ZH Matrix Element
//*
//* (Update Record)
//*    2007/03/29  K.Fujii	Original version in Physsim as generator.
//*    2014/02/17  J.Tian	Modified to calculate Matrix Element only.
//*
//*****************************************************************************

#include "physsim/LCMEZH.h"

#include <sstream>
#include <iomanip>
#ifdef __PHASESPACE__
#define __NODECAY__
#endif

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "physsim/GENNumCon.h"

ClassImp(lcme::LCMEZH)

//-----------------------------------------------------------------------------
// ==============================
//  class LCMEZH
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
namespace lcme{
  LCMEZH::LCMEZH(const char *name, const char *title,
		 Double_t massHiggs,
		 Double_t polE,
		 Double_t polP)
  :LCMEBase(name,title,polE,polP),
   fZModePtr  (0),
   f3Ptr      (0),
   f4Ptr      (0),
   fZBosonPtr (0),
   iAnomalous (kFALSE),
   fA1 (1),
   fA2 (0),
   fA3 (0),
   fPropagator(0)
  {
    //  Constructor of bases.  Default parameter should be initialized here
    //
    cout << "Init LCMEZH " << endl;
    SetMass(massHiggs);
    //  SetBeamPol(polE,polP);
    Initialize();
  }
  // --------------------------
  //  D-tor
  // --------------------------
  LCMEZH::~LCMEZH()
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
  Double_t LCMEZH::GetMatrixElement2()
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
    
    Double_t weightElectron = (1.-GetBeamPolE())/2.;
    Double_t weightPositron = (1.+GetBeamPolP())/2.;
    
    Double_t sigma = 0.;
    sigma += (sigmaLL+sigmaLR)*weightElectron*weightPositron;
    sigma += (sigmaRL+sigmaRR)*(1.-weightElectron)*(1.-weightPositron);
    return (sigma);
  }
  Double_t LCMEZH::GetMatrixElement2(Int_t vHel[])
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
  Double_t LCMEZH::DSigmaDX()
  {
    // --------------------------------------------
    //  Phase space (&& calculate BetaBar)
    // --------------------------------------------
#if 0
    // H rest frame
    ANL4DVector ph(fLortzH);
    fP[0] = ph;
    fM[0] = ph.GetMass();
    
    // Z rest frame
    ANL4DVector pf1(fLortzZf1);
    ANL4DVector pf2(fLortzZf2);
    fP[1] = pf1;
    fP[2] = pf2;
    fM[1] = pf1.GetMass();
    fM[2] = pf2.GetMass();
#endif
    Double_t betaf = Beta(fM[1]*fM[1]/fQ2Z,fM[2]*fM[2]/fQ2Z);
    if (betaf <= 0.) return 0.;
    
    // ZH rest frame
    Double_t betax = Beta(fM[0]*fM[0]/fQ2ZH,fQ2Z/fQ2ZH);
    if (betax <= 0.) return 0.;
    
    // beam
    Double_t eb     = TMath::Sqrt(fQ2ZH)/2.;
    Double_t pb     = TMath::Sqrt((eb-kM_e)*(eb+kM_e));
    Double_t beta_e = pb/eb;
#if 0
    fK[0].SetXYZT(0., 0., pb, eb);
    fK[1].SetXYZT(0., 0.,-pb, eb);
#endif
    
    // --------------------------------------------
    //  Calcuate differential cross section
    // --------------------------------------------
    Double_t s      = fQ2ZH;
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
    ANL4DVector qz = fP[1] + fP[2];
    cerr << " qz = (" 
	 << qz.E () << ","
	 << qz.Px() << ","
	 << qz.Py() << ","
	 << qz.Pz() << ")" << endl;
    ANL4DVector pcm = fP[0] + qz;
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
    if (GetPropagator()) amp   *= GetHiggsPropagator(fM[0]*fM[0]);
    Double_t  amp2  = TMath::Power(abs(amp),2) * color;
    
    // -------------------
    //  Put them together
    // -------------------
    static const Int_t    kNbr  = 2;
    static const Double_t kFact = k2Pi/(TMath::Power(k4Pi,3*kNbr));
    
    Double_t identp = 1.;                            // identical particle factor
    Double_t dPhase = kFact * betax * betaf;         // phase space factor
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
  Complex_t LCMEZH::FullAmplitude()
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
    
    HELScalar  h(fP[0]);
    
    HELFermion f (fP[1], fM[1], fHelFinal [0], +1, kIsOutgoing);
    HELFermion fb(fP[2], fM[2], fHelFinal [1], -1, kIsIncoming);
    HELVector  zf(fb, f, glz, grz, kM_z, gamz);
    
    Complex_t amp = AmpEEtoZH(em, ep, h, zf);
    
    return amp;
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  AmpEEtoZH()
  // --------------------------
  Complex_t LCMEZH::AmpEEtoZH(const HELFermion &em,
			      const HELFermion &ep,
			      const HELScalar  &h,
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
    Double_t mz     = fZBosonPtr->GetMass();
    Double_t gzzh  = kGz*kM_z;
    Double_t  vev    = 2.*mz/kGz;
    Double_t g1 = gzzh * fA1;
    Double_t g2 = fA2/vev;
    Double_t g3 = fA3/vev;
    //    HELVertex amp(zs, zf, h, gzzh);         
    HELVertex amp = iAnomalous? HELVertex(zs, zf, h, g1, g2, g3) : HELVertex(zs, zf, h, gzzh);         
#endif /* end __PHASESPACE__ */
    
    return amp;
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  Initialization
  // --------------------------
  void LCMEZH::Initialize()
  {
    // --------------------------------------------
    //  Initialize Z decay table
    // --------------------------------------------
    if (!fZBosonPtr) fZBosonPtr = new GENPDTZBoson();
    
    // Set Higgs Mass
    //  SetMass(125.);
    
    
    // Set Beam Polarisations
    //  SetBeamPol(0.,0.);
    
    //  cerr << "LCMEZH initialized!" << endl;
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  Set Z decay mode
  // --------------------------
  void LCMEZH::SetZDecayMode(Int_t iDecayMode)
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
  void LCMEZH::SetMomentumFinal(TLorentzVector vLortz[])
  {
    fLortzZf1 = vLortz[0];
    fLortzZf2 = vLortz[1];
    fLortzH   = vLortz[2];
    
    // --------------------------
    //  Calculate variables in rest frame based on input kinematics
    // --------------------------
    // just providing information for verification, not actually being used by calculation of ME
    // Z rest frame
    TLorentzVector fLortzZ = fLortzZf1 + fLortzZf2;
    GetVariablesInRestFrame(fLortzZf1,fLortzZf2,fQ2Z,fCosThetaF,fPhiF);
    
    // ZH rest frame
    //    TLorentzVector fLortzZH = fLortzZ + fLortzH;
    GetVariablesInRestFrame(fLortzZ,fLortzH,fQ2ZH,fCosTheta,fPhi);    

    // feed to internal variables
    // H rest frame
    ANL4DVector ph(fLortzH);
    fP[0] = ph;
    fM[0] = ph.GetMass();
    // Z rest frame
    ANL4DVector pf1(fLortzZf1);
    ANL4DVector pf2(fLortzZf2);
    fP[1] = pf1;
    fP[2] = pf2;
    fM[1] = pf1.GetMass();
    fM[2] = pf2.GetMass();
    // initial beam
    Double_t eb     = TMath::Sqrt(fQ2ZH)/2.;
    Double_t pb     = TMath::Sqrt((eb-kM_e)*(eb+kM_e));
    fK[0].SetXYZT(0., 0., pb, eb);
    fK[1].SetXYZT(0., 0.,-pb, eb);

  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  SetHelicities
  // --------------------------
  void LCMEZH::SetHelicities(Int_t vHel[])
  {
    Int_t iHelI = vHel[0];
    Int_t iHelF = vHel[1];
    
    static const Int_t kNi = 2;
    static const Int_t kIHelComb[kNi][2] = {{-1, +1},
					    {+1, -1}};
    static const Int_t kNf = 2;
    static const Int_t kFHelComb[kNf][2] = {{-1, +1},
					    {+1, -1}};
    
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
  }
}
