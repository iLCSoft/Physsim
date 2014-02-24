//*****************************************************************************
//* =====================
//*  LCMEEEH
//* =====================
//*  
//* (Description)
//*    e+e- --> EEH Matrix Element
//*
//* (Update Record)
//*    2007/03/29  K.Fujii	Original version in Physsim as generator.
//*    2014/02/18  J.Tian	Modified to calculate Matrix Element only.
//*
//*****************************************************************************

#include "LCMEEEH.h"

#include <sstream>
#include <iomanip>
#ifdef __PHASESPACE__
#define __NODECAY__
#endif

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "GENNumCon.h"

ClassImp(lcme::LCMEEEH)

//-----------------------------------------------------------------------------
// ==============================
//  class LCMEEEH
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
namespace lcme{
  LCMEEEH::LCMEEEH(const char *name, const char *title,
		   Double_t massHiggs,
		   Double_t polE,
		   Double_t polP)
  :LCMEBase(name,title,polE,polP)
  {
    //  Constructor of bases.  Default parameter should be initialized here
    //
    cout << "Init LCMEEEH " << endl;
    SetMass(massHiggs);
    //  SetBeamPol(polE,polP);
    Initialize();
  }
  // --------------------------
  //  D-tor
  // --------------------------
  LCMEEEH::~LCMEEEH()
  {
    delete fZBosonPtr;
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  Get Matrix Element Squared
  // --------------------------
  Double_t LCMEEEH::GetMatrixElement2()
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
	Int_t vHel[2] = {i,j};
	sigma += wE*wP*GetMatrixElement2(vHel);
      }
    }
    return (sigma);
  }
  Double_t LCMEEEH::GetMatrixElement2(Int_t vHel[])
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
  Double_t LCMEEEH::DSigmaDX()
  {
    // --------------------------------------------
    //  Phase space (&& calculate BetaBar)
    // --------------------------------------------
    // H rest frame
    ANL4DVector ph(fLortzH);
    fP[0] = ph;
    fM[0] = ph.GetMass();
    
    // EE rest frame
    ANL4DVector pe1(fLortzE1);
    ANL4DVector pe2(fLortzE2);
    fP[1] = pe1;
    fP[2] = pe2;
    fM[1] = kM_e;
    fM[2] = kM_e;
    Double_t betaf = Beta(fM[1]*fM[1]/fQ2EE,fM[2]*fM[2]/fQ2EE);
    if (betaf <= 0.) return 0.;
    
    // EEH rest frame
    Double_t betax = Beta(fM[0]*fM[0]/fQ2EEH,fQ2EE/fQ2EEH);
    if (betax <= 0.) return 0.;
    
    // beam
    Double_t eb     = TMath::Sqrt(fQ2EEH)/2.;
    Double_t pb     = TMath::Sqrt((eb-kM_e)*(eb+kM_e));
    Double_t beta_e = pb/eb;
    fK[0].SetXYZT(0., 0., pb, eb);
    fK[1].SetXYZT(0., 0.,-pb, eb);
    
    // --------------------------------------------
    //  Calcuate differential cross section
    // --------------------------------------------
    Double_t s      = fQ2EEH;

    // -------------------
    //  Amplitude squared
    // -------------------
    Complex_t amp   = FullAmplitude();
    Double_t  amp2  = TMath::Power(abs(amp),2);
    
    // -------------------
    //  Put them together
    // -------------------
    static const Int_t    kNbr  = 2;
    static const Double_t kFact = k2Pi/(TMath::Power(k4Pi,3*kNbr));
    
    Double_t identp = 1.;                            // identical particle factor
    Double_t dPhase = kFact * betax * betaf;         // phase space factor
    Double_t flux   = 1./(2.* s * beta_e);           // beam flux factor
    
    Double_t sigma  = identp * flux * amp2 * dPhase; // in [1/GeV^2]
    sigma *= kGeV2fb;                                // now in [fb]
    
    return sigma;
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  FullAmplitude()
  // --------------------------
  Complex_t LCMEEEH::FullAmplitude()
  {
    static const Bool_t kIsIncoming = kTRUE;
    static const Bool_t kIsOutgoing = kFALSE;
    
    HELFermion em(fK[0], kM_e, fHelInitial[0], +1, kIsIncoming);
    HELFermion ep(fK[1], kM_e, fHelInitial[1], -1, kIsOutgoing);
    
    HELScalar  hs(fP[0]);
    
    HELFermion ne (fP[1], fM[1], fHelFinal [1], +1, kIsOutgoing);
    HELFermion neb(fP[2], fM[2], fHelFinal [2], -1, kIsIncoming);
    
    Complex_t amp = AmpEEtoEEH(em, ep, ne, neb, hs);
    
    return amp;
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  AmpEEtoEEH()
  // --------------------------
  Complex_t LCMEEEH::AmpEEtoEEH(const HELFermion &em,
				 const HELFermion &ep,
				 const HELFermion &e,
				 const HELFermion &eb,
				 const HELScalar  &hs)
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
    // Higgs Production Amplitude
    //---------------------------
    Double_t mz     = fZBosonPtr->GetMass();
    Double_t gamz   = fZBosonPtr->GetWidth();
    Double_t  qe    = -1.;
    Double_t  t3e   = -1./2.;
    Double_t  glze  = -kGz*(t3e - qe*kSin2W);
    Double_t  grze  = -kGz*(    - qe*kSin2W);
    
    HELVector z1(em, e , glze, grze, mz, gamz);
    HELVector z2(eb, ep, glze, grze, mz, gamz);
    
    Double_t  gzzh   = kGz*mz;
    Complex_t amp = HELVertex(z1, z2, hs, gzzh);
#endif /* end __PHASESPACE__ */
    
    return amp;
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  Initialization
  // --------------------------
  void LCMEEEH::Initialize()
  {
    // --------------------------------------------
    //  Initialize Z decay table
    // --------------------------------------------
    if (! fZBosonPtr) fZBosonPtr = new GENPDTZBoson();
    fZBosonPtr->DebugPrint();
    
    // Set Higgs Mass
    //  SetMass(125.);
    
    
    // Set Beam Polarisations
    //  SetBeamPol(0.,0.);
    
    //  cerr << "LCMEEEH initialized!" << endl;
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  Set four-momentum of final states
  // --------------------------
  void LCMEEEH::SetMomentumFinal(TLorentzVector vLortz[])
  {
    fLortzE1  = vLortz[0];
    fLortzE2  = vLortz[1];
    fLortzH   = vLortz[2];
    
    // --------------------------
    //  Calculate variables in rest frame based on input kinematics
    // --------------------------
    // just providing information for verification, not actually being used by calculation of ME
    // EE rest frame
    TLorentzVector fLortzEE = fLortzE1 + fLortzE2;
    GetVariablesInRestFrame(fLortzE1,fLortzE2,fQ2EE,fCosThetaF,fPhiF);
    
    // EEH rest frame
    //    TLorentzVector fLortzEEH = fLortzEE + fLortzH;
    GetVariablesInRestFrame(fLortzEE,fLortzH,fQ2EEH,fCosTheta,fPhi);    
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  SetHelicities
  // --------------------------
  void LCMEEEH::SetHelicities(Int_t vHel[])
  {
    Int_t iHelI = vHel[0];
    Int_t iHelF = vHel[1];
    
    static const Int_t kNi = 4;
    static const Int_t kIHelComb[kNi][2] = {{-1, -1},
					    {-1, +1},
					    {+1, -1},
					    {+1, +1}};
    static const Int_t kNf = 4;
    static const Int_t kFHelComb[kNf][3] = {{0, -1, -1},
					    {0, -1, +1},
					    {0, +1, -1},
					    {0, +1, +1}};
    
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
    if (iHelF == -2) {
      iJCombF = 0;
    }
    else if (iHelF == -1) {
      iJCombF = 1;
    }
    else if (iHelF == 1) {
      iJCombF = 2;
    }
    else if (iHelF == 2) {
      iJCombF = 3;
    }
    
    fHelInitial[0] = kIHelComb[iJCombI][0];
    fHelInitial[1] = kIHelComb[iJCombI][1];
    fHelFinal  [0] = kFHelComb[iJCombF][0];
    fHelFinal  [1] = kFHelComb[iJCombF][1];
    fHelFinal  [2] = kFHelComb[iJCombF][2];
  }
}
