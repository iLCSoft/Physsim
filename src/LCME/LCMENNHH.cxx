//*****************************************************************************
//* =====================
//*  LCMENNHH
//* =====================
//*  
//* (Description)
//*    e+e- --> NNHH Matrix Element
//*
//* (Update Record)
//*    2007/03/29  K.Fujii	Original version in Physsim as generator.
//*    2014/02/17  J.Tian	Modified to calculate Matrix Element only.
//*
//*****************************************************************************

#include "LCMENNHH.h"

#include <sstream>
#include <iomanip>
#ifdef __PHASESPACE__
#define __NODECAY__
#endif

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
#include "GENNumCon.h"

ClassImp(lcme::LCMENNHH)

//-----------------------------------------------------------------------------
// ==============================
//  class LCMENNHH
// ==============================
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
namespace lcme{
  LCMENNHH::LCMENNHH(const char *name, const char *title,
		   Double_t massHiggs,
		   Double_t polE,
		   Double_t polP)
 :LCMEBase(name,title,polE,polP)
  {
    //  Constructor of bases.  Default parameter should be initialized here
    //
    cout << "Init LCMENNHH " << endl;
    SetMass(massHiggs);
    //  SetBeamPol(polE,polP);
    Initialize();
  }
  // --------------------------
  //  D-tor
  // --------------------------
  LCMENNHH::~LCMENNHH()
  {
    delete fWBosonPtr;
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  Get Matrix Element Squared
  // --------------------------
  Double_t LCMENNHH::GetMatrixElement2()
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
  Double_t LCMENNHH::GetMatrixElement2(Int_t vHel[])
  {
    // with initial and final helicities combinations specified
    if (vHel[1] != -1) return 0.;
    SetHelicities(vHel);
    Double_t sigma = DSigmaDX();
    return (sigma);
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  DSigmaDX
  // --------------------------
  Double_t LCMENNHH::DSigmaDX()
  {
    // --------------------------------------------
    //  Phase space (&& calculate BetaBar)
    // --------------------------------------------
    // HH rest frame
    ANL4DVector ph1(fLortzH1);
    ANL4DVector ph2(fLortzH2);
    fP[2] = ph1;
    fM[2] = ph1.GetMass();
    fP[3] = ph2;
    fM[3] = ph2.GetMass();
    Double_t betah = Beta(fM[2]*fM[2]/fQ2HH,fM[3]*fM[3]/fQ2HH);
    if (betah <= 0.) return 0.;
    
    // NN rest frame
    ANL4DVector pn1(fLortzN1);
    ANL4DVector pn2(fLortzN2);
    fP[0] = pn1;
    fP[1] = pn2;
    fM[0] = 0.;
    fM[1] = 0.;
    Double_t betaf = 1.;
    if (betaf <= 0.) return 0.;
    
    // NNHH rest frame
    Double_t betax = Beta(fQ2HH/fQ2NNHH,fQ2NN/fQ2NNHH);
    if (betax <= 0.) return 0.;
    
    // beam
    Double_t eb     = TMath::Sqrt(fQ2NNHH)/2.;
    Double_t pb     = TMath::Sqrt((eb-kM_e)*(eb+kM_e));
    Double_t beta_e = pb/eb;
    fK[0].SetXYZT(0., 0., pb, eb);
    fK[1].SetXYZT(0., 0.,-pb, eb);
    
    // --------------------------------------------
    //  Calcuate differential cross section
    // --------------------------------------------
    Double_t s      = fQ2NNHH;
    // -------------------
    //  Amplitude squared
    // -------------------
    Complex_t amp   = FullAmplitude();
    Double_t  amp2  = TMath::Power(abs(amp),2);
    
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
  Complex_t LCMENNHH::FullAmplitude()
  {
    static const Bool_t kIsIncoming = kTRUE;
    static const Bool_t kIsOutgoing = kFALSE;
    
    HELFermion em(fK[0], kM_e, fHelInitial[0], +1, kIsIncoming);
    HELFermion ep(fK[1], kM_e, fHelInitial[1], -1, kIsOutgoing);
    
    HELScalar  h1(fP[2]);
    HELScalar  h2(fP[3]);
    
    HELFermion ne (fP[0], fM[0], fHelFinal [0], +1, kIsOutgoing);
    HELFermion neb(fP[1], fM[1], fHelFinal [1], -1, kIsIncoming);
    
   Complex_t amp = AmpEEtoNNHH(em, ep, ne, neb, h1, h2);
    
    return amp;
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  AmpEEtoNNHH()
  // --------------------------
  Complex_t LCMENNHH::AmpEEtoNNHH(const HELFermion &em,
				  const HELFermion &ep,
				  const HELFermion &ne,
				  const HELFermion &neb,
				  const HELScalar  &h1,
				  const HELScalar  &h2)
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
    Double_t mw     = fWBosonPtr->GetMass();
    Double_t gamw   = fWBosonPtr->GetWidth();
    Double_t glwf   = -kGw*kSqh;
    Double_t grwf   = 0.;
    
    HELVector w1(em , ne, glwf, grwf, mw, gamw);
    HELVector w2(neb, ep, glwf, grwf, mw, gamw);
    
    Double_t mh     = fMass;
    //    Double_t gamh   = fWidth;
    Double_t gamh   = 0.;
    Double_t gwwh   = kGw*mw;
    Double_t gwwhh  = kGw*kGw/2.;
    Double_t v      = 2*mw/kGw;
    Double_t ghhh   = -3.*mh*mh/v;
    
    // (1) self-coupling
    HELScalar hh(h1, h2, ghhh, mh, gamh);
    Complex_t amp1 = HELVertex(w1, w2, hh, gwwh);
    
    // (2) 4-point
    Complex_t amp2 = HELVertex(w1, w2, h1, h2, gwwhh);
    
    // (3) h1 and h2 from W
    HELVector w1h1(w1, h1, gwwh, mw, gamw);
    Complex_t amp3 = HELVertex(w1h1, w2, h2, gwwh);
    HELVector w1h2(w1, h2, gwwh, mw, gamw);
    amp3 += HELVertex(w1h2, w2, h1, gwwh);
    Complex_t amp = amp1 + amp2 + amp3;
#endif /* end __PHASESPACE__ */
    
    return amp;
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  Initialization
  // --------------------------
  void LCMENNHH::Initialize()
  {
    // --------------------------------------------
    //  Initialize W decay table
    // --------------------------------------------
    if (!fWBosonPtr) fWBosonPtr = new GENPDTWBoson();
    fWBosonPtr->DebugPrint();
    
    
    // Set Higgs Mass
    //  SetMass(125.);
    
    
    // Set Beam Polarisations
    //  SetBeamPol(0.,0.);
    
    //  cerr << "LCMENNHH initialized!" << endl;
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  Set four-momentum of final states
  // --------------------------
  void LCMENNHH::SetMomentumFinal(TLorentzVector vLortz[])
  {
    fLortzN1  = vLortz[0];
    fLortzN2  = vLortz[1];
    fLortzH1  = vLortz[2];
    fLortzH2  = vLortz[3];
    
    // --------------------------
    //  Calculate variables in rest frame based on input kinematics
    // --------------------------
    // just providing information for verification, not actually being used by calculation of ME
    // NN rest frame
    TLorentzVector fLortzNN = fLortzN1 + fLortzN2;
    GetVariablesInRestFrame(fLortzN1,fLortzN2,fQ2NN,fCosThetaF,fPhiF);
    
    // HH rest frame
    TLorentzVector fLortzHH = fLortzH1 + fLortzH2;
    GetVariablesInRestFrame(fLortzH1,fLortzH2,fQ2HH,fCosThetaH,fPhiH);
    
    // NNHH rest frame
    //    TLorentzVector fLortzNNHH = fLortzNN + fLortzHH;
    GetVariablesInRestFrame(fLortzNN,fLortzHH,fQ2NNHH,fCosTheta,fPhi);    
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  SetHelicities
  // --------------------------
  void LCMENNHH::SetHelicities(Int_t vHel[])
  {
    Int_t iHelI = vHel[0];
    
    static const Int_t kNi = 2;
    static const Int_t kIHelComb[kNi][2] = {{-1, +1},
					    {+1, -1}};
    static const Int_t kNf = 1;
    static const Int_t kFHelComb[kNf][4] = {{-1, +1, 0, 0}};
    
    Int_t iJCombI = 0, iJCombF = 0;
    if (iHelI == -1) {
      iJCombI = 0;
    }
    else if (iHelI == 1) {
      iJCombI = 1;
    }
    iJCombF = 0;
    
    fHelInitial[0] = kIHelComb[iJCombI][0];
    fHelInitial[1] = kIHelComb[iJCombI][1];
    fHelFinal  [0] = kFHelComb[iJCombF][0];
    fHelFinal  [1] = kFHelComb[iJCombF][1];
    fHelFinal  [2] = kFHelComb[iJCombF][2];
    fHelFinal  [3] = kFHelComb[iJCombF][3];
  }
}
