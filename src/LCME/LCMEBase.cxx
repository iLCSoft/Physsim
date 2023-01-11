//*****************************************************************************
//* =====================
//*  LCMEBase
//* =====================
//*  
//* (Description)
//*  Base class of e+e- --> XX Matrix Element
//*
//* (Update Record)
//*    2014/02/17  J.Tian	Original version.
//*****************************************************************************

#include "physsim/LCMEBase.h"

//#include <sstream>
//#include <iomanip>

//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
//#include "GENNumCon.h"

ClassImp(lcme::LCMEBase)

//-----------------------------------------------------------------------------
// ==============================
//  class LCMEBase
// ==============================
//_____________________________________________________________________________
namespace lcme{
  // --------------------------
  //  C-tor
  // --------------------------
  LCMEBase::LCMEBase(const char *name, const char *title,
		     Double_t polE,
		     Double_t polP)
  {
    //  Constructor of bases.  Default parameter should be initialized here
    //
    //  cout << "Init LCMEBase " << endl;
    SetBeamPol(polE,polP);
    SetMEType(1);
    SetPropagatorType(1);
    SetHiggsMass(125.);
    SetHiggsWidth(0.0043);
    //  Initialize();
  }
  // --------------------------
  LCMEBase::LCMEBase(const char *name, const char *title,
		     Double_t polE, Double_t polP, Int_t typeME)
  {
    //  Constructor of bases.  Default parameter should be initialized here
    //
    //  cout << "Init LCMEBase " << endl;
    SetBeamPol(polE,polP);
    SetMEType(typeME);
    SetPropagatorType(1);
    SetHiggsMass(125.);
    SetHiggsWidth(0.0043);
    //  Initialize();
  }
  // --------------------------
  //  D-tor
  // --------------------------
  LCMEBase::~LCMEBase()
  {
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  Beta2
  // --------------------------
  Double_t LCMEBase::Beta2(Double_t x1, Double_t x2)
  {
    return 1. - 2*(x1+x2) + (x1-x2)*(x1-x2);
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  Beta
  // --------------------------
  Double_t LCMEBase::Beta(Double_t x1, Double_t x2)
  {
    Double_t beta2 = Beta2(x1,x2);
    if (beta2 < 0.) {
      return 0.;
    }
    else {
      return TMath::Sqrt(beta2);
    }
  }
  
  //_____________________________________________________________________________
  // --------------------------
  //  SetBeamPolarisations
  // --------------------------
  void LCMEBase::SetBeamPol(Double_t polE, Double_t polP) 
  {
    fPolElectron = polE;
    fPolPositron = polP;
  }

  //_____________________________________________________________________________
  // --------------------------
  //  Higgs Propagator
  // --------------------------
  Complex_t LCMEBase::GetHiggsPropagator(Double_t q2)
  {
    Complex_t bw;
    if (GetPropagatorType() == 1) {
      std::cout << "Higgs Propagator (Breit-Wigner) is called: sqrt(Q2) = " << TMath::Sqrt(q2) << std::endl;
      bw = 1./Complex_t(q2-fMass*fMass, fMass*fWidth);
    }
    else if (GetPropagatorType() == 2) {
      std::cout << "Higgs Propagator (Effective Gaussian) is called: sqrt(Q2) = " << TMath::Sqrt(q2) << std::endl;
      bw = TMath::Exp(-(TMath::Sqrt(q2)-fMass)*(TMath::Sqrt(q2)-fMass)/fWidth/fWidth/4.);
    }
      return bw;
  }
  // --------------------------
  Complex_t LCMEBase::GetBosonPropagator(Double_t q2, Double_t mass, Double_t width)
  {
    Complex_t bw;
    if (GetPropagatorType() == 1) {
      std::cout << "Boson Propagator (Breit-Wigner) is called: sqrt(Q2) = " << TMath::Sqrt(q2) << "; Mass = " << mass << "; Width = " << width << std::endl;
      bw = 1./Complex_t(q2-mass*mass, mass*width);
    }
    else if (GetPropagatorType() == 2) {
      std::cout << "Boson Propagator (Effective Gaussian) is called: sqrt(Q2) = " << TMath::Sqrt(q2) << "; Mass = " << mass << "; Width = " << width << std::endl;
      bw = TMath::Exp(-(TMath::Sqrt(q2)-mass)*(TMath::Sqrt(q2)-mass)/width/width/4.);
    }
    return bw;
  }

  //_____________________________________________________________________________
  // --------------------------
  //  Calculate variables in rest frame based on input kinematics
  // --------------------------
  void LCMEBase::GetVariablesInRestFrame(TLorentzVector v1, TLorentzVector v2, 
					 Double_t &q2, Double_t &costheta, Double_t &phi)
  {
    TLorentzVector v0 = v1 + v2;
    q2 = v0.M2();
    if (v0.Vect().Mag()/v0.E() > 1.E-6) {
      // only do boost when beta > 1.E-6
      // redefine the axis
      TVector3    axiszOrig = TVector3(0., 0., 1.);
      TVector3    axis0z    = v0.Vect().Unit();
      TVector3    axis0x    = axis0z.Cross(axiszOrig).Unit();
      TVector3    axis0y    = axis0z.Cross(axis0x);
      // boost to rest frame
      TVector3    bst  = TVector3(0., 0., v0.Vect().Mag()/v0.E());
      TLorentzVector v10  = TLorentzVector(v1.Vect()*axis0x,v1.Vect()*axis0y,v1.Vect()*axis0z,v1.E());
      TLorentzVector v20  = TLorentzVector(v2.Vect()*axis0x,v2.Vect()*axis0y,v2.Vect()*axis0z,v2.E());
      v10.Boost(-bst);
      v20.Boost(-bst);
      costheta = v10.CosTheta();
      phi = v10.Phi();
    }
    else {
      costheta = v1.CosTheta();
      phi = v1.Phi();
    }
  }
}
