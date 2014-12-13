#ifndef LCMEBase_H
#define LCMEBase_H
//*****************************************************************************
//* =====================
//*  LCMEBase
//* =====================
//*  
//* (Description)
//*  Base class of e+e- -> XX Matrix Element
//*
//* (Update Record)
//*    2014/02/17  J.Tian	Original version.
//*****************************************************************************

#include "TNamed.h"
#include "TLorentzVector.h"
#include "HELLib.h"

//_______________________________________________________________________
// =====================
//  class LCMEBase
// =====================
//-----------------------------------------------------------------------
namespace lcme {
  class LCMEBase{
  public:
    // --------------------------------------------------------------------
    //  Member Functions
    // --------------------------------------------------------------------
    // ----------------------
    //  C-tor and D-tor
    // ----------------------
    LCMEBase(const char *name, const char *title,
	     Double_t polE = 0.,
	     Double_t polP = 0.);
    LCMEBase(const char *name, const char *title,
	     Double_t polE,
	     Double_t polP,
	     Int_t    typeME);
    virtual ~LCMEBase();
    
    // ----------------------
    //  Getters and Setters
    // ----------------------
    
    // ----------------------
    //   Base class methods
    // ----------------------
    virtual Double_t GetMatrixElement2() = 0;             // matrix element squared with weighted helicities
    virtual Double_t GetMatrixElement2(Int_t vHel[]) = 0; // matrix element squared with specified helicities
    virtual void     SetMomentumFinal(TLorentzVector vLortz[]) = 0;   // set four-momenta of final states
    void             SetBeamPol(Double_t polE, Double_t polP);    // set beam polarisations
    inline Double_t  GetBeamPolE() {return fPolElectron;};    // get electron beam polarisation
    inline Double_t  GetBeamPolP() {return fPolPositron;};    // get positron beam polarisation
    inline void      SetMEType(Int_t i) {fMEType = i;};   // set fMEType to return matrix element or different cross section
    inline Int_t     GetMEType() {return fMEType;};   // get fMEType
    inline Int_t     GetPropagatorType() {return fPropagatorType;};   // get fPropatatorType
    inline void      SetPropagatorType(Int_t i) {fPropagatorType = i;};   // set fPropatatorType

    // ----------------------
    //   Utility methods
    // ----------------------
    virtual void Initialize() = 0;        // Bases initialization
    virtual Double_t  DSigmaDX     () = 0;
    virtual Complex_t FullAmplitude() = 0;
    Double_t Beta2(Double_t x1, Double_t x2);
    Double_t Beta(Double_t x1, Double_t x2);
    void     GetVariablesInRestFrame(TLorentzVector v1, TLorentzVector v2, 
				     Double_t &q2, Double_t &costheta, Double_t &phi);
    virtual void  SetHelicities(Int_t vHel[]) = 0;
    inline void     SetHiggsMass(Double_t mass) {fMass = mass;};
    inline Double_t GetHiggsMass() {return fMass;};
    inline void     SetHiggsWidth(Double_t width) {fWidth = width;};
    inline Double_t GetHiggsWidth() {return fWidth;};
    Complex_t GetHiggsPropagator(Double_t q2);
    Complex_t GetBosonPropagator(Double_t q2, Double_t mass, Double_t width);
    
    // --------------------------------------------------------------------
    //  Data Members
    // --------------------------------------------------------------------
    // ----------------
    //  Beam polarisations
    Double_t fPolElectron, fPolPositron;
    
    // ----------------
    //  Higgs
    // ----------------
    Double_t       fMass;           // m_h    : mass  of H
    Double_t       fWidth;          // w_h    : total width of H
    
    Int_t fMEType; // 1: return differential cross section; 2: matrix element squared
    Int_t fPropagatorType; // 1: Relatistic Breit-Wigner; 2: Effective Gaussian with Detector Resolution

    ClassDef(LCMEBase, 1) // Base Matrix Element for e+e- -> XX process
  };
}
#endif
