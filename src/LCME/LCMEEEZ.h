#ifndef LCMEEEZ_H
#define LCMEEEZ_H
//*****************************************************************************
//* =====================
//*  LCMEEEZ
//* =====================
//*  
//* (Description)
//*    e+e- -> EEZ Matrix Element
//*
//* (Update Record)
//*    2012/03/30  K.Fujii	Original version.
//*    2014/04/08  J.Tian       modified for MEM in Marlin
//*****************************************************************************

#include "TNamed.h"
#include "TLorentzVector.h"

#include "HELLib.h"
#include "GENLib.h"

#include "LCMEBase.h"
//_______________________________________________________________________
// =====================
//  class LCMEEEZ
// =====================
//-----------------------------------------------------------------------
namespace lcme{
  class LCMEEEZ : public LCMEBase {
  public:
    // --------------------------------------------------------------------
    //  Member Functions
    // --------------------------------------------------------------------
    // ----------------------
    //  C-tor and D-tor
    // ----------------------
    LCMEEEZ(const char *name, const char *title,
	    Double_t polE = 0.,
	    Double_t polP = 0.,
	    Bool_t   iNoDecay = 0);
    virtual ~LCMEEEZ();
    
    // ----------------------
    //  Getters and Setters
    // ----------------------
    Double_t GetQ2EEZ    ()           const { return fQ2EEZ;      }
    Double_t GetCosTheta ()           const { return fCosTheta;  }
    Double_t GetPhi      ()           const { return fPhi;       }
    Double_t GetQ2EE     ()           const { return fQ2EE;       }
    Double_t GetCosThetaF()           const { return fCosThetaF; }
    Double_t GetPhiF     ()           const { return fPhiF;      }
    
    // ----------------------
    //   Base class methods
    // ----------------------
    Double_t GetMatrixElement2();                 // matrix element squared with weighted helicities
    Double_t GetMatrixElement2(Int_t vHel[]);     // matrix element squared with specified helicities
    void     SetMomentumFinal(TLorentzVector vLortz[]); // set four-momenta of final states
    void     SetMomentumFinalNoDecay(TLorentzVector vLortz[]); // set four-momenta of final states
    void     SetZDecayMode(Int_t iDecayMode = 1);             // set Z decay mode
    
    // ----------------------
    //   Utility methods
    // ----------------------
  private:
    void Initialize();        // Bases initialization
    Double_t  DSigmaDX     ();
    Complex_t FullAmplitude();
    Complex_t AmpEEtoEEZ (const HELFermion &em,
                          const HELFermion &ep,
                          const HELFermion &e,
                          const HELFermion &eb,
                          const HELVector  &zf);
    void     SetHelicities(Int_t vHel[]);
    void     SetHelicitiesNoDecay(Int_t vHel[]);
    
  private:
    // --------------------------------------------------------------------
    //  Data Members
    // --------------------------------------------------------------------
    // ----------------
    //  Lorentz Vector of final state particle
    TLorentzVector fLortzE1,fLortzE2;   // two electrons
    TLorentzVector fLortzZ;             // Z
    TLorentzVector fLortzZf1,fLortzZf2; // two fermions from Z decay
    
    // ----------------
    //  Particle Data
    // ----------------
    Int_t    fZDecayMode;           // Z decay mode;
    GENPDTZBoson *fZBosonPtr;       //! PD table entry of "Z"
    GENDecayMode  *fZModePtr;       // point to Z decay mode 
    GENPDTEntry   *f3Ptr;           // point to 1st Z daughter
    GENPDTEntry   *f4Ptr;           // point to 2nd Z daughter
    
    // ----------------
    //  Event info
    // ----------------
    Int_t          fHelInitial[2];  // initial state helicities
    Int_t          fHelFinal  [4];  // final   state helicities
    ANL4DVector    fK[2];           // [0,1] = [e-, e+]
    ANL4DVector    fP[4];           // [0,1,2,3] = [e-, e+,  f, fb]
    Double_t       fM[4];           // [0,1,2,3] = [me, me, mf, mf]
    
    Double_t       fQ2EEZ;          // q^2 of EEZ system
    Double_t       fQ2EE;           // q^2 of final state EE
    Double_t       fCosTheta;       // cos(theta_Z) in cm  frame
    Double_t       fPhi;            // phi_Z        in cm  frame
    Double_t       fCosThetaE;      // cos(theta_E) in EE  frame
    Double_t       fPhiE;           // phi_E        in EE  frame
    Double_t       fQ2Z;            // q^2 of final state Z
    Double_t       fCosThetaF;      // cos(theta_f) in Z   frame
    Double_t       fPhiF;           // phi_f        in Z   frame
    Double_t       fPhi1;           // phi_1 for e- -> e-
    Double_t       fPhi2;           // phi_2 for e+ -> e+
    Double_t       fSh1;            // sh1=sin(theta/2)
    Double_t       fCh1;            // ch1=cos(theta/2)
    Double_t       fSh2;            // sh2
    Double_t       fCh2;            // ch2
    Bool_t         fNoDecay;        // z decays or not
    
    ClassDef(LCMEEEZ, 1) // Matrix Element for e+e- -> EEZ process
  };
}
#endif
