#ifndef LCMEEEH_H
#define LCMEEEH_H
//*****************************************************************************
//* =====================
//*  LCMEEEH
//* =====================
//*  
//* (Description)
//*    e+e- -> EEH Matrix Element
//*
//* (Update Record)
//*    2007/03/29  K.Fujii	Original version.
//*    2014/02/18  J.Tian       modified for MEM in Marlin
//*****************************************************************************

#include "TNamed.h"
#include "TLorentzVector.h"

#include "HELLib.h"
#include "GENLib.h"

#include "LCMEBase.h"
//_______________________________________________________________________
// =====================
//  class LCMEEEH
// =====================
//-----------------------------------------------------------------------
namespace lcme{
  class LCMEEEH : public LCMEBase {
  public:
    // --------------------------------------------------------------------
    //  Member Functions
    // --------------------------------------------------------------------
    // ----------------------
    //  C-tor and D-tor
    // ----------------------
    LCMEEEH(const char *name, const char *title,
	    Double_t massHiggs = 125.,
	    Double_t polE = 0.,
	    Double_t polP = 0.);
    virtual ~LCMEEEH();
    
    // ----------------------
    //  Getters and Setters
    // ----------------------
    Double_t GetMass     ()           const { return fMass;      }
    Double_t GetQ2EEH    ()           const { return fQ2EEH;      }
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
    void     SetMass     (Double_t m      ) { fMass      = m;    }
    void     SetMomentumFinal(TLorentzVector vLortz[]); // set four-momenta of final states
    
    // ----------------------
    //   Utility methods
    // ----------------------
  private:
    void Initialize();        // Bases initialization
    Double_t  DSigmaDX     ();
    Complex_t FullAmplitude();
    Complex_t AmpEEtoEEH (const HELFermion &em,
                          const HELFermion &ep,
                          const HELFermion &e,
                          const HELFermion &eb,
                          const HELScalar  &hs);
    void     SetHelicities(Int_t vHel[]);
    
  private:
    // --------------------------------------------------------------------
    //  Data Members
    // --------------------------------------------------------------------
    // ----------------
    //  Lorentz Vector of final state particle
    TLorentzVector fLortzE1,fLortzE2; // two electrons
    TLorentzVector fLortzH; // Higgs
    
    // ----------------
    //  Particle Data
    // ----------------
    GENPDTZBoson *fZBosonPtr;       //! PD table entry of "Z"
    
    // ----------------
    //  Event info
    // ----------------
    Int_t          fHelInitial[2];  // initial state helicities
    Int_t          fHelFinal  [2];  // final   state helicities
    ANL4DVector    fK[2];           // [0,1] = [e-, e+]
    ANL4DVector    fP[3];           // [0,1,2] = [h, e-, e+]
    Double_t       fM[3];           // [0,1,2] = [mh,me , me ]
    
    Double_t       fMass;           // m_h    : mass  of H
    Double_t       fQ2EEH;          // q^2 of EEHH system
    Double_t       fQ2EE;           // q^2 of final state EE
    Double_t       fCosTheta;       // cos(theta_x) in cm  frame
    Double_t       fPhi;            // phi_x        in cm  frame
    Double_t       fCosThetaF;      // cos(theta_N) in EE  frame
    Double_t       fPhiF;           // phi_N        in EE  frame
    
    
    ClassDef(LCMEEEH, 1) // Matrix Element for e+e- -> EEH process
  };
}
#endif
