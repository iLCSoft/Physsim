#ifndef LCMENNHH_H
#define LCMENNHH_H
//*****************************************************************************
//* =====================
//*  LCMENNHH
//* =====================
//*  
//* (Description)
//*    e+e- -> NNHH Matrix Element
//*
//* (Update Record)
//*    2007/03/29  K.Fujii	Original version.
//*    2014/02/17  J.Tian       modified for MEM in Marlin
//*****************************************************************************

#include "TNamed.h"
#include "TLorentzVector.h"

#include "HELLib.h"
#include "GENLib.h"

#include "LCMEBase.h"
//_______________________________________________________________________
// =====================
//  class LCMENNHH
// =====================
//-----------------------------------------------------------------------
namespace lcme{
  class LCMENNHH : public LCMEBase {
  public:
    // --------------------------------------------------------------------
    //  Member Functions
    // --------------------------------------------------------------------
    // ----------------------
    //  C-tor and D-tor
    // ----------------------
    LCMENNHH(const char *name, const char *title,
	     Double_t massHiggs = 125.,
	     Double_t polE = 0.,
	     Double_t polP = 0.);
    virtual ~LCMENNHH();
    
    // ----------------------
    //  Getters and Setters
    // ----------------------
    Double_t GetMass     ()           const { return fMass;      }
    Double_t GetQ2NNHH   ()           const { return fQ2NNHH;     }
    Double_t GetCosTheta ()           const { return fCosTheta;  }
    Double_t GetPhi      ()           const { return fPhi;       }
    Double_t GetQ2NN     ()           const { return fQ2NN;       }
    Double_t GetCosThetaF()           const { return fCosThetaF; }
    Double_t GetPhiF     ()           const { return fPhiF;      }
    Double_t GetQ2HH     ()           const { return fQ2HH;      }
    Double_t GetCosThetaH()           const { return fCosThetaH; }
    Double_t GetPhiH     ()           const { return fPhiH;      }
    
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
    Complex_t AmpEEtoNNHH (const HELFermion &em,
                           const HELFermion &ep,
                           const HELFermion &ne,
                           const HELFermion &neb,
                           const HELScalar  &h1,
                           const HELScalar  &h2);
    void     SetHelicities(Int_t vHel[]);
    
  private:
    // --------------------------------------------------------------------
    //  Data Members
    // --------------------------------------------------------------------
    //  Lorentz Vector of final state particle
    TLorentzVector fLortzN1,fLortzN2; // two neutrinos
    TLorentzVector fLortzH1,fLortzH2; // two H
    
    // ----------------
    //  Particle Data
    // ----------------
    GENPDTWBoson *fWBosonPtr;       //! PD table entry of "W"
    
    // ----------------
    //  Event info
    // ----------------
    Int_t          fHelInitial[2];  // initial state helicities
    Int_t          fHelFinal  [4];  // final   state helicities
    ANL4DVector    fK[2];           // [0,1] = [e-, e+]
    ANL4DVector    fP[4];           // [0,1,2,3] = [ne, neb, h1, h2]
    Double_t       fM[4];           // [0,1,2,3] = [mn, mn , mh, mh]
    
    Double_t       fMass;           // m_h    : mass  of H
    Double_t       fQ2NNHH;         // q^2 of NNHH system
    Double_t       fQ2NN;           // q^2 of final state NN
    Double_t       fQ2HH;           // q^2 of final state HH
    Double_t       fCosTheta;       // cos(theta_x) in cm  frame
    Double_t       fPhi;            // phi_x        in cm  frame
    Double_t       fCosThetaF;      // cos(theta_N) in NN  frame
    Double_t       fPhiF;           // phi_N        in NN  frame
    Double_t       fCosThetaH;      // cos(theta_H) in HH  frame
    Double_t       fPhiH;           // phi_H        in HH  frame
    
    
    ClassDef(LCMENNHH, 1) // Matrix Element for e+e- -> NNHH process
  };
}
#endif
