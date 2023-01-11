#ifndef LCMEZHH_H
#define LCMEZHH_H
//*****************************************************************************
//* =====================
//*  LCMEZHH
//* =====================
//*  
//* (Description)
//*    e+e- -> ZHH Matrix Element
//*
//* (Update Record)
//*    2007/03/29  K.Fujii	Original version.
//*    2013/10/03  J.Tian       modified for MEM in Marlin
//*****************************************************************************

#include "TNamed.h"
#include "TLorentzVector.h"

#include "physsim/HELLib.h"
#include "physsim/GENLib.h"

#include "physsim/LCMEBase.h"
//_______________________________________________________________________
// =====================
//  class LCMEZHH
// =====================
//-----------------------------------------------------------------------
namespace lcme{
  class LCMEZHH : public LCMEBase {
  public:
    // --------------------------------------------------------------------
    //  Member Functions
    // --------------------------------------------------------------------
    // ----------------------
    //  C-tor and D-tor
    // ----------------------
    LCMEZHH(const char *name, const char *title,
	    Double_t massHiggs = 125.,
	    Double_t polE = 0.,
	    Double_t polP = 0.);
    virtual ~LCMEZHH();
    
    // ----------------------
    //  Getters and Setters
    // ----------------------
    Double_t GetMass     ()           const { return fMass;      }
    Double_t GetQ2ZHH    ()           const { return fQ2ZHH;     }
    Double_t GetCosTheta ()           const { return fCosTheta;  }
    Double_t GetPhi      ()           const { return fPhi;       }
    Double_t GetQ2Z      ()           const { return fQ2Z;       }
    Double_t GetCosThetaF()           const { return fCosThetaF; }
    Double_t GetPhiF     ()           const { return fPhiF;      }
    Double_t GetQ2HH     ()           const { return fQ2HH;      }
    Double_t GetCosThetaH()           const { return fCosThetaH; }
    Double_t GetPhiH     ()           const { return fPhiH;      }
    Bool_t   GetPropagator ()         const { return fPropagator;  }
    
    // ----------------------
    //   Base class methods
    // ----------------------
    Double_t GetMatrixElement2();                 // matrix element squared with weighted helicities
    Double_t GetMatrixElement2(Int_t vHel[]);     // matrix element squared with specified helicities
    void     SetMass     (Double_t m      ) { fMass      = m;    }
    void     SetMomentumFinal(TLorentzVector vLortz[]); // set four-momenta of final states
    void     SetZDecayMode(Int_t iDecayMode = 1);             // set Z decay mode
    void     SetPropagator(Bool_t i) {fPropagator = i;};
    
    // ----------------------
    //   Utility methods
    // ----------------------
  private:
    void Initialize();        // Bases initialization
    Double_t  DSigmaDX     ();
    Complex_t FullAmplitude();
    Complex_t AmpEEtoZHH   (const HELFermion &em,
			    const HELFermion &ep,
			    const HELScalar  &h1,
			    const HELScalar  &h2,
			    const HELVector  &zf);
    void     SetHelicities(Int_t vHel[]);
    
  private:
    // --------------------------------------------------------------------
    //  Data Members
    // --------------------------------------------------------------------
    //  Lorentz Vector of final state particle
    TLorentzVector fLortzZf1,fLortzZf2; // two fermions from Z
    TLorentzVector fLortzH1,fLortzH2; // two H
    
    // ----------------
    //  Z decay mode
    // ----------------
    Int_t    fZDecayMode;           // Z decay mode;
    GENDecayMode  *fZModePtr;       // point to Z decay mode 
    GENPDTEntry   *f3Ptr;           // point to 1st Z daughter
    GENPDTEntry   *f4Ptr;           // point to 2nd Z daughter
    
    // ----------------
    //  Particle Data
    // ----------------
    GENPDTZBoson *fZBosonPtr;       //! PD table entry of "Z"
    
    // ----------------
    //  Event info
    // ----------------
    Int_t          fHelInitial[2];  // initial state helicities
    Int_t          fHelFinal  [4];  // final   state helicities
    ANL4DVector    fK[2];           // [0,1] = [e-, e+]
    ANL4DVector    fP[4];           // [0,1,2,3] = [h1, h2, fz1, fz2]
    Double_t       fM[4];           // [0,1,2,4] = [mh, mh, m3 , m4 ]
    
    Double_t       fMass;           // m_h    : mass  of H
    Double_t       fQ2ZHH;          // q^2 of ZHH system
    Double_t       fQ2Z;            // q^2 of final state Z
    Double_t       fQ2HH;           // q^2 of final state HH
    Double_t       fCosTheta;       // cos(theta_x) in cm  frame
    Double_t       fPhi;            // phi_x        in cm  frame
    Double_t       fCosThetaF;      // cos(theta_f) in Z   frame
    Double_t       fPhiF;           // phi_f        in Z   frame
    Double_t       fCosThetaH;      // cos(theta_H) in HH  frame
    Double_t       fPhiH;           // phi_H        in HH  frame
    
    Bool_t         fPropagator;     // including B-W of Higgs
    
    ClassDef(LCMEZHH, 1) // Matrix Element for e+e- -> ZHH process
  };
}
#endif
