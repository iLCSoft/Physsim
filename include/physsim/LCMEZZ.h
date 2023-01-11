#ifndef LCMEZZ_H
#define LCMEZZ_H
//*****************************************************************************
//* =====================
//*  LCMEZZ
//* =====================
//*  
//* (Description)
//*    e+e- -> ZZ Matrix Element
//*
//* (Update Record)
//*    2012/03/30  K.Fujii	Original version.
//*    2014/10/28  J.Tian       modified for MEM in Marlin
//*****************************************************************************

#include "TNamed.h"
#include "TLorentzVector.h"

#include "physsim/HELLib.h"
#include "physsim/GENLib.h"

#include "physsim/LCMEBase.h"
//_______________________________________________________________________
// =====================
//  class LCMEZZ
// =====================
//-----------------------------------------------------------------------
namespace lcme{
  class LCMEZZ : public LCMEBase {
  public:
    // --------------------------------------------------------------------
    //  Member Functions
    // --------------------------------------------------------------------
    // ----------------------
    //  C-tor and D-tor
    // ----------------------
    LCMEZZ(const char *name, const char *title,
	    Double_t polE = 0.,
	    Double_t polP = 0.,
	    Int_t   iNoDecay = 0);
    virtual ~LCMEZZ();
    
    // ----------------------
    //  Getters and Setters
    // ----------------------
    Double_t GetQ2ZZ       ()           const { return fQ2ZZ;        }
    Double_t GetCosTheta   ()           const { return fCosTheta;    }
    Double_t GetPhi        ()           const { return fPhi;         }
    Double_t GetQ2Z1       ()           const { return fQ2Z1;        }
    Double_t GetCosThetaZ1F()           const { return fCosThetaZ1F; }
    Double_t GetPhiZ1F     ()           const { return fPhiZ1F;      }
    Double_t GetQ2Z2       ()           const { return fQ2Z2;        }
    Double_t GetCosThetaZ2F()           const { return fCosThetaZ2F; }
    Double_t GetPhiZ2F     ()           const { return fPhiZ2F;      }
    Bool_t   GetNoDecay    ()           const { return fNoDecay;     }
    Double_t GetPropagator ()           const { return fPropagator;  }
    Double_t GetZ1Width    ()           const { return fZ1Width;     }
    Double_t GetZ2Width    ()           const { return fZ2Width;     }
    
    // ----------------------
    //   Base class methods
    // ----------------------
    Double_t GetMatrixElement2();                 // matrix element squared with weighted helicities
    Double_t GetMatrixElement2(Int_t vHel[]);     // matrix element squared with specified helicities
    void     SetMomentumFinal(TLorentzVector vLortz[]); // set four-momenta of final states
    void     SetMomentumFinalNoDecay(TLorentzVector vLortz[]); // set four-momenta of final states
    void     SetMomentumFinalNoDecay2(TLorentzVector vLortz[]); // set four-momenta of final states
    void     SetZDecayMode(Int_t iDecayMode);             // set Z decay mode
    void     SetZDecayMode(Int_t iDecayMode1, Int_t iDecayMode2);             // set Z decay mode
    void     SetNoDecay(Int_t i) {fNoDecay = i;};
    void     SetPropagator(Bool_t i) {fPropagator = i;};
    void     SetZ1Width(Double_t w) { fZ1Width = w;};
    void     SetZ2Width(Double_t w) { fZ2Width = w;};
    
    // ----------------------
    //   Utility methods
    // ----------------------
  private:
    void Initialize();        // Bases initialization
    Double_t  DSigmaDX     ();
    Complex_t FullAmplitude();
    Complex_t AmpEEtoZZ (const HELFermion &em,
			 const HELFermion &ep,
			 const HELVector  &z1,
			 const HELVector  &z2);
    void     SetHelicities(Int_t vHel[]);
    void     SetHelicitiesNoDecay(Int_t vHel[]);
    void     SetHelicitiesNoDecay2(Int_t vHel[]);
    
  private:
    // --------------------------------------------------------------------
    //  Data Members
    // --------------------------------------------------------------------
    // ----------------
    //  Lorentz Vector of final state particle
    TLorentzVector fLortzZ1,fLortzZ2;   // Z1,Z2
    TLorentzVector fLortzZ1f1,fLortzZ1f2; // two fermions from Z1 decay
    TLorentzVector fLortzZ2f1,fLortzZ2f2; // two fermions from Z2 decay
    
    // ----------------
    //  Particle Data
    // ----------------
    Int_t    fZ1DecayMode;           // Z1 decay mode;
    Int_t    fZ2DecayMode;           // Z2 decay mode;
    GENPDTZBoson  *fZ1BosonPtr;     //! PD table entry of "Z1"
    GENDecayMode  *fZ1ModePtr;      // pointer to Z1 decay mode
    GENPDTEntry   *f1Ptr;           // point to 1st Z1 daughter
    GENPDTEntry   *f2Ptr;           // point to 2nd Z1 daughter
    GENPDTZBoson  *fZ2BosonPtr;     //! PD table entry of "Z2"
    GENDecayMode  *fZ2ModePtr;      // pointer to Z2 decay mode
    GENPDTEntry   *f3Ptr;           // point to 1st Z2 daughter
    GENPDTEntry   *f4Ptr;           // point to 2nd Z2 daughter
    
    // ----------------
    //  Event info
    // ----------------
    Int_t          fHelInitial[2];  // initial state helicities
    Int_t          fHelFinal  [4];  // final   state helicities
    ANL4DVector    fK[2];           // [0,1] = [e-, e+]
    ANL4DVector    fP[4];           // [0,1,2,3] = [e-, e+,  f, fb]
    Double_t       fM[4];           // [0,1,2,3] = [me, me, mf, mf]
    
    Double_t       fQ2ZZ;           // q^2 of ZZ system
    Double_t       fCosTheta;       // cos(theta_x) in cm  frame
    Double_t       fPhi;            // phi_x        in cm  frame
    Double_t       fQ2Z1;          // q^2 of final state Z1
    Double_t       fCosThetaZ1F;    // cos(theta_f) in Z1  frame
    Double_t       fPhiZ1F;         // phi_f        in Z1  frame
    Double_t       fQ2Z2;          // q^2 of final state Z2
    Double_t       fCosThetaZ2F;    // cos(theta_f) in Z2  frame
    Double_t       fPhiZ2F;         // phi_f        in Z2  frame

    Int_t          fNoDecay;        // 0: both decay; 1: one Z decay; 2: both not decay
    Bool_t         fPropagator;     // including B-W if Z not decay
    Double_t       fZ1Width;        // width of Z (set according to measurement)
    Double_t       fZ2Width;        // width of Z (set according to measurement)
    
    ClassDef(LCMEZZ, 1) // Matrix Element for e+e- -> ZZ process
  };
}
#endif
